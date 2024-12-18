// Creates the uniqueness and bucketing data for bundy.
package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"iter"
	"path/filepath"
	"regexp"
	"sort"
	"time"

	"github.com/fluhus/biostuff/formats/fasta"
	"github.com/fluhus/biostuff/formats/sam"
	"github.com/fluhus/bundy/bowtie"
	"github.com/fluhus/bundy/common"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/hashx"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
)

const (
	qualThresh  = 2
	readStep    = 4
	bucket2Size = 1000

	useFastBowtie = false // Experiment: use bowtie's fast mapping.
)

var (
	refFile  = flag.String("r", "", "Bowtie reference")
	inGlob   = flag.String("i", "", "Input file glob pattern")
	outFile  = flag.String("o", "", "Output file")
	readLen  = flag.Int("l", 0, "Read length")
	part     = flag.Int("p", 1, "Part number")
	nparts   = flag.Int("np", 1, "Total number of parts")
	nthreads = flag.Int("t", 1, "Number of threads")

	inFiles []string
)

func main() {
	common.Die(parseArgs())

	fmt.Println("Found", len(inFiles), "input files")
	fmt.Println("Read length:", *readLen)
	fmt.Println("Read step:", readStep)
	fmt.Println("Qual:", qualThresh)
	fmt.Printf("Part: %d/%d\n", *part, *nparts)

	fmt.Println("Starting")
	t := time.Now()
	fa := makeFasta()
	var sams iter.Seq2[*sam.SAM, error]
	if useFastBowtie {
		sams = bowtie.MapReader(fa, *refFile, *nthreads,
			"-F", fmt.Sprintf("%d,%d", *readLen, readStep), "--very-fast")
	} else {
		sams = bowtie.MapReader(fa, *refFile, *nthreads,
			"-F", fmt.Sprintf("%d,%d", *readLen, readStep))
	}
	common.Die(checkSam(sams))
	fmt.Println("Took", time.Since(t))
	fmt.Println("Done")
}

// Parses and checks arguments.
func parseArgs() error {
	flag.Parse()
	if *readLen <= 0 {
		return fmt.Errorf("bad read length (-l): %d", *readLen)
	}
	if inFiles, _ = filepath.Glob(*inGlob); len(inFiles) == 0 {
		return fmt.Errorf("no input files found (-i)")
	}
	if *outFile == "" {
		return fmt.Errorf("no output file selected (-o)")
	}
	if *refFile == "" {
		return fmt.Errorf("no reference selected (-r)")
	}
	return nil
}

// Generates a fasta subset stream from the input genomes.
func makeFasta() io.Reader {
	r, w := io.Pipe()
	go func() {
		for _, f := range inFiles {
			if err := makeFastaFile(f, w); err != nil {
				w.CloseWithError(err)
			}
		}
		w.Close()
	}()
	return r
}

// Writes fastq reads from the given fasta to the given writer.
func makeFastaFile(file string, w io.Writer) error {
	for fa, err := range fasta.File(file) {
		if err != nil {
			return err
		}
		if hashx.Bytes(fa.Name)%uint64(*nparts) != uint64(*part-1) {
			continue
		}
		txt, _ := fa.MarshalText()
		if _, err := w.Write(txt); err != nil {
			return err
		}
	}
	return nil
}

// Aggregates mapping results and creates the data for bundy.
func checkSam(sams iter.Seq2[*sam.SAM, error]) error {
	all := map[string]int{}
	ok := map[string]int{}
	okPos := snm.NewDefaultMap(func(s string) sets.Set[int] {
		return sets.Set[int]{}
	})

	pt := ptimer.NewMessage("loading reference")
	fre := regexp.MustCompile(`^(.*)_(\d+)$`)

	for sm, err := range sams {
		if err != nil {
			return err
		}
		if pt.N == 0 {
			pt.Done()
			pt = ptimer.NewMessage("{} reads")
			pt.Inc()
		}
		pt.Inc()

		var splt []string
		match := fre.FindStringSubmatch(sm.Qname)
		splt = []string{match[2], match[1]}
		rname := splt[1]
		all[rname]++
		if sm.Flag == sam.FlagUnmapped {
			continue
		}
		if sm.Mapq < qualThresh {
			continue
		}

		ok[rname]++
		okPosSet := okPos.Get(rname)
		okPosSet.Add(sm.Pos)
	}
	pt.Done()

	mulByReadStep(all)
	mulByReadStep(ok)

	f, err := aio.Create(*outFile)
	if err != nil {
		return err
	}
	j := json.NewEncoder(f)
	for k := range all {
		j.Encode(map[string]any{
			"name":    k,
			"all":     all[k],
			"ok":      ok[k],
			"buckets": posToBuckets(all[k], okPos.Get(k)),
		})
	}
	f.Close()
	fmt.Println("Wrote to:", *outFile)
	return nil
}

// Multiplies raw counts by read step to simulate real counts.
func mulByReadStep(m map[string]int) {
	for k := range m {
		m[k] *= readStep
	}
}

type bucketOKs struct {
	Buckets []int
	OK      []int
}

// Creates bucket positions from all the positions that were mapped to.
func posToBuckets(all int, okPos sets.Set[int]) *bucketOKs {
	nBuckets := max(1, gnum.Idiv(all, bucket2Size))
	var buckets []int
	for i := 1; i < nBuckets; i++ {
		buckets = append(buckets, gnum.Idiv(all*i, nBuckets))
	}
	ok := make([]int, nBuckets)
	for pos := range okPos {
		i := sort.SearchInts(buckets, pos)
		ok[i] += readStep
	}
	return &bucketOKs{buckets, ok}
}
