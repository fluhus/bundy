// Calculates relative abundance.
package main

import (
	"cmp"
	"flag"
	"fmt"
	"io"
	"iter"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"runtime/debug"
	"sort"
	"strings"

	"github.com/fluhus/biostuff/formats/fastq"
	"github.com/fluhus/biostuff/formats/sam"
	"github.com/fluhus/bundy/bowtie"
	"github.com/fluhus/bundy/common"
	"github.com/fluhus/gostuff/aio"
	"github.com/fluhus/gostuff/gnum"
	"github.com/fluhus/gostuff/jio"
	"github.com/fluhus/gostuff/ptimer"
	"github.com/fluhus/gostuff/sets"
	"github.com/fluhus/gostuff/snm"
	"golang.org/x/exp/maps"
)

const (
	qualThresh     = 30   // Quality threshold for first pass.
	qualThresh2    = 2    // Quality threshold for second pass.
	denseSumRatio  = 20   // Dense sum ratio for first pass.
	denseSumRatio2 = 20   // Dense sum ratio for second pass.
	minNZ          = 0.01 // Minimal coverage for first pass.
	minNZ2         = 0.66 // Minimal coverage for second pass.
	maxBinomialErr = 0.05 // Disqualify buckets with this binomial error.
)

// Debug stuff.
const (
	printConsts     = false
	printCumQuals   = false
	denseSumNonZero = false
	printSpecies    = ""
	printWhiteList  = false
	assertNZNonNeg  = true
	printRawCounts  = false
	printNReads     = true
)

var (
	inFile        = flag.String("i", "", "Input fastq file")
	inFile2       = flag.String("i2", "", "Second input fastq file for paired-end")
	outFile       = flag.String("o", "", "Output TSV file")
	refFile       = flag.String("r", "", "Bowtie reference")
	usedFile      = flag.String("u", "", "Print USED reads to this fastq")
	unusedFile    = flag.String("uu", "", "Print UNUSED reads to this fastq")
	usedSAMFile   = flag.String("us", "", "Print USED reads to this SAM")
	unusedSAMFile = flag.String("uus", "", "Print UNUSED reads to this SAM")
	oksGlob       = flag.String("x", "", "Bundyx data files glob")
	threads       = flag.Int("t", 1, "Number of bowtie2 threads")
	toJSON        = flag.Bool("j", false, "Output JSON instead of TSV")
	namePat       = flag.String("n", ".*",
		"Pattern by which to group contigs of the same species")

	ignoreLength = flag.Bool("ignlen", false, "Ignore genome lengths in normalization")
	fast         = flag.Bool("fast", false, "Quick run, loses some accuracy")
	interleaved  = flag.Bool("interleaved", false, "Input fasta has interleaved paired-end reads")
	diskMode     = flag.String("diskmode", "", "Write intermediate data to a `directory` rather than to RAM")

	speciesToPrint = createPrintSpeciesMap()
)

func main() {
	if printConsts {
		fmt.Fprintln(os.Stderr, "**********")
		fmt.Fprintln(os.Stderr, "qualThresh", qualThresh, qualThresh2)
		fmt.Fprintln(os.Stderr, "denseSumRatio", denseSumRatio, denseSumRatio2)
		fmt.Fprintln(os.Stderr, "printCumQuals", printCumQuals)
		fmt.Fprintln(os.Stderr, "denseSumNonZero", denseSumNonZero)
		fmt.Fprintln(os.Stderr, "printSpecies", printSpecies)
		fmt.Fprintln(os.Stderr, "minNZ", minNZ, minNZ2)
		fmt.Fprintln(os.Stderr, "printWhiteList", printWhiteList)
		fmt.Fprintln(os.Stderr, "maxBinomialErr", maxBinomialErr)
		fmt.Fprintln(os.Stderr, "assertNZNonNeg", assertNZNonNeg)
		fmt.Fprintln(os.Stderr, "**********")
		fmt.Fprintln(os.Stderr)
	}

	flag.Parse()
	debug.SetGCPercent(20)

	if *inFile == "" {
		common.Die(fmt.Errorf("no input file"))
	}
	if *outFile == "" {
		common.Die(fmt.Errorf("no output file"))
	}
	if _, err := os.Stat(*inFile); err != nil {
		common.Die(fmt.Errorf("unable to access input file: %w", err))
	}
	var err error
	nameRE, err = regexp.Compile(*namePat)
	common.Die(err)

	fmt.Fprintln(os.Stderr, "Running with:")
	fmt.Fprintln(os.Stderr, "\tRef:\t", shortenString(*refFile, 70))
	fmt.Fprintln(os.Stderr, "\tOKs:\t", shortenString(*oksGlob, 70))
	fmt.Fprintln(os.Stderr, "\tRegex:\t", nameRE)
	fmt.Fprintln(os.Stderr)

	fmt.Fprintln(os.Stderr, "Loading bundyx data")
	pt := ptimer.New()
	var entries map[string]*contigEntry
	entries, err = loadBuckets()
	pt.Done()
	common.Die(err)

	fmt.Fprintln(os.Stderr, "Mapping")
	all, unmapped, lowq := 0, 0, 0
	pt = ptimer.NewMessage("Loading reference")
	quals := map[int]int{}
	samw, err := samWriter()
	common.Die(err)

	var sams iter.Seq2[*sam.SAM, error]
	args := common.If(*fast, []string{"--very-fast"}, nil)
	if *inFile2 == "" {
		if *interleaved {
			sams = bowtie.MapInt(*inFile, *refFile, *threads, args...)
		} else {
			sams = bowtie.Map(*inFile, *refFile, *threads, args...)
		}
	} else {
		sams = bowtie.Map2(*inFile, *inFile2, *refFile, *threads, args...)
	}

	nreads := 0 // TODO(amit): Consider whether we need this.
	for sm, err := range sams {
		common.Die(err)
		if pt.N == 0 {
			pt.Done()
			pt = ptimer.NewMessage("{} reads processed")
			pt.Inc()
		}
		pt.Inc()

		txt, _ := sm.MarshalText()
		samw.Write(txt)

		all++
		if sm.Flag == sam.FlagUnmapped {
			unmapped++
			continue
		}
		quals[sm.Mapq]++
		if sm.Mapq < qualThresh {
			lowq++
			continue
		}

		nreads++
		entries[sm.Rname].addPos(sm.Pos)
	}
	pt.Done()
	if printNReads {
		fmt.Println("NReads:", nreads)
	}
	samw.Close()

	fmt.Fprintf(os.Stderr, "Mapped OK %v | Low quality %v | Unmapped %v\n",
		common.Percf(all-unmapped-lowq, all, 0),
		common.Percf(lowq, all, 0),
		common.Percf(unmapped, all, 0))

	var wl sets.Set[string]
	{ // TODO(amit): Make this a function?
		wl = sets.FromKeys(
			entriesToAbundances(entries, denseSumRatio, minNZ, 0))
		fmt.Fprintln(os.Stderr, "Found", len(wl), "candidate genomes")
		if printWhiteList {
			fmt.Fprintln(os.Stderr, wl)
		}
		for _, e := range entries {
			e.counts = nil
		}

		all, unmapped, lowq := 0, 0, 0
		pt = ptimer.NewMessage("{} reads processed")
		quals := map[int]int{}

		nreads := 0 // TODO(amit): Consider whether we need this.
		for sm, err := range samReader() {
			common.Die(err)
			pt.Inc()
			all++
			if sm.Flag == sam.FlagUnmapped {
				unmapped++
				continue
			}
			quals[sm.Mapq]++
			if sm.Mapq < qualThresh2 {
				lowq++
				continue
			}
			nreads++
			entries[sm.Rname].addPos(sm.Pos)
		}
		pt.Done()
		if printNReads {
			fmt.Println("NReads:", nreads)
		}
		fmt.Fprintf(os.Stderr, "Mapped OK %v | Low quality %v | Unmapped %v\n",
			common.Percf(all-unmapped-lowq, all, 0),
			common.Percf(lowq, all, 0),
			common.Percf(unmapped, all, 0))
	}

	abnd := entriesToAbundances(entries, denseSumRatio2, minNZ2, maxBinomialErr)
	if printRawCounts {
		abnd = entriesToRawCounts(entries)
	}
	abnd = snm.FilterMap(abnd, func(s string, f float64) bool {
		return wl.Has(s)
	})
	if !printRawCounts {
		toSum1(abnd)
	}
	fmt.Fprintln(os.Stderr, "Grouped to", len(abnd), "genomes")

	fmt.Fprintln(os.Stderr, "Saving")
	if *toJSON {
		common.Die(jio.Save(*outFile, abnd))
	} else {
		of, err := aio.Create(*outFile)
		common.Die(err)
		keys := snm.SortedFunc(maps.Keys(abnd), func(a, b string) int {
			return cmp.Compare(abnd[b], abnd[a])
		})
		for _, k := range keys {
			if printRawCounts {
				fmt.Fprintf(of, "%s\t%d\n", k, int(abnd[k]))
			} else {
				fmt.Fprintf(of, "%s\t%g\n", k, abnd[k])
			}
		}
		of.Close()
	}

	// Debug stats printing.
	if printCumQuals {
		cumquals := map[int]float64{}
		keys := snm.Sorted(maps.Keys(quals))
		sum := float64(gnum.Sum(maps.Values(quals)))
		cumquals[keys[0]] = float64(quals[keys[0]]) / sum
		for i, k := range keys[1:] {
			cumquals[k] = cumquals[keys[i]] + float64(quals[k])/sum
		}
		for _, k := range keys {
			fmt.Fprintf(os.Stderr, "%v:%.2f ", k, cumquals[k])
		}
		fmt.Fprintln(os.Stderr)
	}

	// Dump used/unused reads.
	if cmp.Or(*usedFile, *unusedFile, *usedSAMFile, *unusedSAMFile) != "" {
		fmt.Fprintln(os.Stderr, "Dumping used/unused reads")
		var uout, uuout, usout, uusout io.WriteCloser
		if *usedFile != "" {
			uout, err = aio.Create(*usedFile)
			common.Die(err)
		}
		if *unusedFile != "" {
			uuout, err = aio.Create(*unusedFile)
			common.Die(err)
		}
		if *usedSAMFile != "" {
			usout, err = aio.Create(*usedSAMFile)
			common.Die(err)
		}
		if *unusedSAMFile != "" {
			uusout, err = aio.Create(*unusedSAMFile)
			common.Die(err)
		}
		nused := 0
		pt = ptimer.NewMessage("{} reads processed")

		for sm, err := range samReader() {
			common.Die(err)
			pt.Inc()
			if !samMapped(sm) || sm.Mapq < qualThresh2 ||
				abnd[nameRE.FindString(sm.Rname)] == 0 { // Unused.
				if uuout != nil {
					common.Die(writeSamAsFastq(sm, uuout))
				}
				if uusout != nil {
					txt, _ := sm.MarshalText()
					_, err := uusout.Write(txt)
					common.Die(err)
				}
			} else { // Used.
				nused++
				if uout != nil {
					common.Die(writeSamAsFastq(sm, uout))
				}
				if usout != nil {
					txt, _ := sm.MarshalText()
					_, err := usout.Write(txt)
					common.Die(err)
				}
			}
		}
		pt.Done()
		fmt.Fprintf(os.Stderr, "Used %v of the reads\n", common.Percf(nused, all, 1))
		common.Die(closeAll(uout, uuout, usout, uusout))
	}

	fmt.Fprintln(os.Stderr, "Done")
}

var nameRE *regexp.Regexp

// Returns the sum of a, discarding some outliers.
func fDenseSum(a []float64, ratio int, nz float64) float64 {
	if assertNZNonNeg && nz < 0 { // Debug assert.
		panic(fmt.Sprintf("negative nz: %f", nz))
	}
	// No use for a window. For len=2 it will return the lower value.
	if len(a) <= 1 {
		return gnum.Sum(a)
	}

	// Actual dense sum.
	sort.Float64s(a)
	if nz != 0 {
		i := len(a) - 1 - iround(float64(len(a)-1)*nz)
		if a[i] == 0 {
			return 0 // Too many zeros.
		}
	}
	n := len(a)
	if denseSumNonZero {
		a = snm.FilterSlice(a, func(f float64) bool { return f > 0 })
		if len(a) <= 1 {
			return gnum.Sum(a) * float64(n)
		}
	}

	winlen := len(a)
	if ratio > 1 {
		winlen = gnum.Idiv(len(a)*(ratio-1), ratio)
		if winlen == len(a) {
			winlen--
		}
	}

	min, pos := math.Inf(1), 0
	for i := range a[winlen-1:] { // We have winlen+1 windows.
		diff := a[i+winlen-1] - a[i]
		if diff < min {
			min, pos = diff, i
		}
	}
	return gnum.Sum(a[pos:pos+winlen]) * float64(n) / float64(winlen)
}

// Shortens a string for display.
func shortenString(s string, n int) string {
	if len(s) <= n {
		return s
	}
	const filler = " ... "
	pre := (n - len(filler)) / 2
	suf := (n - len(filler) + 1) / 2
	return s[:pre] + filler + s[len(s)-suf:]
}

type contigEntry struct {
	ok      int        // Total OK mappings.
	all     int        // All positions.
	buckets *bucketOKs // Per-bucket information.
	counts  []int      // Mapping counts.
	sum     float64    // Dense sum.
}

type bucketOKs struct {
	pos []int // Boundries of buckets.
	ok  []int // OK mappings per bucket.
}

// Adds one count to the bucket at the given position.
func (e *contigEntry) addPos(pos int) {
	if e.counts == nil {
		e.counts = make([]int, len(e.buckets.ok))
	}
	bucket := sort.SearchInts(e.buckets.pos, pos)
	e.counts[bucket]++
}

// Loads and returns bundyx buckets.
func loadBuckets() (map[string]*contigEntry, error) {
	type entry struct {
		OK      int
		All     int
		Name    string
		Buckets struct {
			Buckets []int
			OK      []int
		}
	}

	result := map[string]*contigEntry{}
	files, _ := filepath.Glob(*oksGlob)
	if len(files) == 0 {
		return nil, fmt.Errorf("no OKs files found")
	}

	for _, file := range files {
		for e, err := range jio.Iter[entry](file) {
			if err != nil {
				return nil, err
			}
			result[e.Name] = &contigEntry{
				ok: e.OK, all: e.All,
				buckets: &bucketOKs{e.Buckets.Buckets, e.Buckets.OK},
			}
		}
	}
	return result, nil
}

// Returns a float as an integer, rounded to the nearest whole.
func iround(f float64) int {
	return int(math.Round(f))
}

// Populates the species debug map according to the constant.
func createPrintSpeciesMap() map[string]bool {
	if printSpecies == "" {
		return nil
	}
	m := map[string]bool{}
	for _, s := range strings.Split(printSpecies, ",") {
		m[s] = true
	}
	return m
}

// Returns a map from species to relative abundance.
func entriesToAbundances(m map[string]*contigEntry, ratio int, mz float64, maxBinom float64) map[string]float64 {
	abnd := map[string]float64{}
	aggEntries := snm.NewDefaultMap(func(s string) *contigEntry {
		return &contigEntry{counts: []int{0}}
	})
	aggBuckets := map[string][]float64{} // TODO(amit): See if I can get rid of this.
	aggOK := map[string]int{}            // Original OKs before subtracting in normCounts, for debugging.
	for s, e := range m {
		match := nameRE.FindString(s)
		agg := aggEntries.Get(match)
		agg.all += e.all
		agg.ok += e.ok
		aggOK[match] += len(e.buckets.ok)
		bucketSize := e.all / len(e.buckets.ok)
		normCounts := snm.Slice(len(e.buckets.ok), func(i int) float64 {
			if e.buckets.ok[i] < bucketSize/10 {
				agg.ok -= e.buckets.ok[i]
				agg.all -= bucketSize
				return math.NaN()
			}
			cnt := 0
			if len(e.counts) > 0 {
				cnt = e.counts[i]
			}
			return float64(cnt) * float64(bucketSize) / float64(e.buckets.ok[i])
		})

		normCounts = snm.FilterSlice(normCounts, func(f float64) bool {
			return !math.IsNaN(f) && !math.IsInf(f, 0)
		})
		aggBuckets[match] = append(aggBuckets[match], normCounts...)
	}
	filteredBinom := 0
	for s, e := range aggEntries.M {
		if printSpecies != "" && speciesToPrint[s] {
			toPrint := fmt.Sprintf("%.0f", aggBuckets[s])
			nzeros := len(snm.FilterSlice(aggBuckets[s], func(f float64) bool {
				return f == 0
			}))
			pzeros := float64(nzeros) / float64(len(aggBuckets[s])) * 100
			nzperc := fmt.Sprintf("%.0f%%", 100-pzeros)
			binerr := fmt.Sprintf("%.2f", binomialError(aggBuckets[s], aggOK[s]))
			fmt.Fprintln(os.Stderr, s, len(aggBuckets[s]), aggOK[s], nzperc, binerr, toPrint)
		}
		e.sum = fDenseSum(aggBuckets[s], ratio, mz)
		if e.sum > 0 && maxBinom > 0 && binomialError(aggBuckets[s], aggOK[s]) > maxBinom {
			e.sum = 0
			filteredBinom++
		}
		if e.sum == 0 || e.ok == 0 {
			continue
		}
		if !*ignoreLength {
			abnd[s] = e.sum / float64(e.all)
		} else {
			abnd[s] = e.sum * float64(e.all) / float64(e.ok)
		}
	}
	if maxBinom > 0 {
		fmt.Fprintln(os.Stderr, "Filtered binom:", filteredBinom)
	}
	toSum1(abnd)
	return abnd
}

// Returns a count map for the given entries.
func entriesToRawCounts(m map[string]*contigEntry) map[string]float64 {
	abnd := map[string]float64{}
	for s, e := range m {
		match := nameRE.FindString(s)
		abnd[match] += float64(gnum.Sum(e.counts))
	}
	return abnd
}

// Normalizes m's values to sum up to 1.
func toSum1(m map[string]float64) {
	sum := gnum.Sum(maps.Values(m))
	for k := range m {
		m[k] /= sum
	}
}

// Writes the given SAM line as a fastq entry.
func writeSamAsFastq(sm *sam.SAM, w io.Writer) error {
	fq := fastq.Fastq{
		Name:     []byte(sm.Qname),
		Sequence: []byte(sm.Seq),
		Quals:    []byte(sm.Qual),
	}
	txt, _ := fq.MarshalText()
	_, err := w.Write(txt)
	return err
}

// Returns the binomial STD for the given params.
func binomialError(abnd []float64, nn int) float64 {
	n := len(abnd)
	k := gnum.Sum(abnd)
	q := float64(nn-n) / float64(nn)
	return math.Sqrt(q / k)
}

// Returnes true if this SAM entry was mapped properly,
// depending on whether it's single or paried end.
func samMapped(sm *sam.SAM) bool {
	paired := *inFile2 != "" || *interleaved
	return (paired && sm.Flag&sam.FlagEach > 0) ||
		(!paired && sm.Flag&sam.FlagUnmapped == 0)
}

// Closes non-nil writers.
func closeAll(w ...io.WriteCloser) error {
	var err error
	for _, w := range w {
		if w != nil {
			err = w.Close()
		}
	}
	return err
}
