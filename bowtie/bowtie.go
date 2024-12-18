// Package bowtie provides functionality for running bowtie2.
package bowtie

import (
	"bytes"
	"fmt"
	"io"
	"iter"
	"os/exec"

	"github.com/fluhus/biostuff/formats/sam"
)

const (
	exe = "bowtie2"
)

// Map runs bowtie on the given fastq file and returns a real-time iterator
// over the resulting SAM lines.
func Map(fq, ref string, threads int, args ...string) iter.Seq2[*sam.SAM, error] {
	return func(yield func(*sam.SAM, error) bool) {
		allArgs := []string{
			"-t", "--no-head", "-p", fmt.Sprint(threads),
			"-x", ref, "-U", fq}
		cmd := exec.Command(exe, append(allArgs, args...)...)
		stderr := bytes.NewBuffer(nil)
		cmd.Stderr = stderr
		r, _ := cmd.StdoutPipe()
		c := make(chan error, 1)
		go func() {
			err := cmd.Run()
			if err != nil {
				err = fmt.Errorf("%w\n%s", err, stderr.Bytes())
			}
			c <- err
		}()
		for sm, err := range sam.Reader(r) {
			if !yield(sm, err) {
				return
			}
		}
		if err := <-c; err != nil {
			yield(nil, err)
		}
	}
}

// MapInt runs bowtie on the given interleaved pairs fastq file
// and returns a real-time iterator over the resulting SAM lines.
func MapInt(fq, ref string, threads int, args ...string) iter.Seq2[*sam.SAM, error] {
	return func(yield func(*sam.SAM, error) bool) {
		allArgs := []string{
			"-t", "--no-head", "-p", fmt.Sprint(threads),
			"-x", ref, "--interleaved", fq}
		cmd := exec.Command(exe, append(allArgs, args...)...)
		stderr := bytes.NewBuffer(nil)
		cmd.Stderr = stderr
		r, _ := cmd.StdoutPipe()
		c := make(chan error, 1)
		go func() {
			err := cmd.Run()
			if err != nil {
				err = fmt.Errorf("%w\n%s", err, stderr.Bytes())
			}
			c <- err
		}()
		for sm, err := range sam.Reader(r) {
			if !yield(sm, err) {
				return
			}
		}
		if err := <-c; err != nil {
			yield(nil, err)
		}
	}
}

// Map2 runs bowtie on the given paired-end fastq files and returns a real-time
// iterator over the resulting SAM lines.
func Map2(fq1, fq2, ref string, threads int, args ...string) iter.Seq2[*sam.SAM, error] {
	return func(yield func(*sam.SAM, error) bool) {
		allArgs := []string{
			"-t", "--no-head", "-p", fmt.Sprint(threads),
			"-x", ref, "-1", fq1, "-2", fq2}
		cmd := exec.Command(exe, append(allArgs, args...)...)
		stderr := bytes.NewBuffer(nil)
		cmd.Stderr = stderr
		r, _ := cmd.StdoutPipe()
		c := make(chan error, 1)
		go func() {
			err := cmd.Run()
			if err != nil {
				err = fmt.Errorf("%w\n%s", err, stderr.Bytes())
			}
			c <- err
		}()
		for sm, err := range sam.Reader(r) {
			if !yield(sm, err) {
				return
			}
		}
		if err := <-c; err != nil {
			yield(nil, err)
		}
	}
}

// MapReader runs bowtie on the given fastq stream
// and returns a real-time iterator over the resulting SAM lines.
func MapReader(fq io.Reader, ref string, threads int, args ...string,
) iter.Seq2[*sam.SAM, error] {
	return func(yield func(*sam.SAM, error) bool) {
		allArgs := []string{
			"-t", "--no-head", "-p", fmt.Sprint(threads),
			"-x", ref, "-U", "-"}
		cmd := exec.Command(exe, append(allArgs, args...)...)
		stderr := bytes.NewBuffer(nil)
		cmd.Stderr = stderr
		cmd.Stdin = fq
		r, _ := cmd.StdoutPipe()
		c := make(chan error, 1)
		go func() {
			err := cmd.Run()
			if err != nil {
				err = fmt.Errorf("%w\n%s", err, stderr.Bytes())
			}
			c <- err
		}()
		for sm, err := range sam.Reader(r) {
			if !yield(sm, err) {
				return
			}
		}
		if err := <-c; err != nil {
			yield(nil, err)
		}
	}
}
