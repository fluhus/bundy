// An abstraction over SAM read/write (memory vs disk).

package main

import (
	"io"
	"iter"

	"github.com/fluhus/biostuff/formats/sam"
	"github.com/fluhus/bundy/mybuf"
	"github.com/fluhus/gostuff/aio"
)

// Returns a writer based on the diskmode argument.
func samWriter() (io.WriteCloser, error) {
	if *diskMode != "" {
		return aio.Create(*diskMode)
	}
	return sambuf, nil
}

// Returns a reader based on the diskmode argument.
func samReader() iter.Seq2[*sam.SAM, error] {
	if *diskMode != "" {
		return sam.File(*diskMode)
	}
	return sam.Reader(sambuf.Reader())
}

// A buffer for memory mode.
var sambuf = &mybuf.Buffer{}
