package main

import (
	"io"
	"iter"

	"github.com/fluhus/biostuff/formats/sam"
	"github.com/fluhus/bundy/mybuf"
	"github.com/fluhus/gostuff/aio"
)

func samWriter() (io.WriteCloser, error) {
	if *diskMode != "" {
		return aio.Create(*diskMode)
	}
	return sambuf, nil
}

func samReader() iter.Seq2[*sam.SAM, error] {
	if *diskMode != "" {
		return sam.IterFile(*diskMode)
	}
	return sam.NewReader(sambuf.Reader()).Iter()
}

var sambuf = &mybuf.Buffer{}
