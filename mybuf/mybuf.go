package mybuf

import (
	"bytes"
	"fmt"
	"io"
	"slices"

	"github.com/klauspost/compress/zstd"
)

type Buffer struct {
	data []byte
	buf  *bytes.Buffer
	zw   io.WriteCloser
}

func (b *Buffer) Write(p []byte) (n int, err error) {
	if b.zw == nil {
		b.buf = bytes.NewBuffer(nil)
		b.zw, _ = zstd.NewWriter(b.buf, zstd.WithEncoderConcurrency(1),
			zstd.WithEncoderLevel(1))
	}
	return b.zw.Write(p)
}

func (b *Buffer) Close() error {
	if b.zw != nil {
		if err := b.zw.Close(); err != nil {
			return err
		}
		b.data = slices.Clip(b.buf.Bytes())
		b.zw = nil
		b.buf = nil
	}
	return nil
}

func (b *Buffer) Reader() io.Reader {
	if b.zw != nil {
		return &errReader{fmt.Errorf("called Read without closing the writer first")}
	}
	r, err := zstd.NewReader(bytes.NewBuffer(b.data), zstd.WithDecoderConcurrency(1))
	if err != nil {
		return &errReader{err}
	}
	return r
}

type errReader struct {
	err error
}

func (r *errReader) Read([]byte) (int, error) {
	return 0, r.err
}
