package mybuf

import (
	"io"
	"testing"
)

func TestBuffer(t *testing.T) {
	b := &Buffer{}
	for _, x := range []string{"hi", ",", "hello"} {
		if _, err := b.Write([]byte(x)); err != nil {
			t.Fatalf("Buffer.Write(%q) failed: %v", x, err)
		}
	}
	if err := b.Close(); err != nil {
		t.Fatalf("Buffer.Close() failed: %v", err)
	}
	want := "hi,hello"
	got, err := io.ReadAll(b.Reader())
	if err != nil {
		t.Fatalf("Buffer.Read() failed: %v", err)
	}
	if string(got) != want {
		t.Fatalf("Buffer.Read()=%q, want %q", got, want)
	}
	if err := b.Close(); err != nil {
		t.Fatalf("Buffer.Close() failed: %v", err)
	}
	got, err = io.ReadAll(b.Reader())
	if err != nil {
		t.Fatalf("Buffer.Read() failed: %v", err)
	}
	if string(got) != want {
		t.Fatalf("Buffer.Read()=%q, want %q", got, want)
	}
}
