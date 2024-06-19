package common

import (
	"fmt"
	"os"
)

// Die prints the error and exits if the error is non-nil.
func Die(err error) {
	if err != nil {
		fmt.Fprintln(os.Stderr, "ERROR:", err)
		os.Exit(2)
	}
}

// Perc returns a/b in %.
func Perc(a, b int) float64 {
	return 100 * float64(a) / float64(b)
}

// Percf returns a/b in the format "x%" with the given precision.
func Percf(a, b, p int) string {
	return fmt.Sprintf(fmt.Sprintf("%%.%df%%%%", p), Perc(a, b))
}
