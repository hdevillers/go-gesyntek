package kmer

import "gonum.org/v1/gonum/mat"

const (
	MaxKSmall    int = 10
	MaxK32Bits   int = 15
	MaxKPrintAll int = 10
	MaxKAbsolute int = 15
)

type KCount interface {
	Count(*[]byte) error
	GetSkippedBases() int
	GetSkippedDegeneratedBases() int
	GetSkippedTooShortBases() int
	GetCounts() *mat.Dense
	GetKmers() *[][]uint64
	GetNKmers() int
	NeedToMerge() bool
}
