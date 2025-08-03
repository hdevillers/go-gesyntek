package kmer

import "gonum.org/v1/gonum/mat"

const (
	MaxKSmall    int = 8
	MaxK32Bits   int = 15
	MaxK64Bits   int = 31
	MaxKPrintAll int = 10
	MaxKAbsolute int = 31
)

type KCount interface {
	Count(*[]byte) error
	MergeKmers(*[][]uint64) error
	GetSkippedBases() int
	GetSkippedDegeneratedBases() int
	GetSkippedTooShortBases() int
	GetCounts() *mat.Dense
	GetKmers() *[][]uint64
	GetNKmers() int
	IsCanonical() bool
}
