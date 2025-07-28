package kmer

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
	GetCounts() *[]uint32
	GetKmers() *[]uint32
	GetNKmers() int
	NeedToMerge() bool
}
