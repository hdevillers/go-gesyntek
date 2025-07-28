package kmer

type KDist interface {
	Compute(*[]uint32, *[]uint32) error
	GetDistance() float64
	NeedSelfComparison() bool
}
