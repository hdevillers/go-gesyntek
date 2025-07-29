package kmer

import "gonum.org/v1/gonum/mat"

type KDist interface {
	Compute(*mat.Dense, *mat.Dense) error
	GetDistance() float64
	NeedSelfComparison() bool
}
