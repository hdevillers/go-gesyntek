package kmer

import (
	"errors"
	"math"

	"gonum.org/v1/gonum/mat"
)

type KDistEuclidean struct {
	Dist float64
}

func NewKDistEuclidean() *KDistEuclidean {
	var kde KDistEuclidean
	kde.Dist = float64(0.0)
	return &kde
}

func (kde *KDistEuclidean) Compute(a *mat.Dense, b *mat.Dense) error {
	aLen, _ := (*a).Dims()
	bLen, _ := (*b).Dims()
	if aLen != bLen {
		return errors.New("cannot compare vectors of kmer counts with different lengths")
	}

	// Compute the differences
	sub1 := mat.NewDense(1, 1, nil)
	sub1.Reset()
	sub1.Sub(a, b)

	// Compute the square
	sub2 := mat.NewDense(1, 1, nil)
	sub2.Reset()
	sub2.MulElem(sub1, sub1)

	// Compute the sum
	tot1 := mat.Sum(sub2)

	// Compute the square root
	kde.Dist = math.Sqrt(tot1)

	return (nil)
}

func (kde *KDistEuclidean) GetDistance() float64 {
	return kde.Dist
}

func (kde *KDistEuclidean) NeedSelfComparison() bool {
	return false
}
