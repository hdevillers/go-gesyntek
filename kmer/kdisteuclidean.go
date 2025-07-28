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

func (kde *KDistEuclidean) Compute(a *[]uint32, b *[]uint32) error {
	aLen := len(*a)
	bLen := len(*b)
	if aLen != bLen {
		return errors.New("cannot compare vectors of kmer counts with different lengths")
	}

	tmp1 := make([]float64, aLen)
	for i := range aLen {
		tmp1[i] = float64((*a)[i])
	}
	va := mat.NewDense(aLen, 1, tmp1)

	tmp2 := make([]float64, bLen)
	for i := range aLen {
		tmp2[i] = float64((*b)[i])
	}
	vb := mat.NewDense(bLen, 1, tmp2)

	sub1 := mat.NewDense(1, 1, nil)
	sub1.Reset()
	sub1.Sub(va, vb)

	sub2 := mat.NewDense(1, 1, nil)
	sub2.Reset()
	sub2.MulElem(sub1, sub1)

	tot1 := mat.Sum(sub2)

	kde.Dist = math.Sqrt(tot1)

	return (nil)
}

func (kde *KDistEuclidean) GetDistance() float64 {
	return kde.Dist
}

func (kde *KDistEuclidean) NeedSelfComparison() bool {
	return false
}
