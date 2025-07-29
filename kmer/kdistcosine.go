package kmer

import (
	"errors"
	"math"

	"gonum.org/v1/gonum/mat"
)

type KDistCosine struct {
	Dist float64
}

func NewKDistCosine() *KDistCosine {
	var kde KDistCosine
	kde.Dist = float64(0.0)
	return &kde
}

func (kde *KDistCosine) Compute(a *mat.Dense, b *mat.Dense) error {
	aLen, _ := (*a).Dims()
	bLen, _ := (*b).Dims()
	if aLen != bLen {
		return errors.New("cannot compare vectors of kmer counts with different lengths")
	}

	sub1 := mat.NewDense(1, 1, nil)
	sub1.Reset()
	sub1.MulElem(a, b)
	sumXY := mat.Sum(sub1)

	sub1.Reset()
	sub1.MulElem(a, a)
	sumX := math.Sqrt(mat.Sum(sub1))

	sub1.Reset()
	sub1.MulElem(b, b)
	sumY := math.Sqrt(mat.Sum(sub1))

	kde.Dist = 1 - (sumXY / (sumX * sumY))

	return (nil)
}

func (kde *KDistCosine) GetDistance() float64 {
	return kde.Dist
}

func (kde *KDistCosine) NeedSelfComparison() bool {
	return false
}
