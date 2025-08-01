package kmer

import (
	"errors"
	"math"

	"gonum.org/v1/gonum/mat"
)

type KDistMash struct {
	Dist    float64
	K       int
	Epsilon float64
}

func NewKDistMash(k int) *KDistMash {
	var kdm KDistMash
	kdm.Dist = float64(0.0)
	kdm.K = k
	kdm.Epsilon = 1e-7
	return &kdm
}

func isNotNull(i, j int, v float64) float64 {
	if v != 0 {
		return 1.0
	}
	return 0.0
}

func (kdm *KDistMash) Compute(a *mat.Dense, b *mat.Dense) error {
	aLen, _ := (*a).Dims()
	bLen, _ := (*b).Dims()
	if aLen != bLen {
		return errors.New("cannot compare vectors of kmer counts with different lengths")
	}

	// Compute Jaccard index
	// First identify non-nul counts
	aNotNul := mat.NewDense(aLen, 1, nil)
	aNotNul.Apply(isNotNull, a)
	bNotNul := mat.NewDense(bLen, 1, nil)
	bNotNul.Apply(isNotNull, b)

	// Intersection
	abInter1 := mat.NewDense(aLen, 1, nil)
	abInter1.MulElem(aNotNul, bNotNul)
	abInter2 := mat.Sum(abInter1)

	// Union
	abUnion1 := mat.NewDense(aLen, 1, nil)
	abUnion1.Add(aNotNul, bNotNul)
	abUnion2 := mat.NewDense(aLen, 1, nil)
	abUnion2.Apply(isNotNull, abUnion1)
	abUnion3 := mat.Sum(abUnion2)

	// The J index
	J := abInter2 / abUnion3

	// Mash distance
	tmp := (2.0*J)/(1.0+J) + kdm.Epsilon
	kdm.Dist = -1.0 / float64(kdm.K) * math.Log(tmp)

	return nil
}

func (kdm *KDistMash) GetDistance() float64 {
	return kdm.Dist
}

func (kdm *KDistMash) NeedSelfComparison() bool {
	return false
}
