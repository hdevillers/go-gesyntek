package kmer

import (
	"slices"

	"gonum.org/v1/gonum/mat"
)

type KCount31 struct {
	K         int
	Convert   []uint64
	Fwd       uint64
	Bwd       uint64
	Kmers     [][]uint64
	Counts    mat.Dense
	SkipDeg   int
	SkipShort int
}

func NewKCount31(K int) *KCount31 {
	var kcs KCount31

	kcs.K = K
	kcs.Fwd = uint64((32 - K + 1) * 2)
	kcs.Bwd = uint64((32 - K) * 2)
	kcs.SkipDeg = 0
	kcs.SkipShort = 0
	kcs.Counts = *mat.NewDense(1, 1, nil)
	kcs.Kmers = make([][]uint64, 1)
	kcs.Convert = make([]uint64, 256)

	// Set up conversion
	kcs.Convert['C'] = uint64(1)
	kcs.Convert['c'] = uint64(1)
	kcs.Convert['G'] = uint64(2)
	kcs.Convert['g'] = uint64(2)
	kcs.Convert['T'] = uint64(3)
	kcs.Convert['t'] = uint64(3)

	return &kcs
}

func (kcs *KCount31) Count(seq *[]byte) error {
	// Split the input sequence to keep only countable words
	seqSpl := NewKSplit(kcs.K)
	err := seqSpl.SplitSeq(seq)
	if err != nil {
		return err
	}

	// Retrieve the number of skipped bases
	kcs.SkipDeg = seqSpl.NSkipped
	kcs.SkipShort = seqSpl.NTooShort

	// Count the number of kmers from seqSpl
	nKmers := 0
	for i := range len(seqSpl.SeqSplit) {
		nKmers += len(seqSpl.SeqSplit[i]) - kcs.K + 1
	}

	// Enumerate kmers
	lab := make([]uint64, nKmers)
	iLab := 0
	for iSeq := range len(seqSpl.SeqSplit) {
		// Initialize the first word
		w := uint64(0)
		for i := 0; i < kcs.K; i++ {
			w = (w << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
		}

		// Store the first kmer
		lab[iLab] = w
		iLab++

		// Continue with the following kmer(s)
		for i := kcs.K; i < len(seqSpl.SeqSplit[iSeq]); i++ {
			w = (w<<kcs.Fwd)>>kcs.Bwd | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
			lab[iLab] = w
			iLab++
		}
	}

	// Sort kmers
	slices.Sort(lab)

	// Count kmers
	tmpLab := make([]uint64, nKmers)
	tmpCnt := make([]float64, nKmers)
	iKmer := 0
	for i := 0; i < nKmers; {
		tmpLab[iKmer] = lab[i]
		cnt := 1
		for j := i + 1; j < nKmers && lab[i] == lab[j]; {
			cnt++
			j++
		}
		tmpCnt[iKmer] = float64(cnt)
		iKmer++
		i += cnt
	}

	// Copy labels and counts
	kcs.Kmers[0] = make([]uint64, iKmer)
	copy(kcs.Kmers[0], tmpLab[0:iKmer])
	kcs.Counts = *mat.NewDense(iKmer, 1, nil)
	kcs.Counts.SetCol(0, tmpCnt[0:iKmer])

	return nil
}

func (kcs *KCount31) MergeKmers(kl *[][]uint64) error {
	nIn := len(kcs.Kmers[0])
	nOut := len((*kl)[0])
	tmpLab := make([]uint64, nIn+nOut)
	tmpCnt := make([]float64, nIn+nOut)
	tmpI := 0
	i := 0
	j := 0
	for i < nIn && j < nOut {
		if kcs.Kmers[0][i] == (*kl)[0][j] {
			tmpLab[tmpI] = kcs.Kmers[0][i]
			tmpCnt[tmpI] = kcs.Counts.At(i, 0)
			tmpI++
			i++
			j++
		} else if kcs.Kmers[0][i] < (*kl)[0][j] {
			tmpLab[tmpI] = kcs.Kmers[0][i]
			tmpCnt[tmpI] = kcs.Counts.At(i, 0)
			tmpI++
			i++
		} else {
			tmpLab[tmpI] = (*kl)[0][j]
			tmpCnt[tmpI] = 0.0
			tmpI++
			j++
		}
	}
	for i < nIn {
		tmpLab[tmpI] = kcs.Kmers[0][i]
		tmpCnt[tmpI] = kcs.Counts.At(i, 0)
		tmpI++
		i++
	}
	for j < nOut {
		tmpLab[tmpI] = (*kl)[0][j]
		tmpCnt[tmpI] = 0.0
		tmpI++
		j++
	}

	// Copy new data
	kcs.Kmers[0] = make([]uint64, tmpI)
	copy(kcs.Kmers[0], tmpLab[0:tmpI])
	kcs.Counts = *mat.NewDense(tmpI, 1, nil)
	kcs.Counts.SetCol(0, tmpCnt[0:tmpI])

	return nil
}

func (kcs *KCount31) GetSkippedDegeneratedBases() int {
	return kcs.SkipDeg
}

func (kcs *KCount31) GetSkippedTooShortBases() int {
	return kcs.SkipShort
}

func (kcs *KCount31) GetSkippedBases() int {
	return kcs.SkipDeg + kcs.SkipShort
}

func (kcs *KCount31) GetKmers() *[][]uint64 {
	return &kcs.Kmers
}

func (kcs *KCount31) GetCounts() *mat.Dense {
	return &kcs.Counts
}

func (ksc *KCount31) GetNKmers() int {
	return len(ksc.Kmers[0])
}
