package kmer

import (
	"cmp"
	"slices"

	"gonum.org/v1/gonum/mat"
)

type KLab64 struct {
	w1 uint64
	w2 uint64
}

func MinKlab64(a, b KLab64) KLab64 {
	if a.w1 < b.w1 {
		return a
	} else if a.w1 == b.w1 {
		if a.w2 < b.w2 {
			return a
		} else {
			return b
		}
	}
	return b
}

type KCount64 struct {
	K         int
	SubK      int
	Canonical bool
	Convert   []uint64
	RcConvert [][]uint64
	Fwd       uint64
	Bwd       uint64
	Kmers     [][]uint64
	Counts    mat.Dense
	SkipDeg   int
	SkipShort int
}

func NewKCount64(K int, c bool) *KCount64 {
	var kcs KCount64

	kcs.K = K
	kcs.SubK = K - 32
	kcs.Canonical = c
	kcs.Fwd = uint64((64 - K + 1) * 2)
	kcs.Bwd = uint64((64 - K) * 2)
	kcs.SkipDeg = 0
	kcs.SkipShort = 0
	kcs.Counts = *mat.NewDense(1, 1, nil)
	kcs.Kmers = make([][]uint64, 2)
	kcs.Convert = make([]uint64, 256)

	// Set up conversion
	kcs.Convert['C'] = uint64(1)
	kcs.Convert['c'] = uint64(1)
	kcs.Convert['G'] = uint64(2)
	kcs.Convert['g'] = uint64(2)
	kcs.Convert['T'] = uint64(3)
	kcs.Convert['t'] = uint64(3)

	// Set reverse complement conversion if necessary
	if c {
		kcs.RcConvert = make([][]uint64, 3)
		kcs.RcConvert[0] = make([]uint64, 256)
		kcs.RcConvert[1] = make([]uint64, 256)
		kcs.RcConvert[2] = make([]uint64, 256)
		n := 2*(K-32) - 2
		kcs.RcConvert[0]['C'] = uint64(2)
		kcs.RcConvert[0]['c'] = uint64(2)
		kcs.RcConvert[0]['G'] = uint64(1)
		kcs.RcConvert[0]['g'] = uint64(1)
		kcs.RcConvert[0]['A'] = uint64(3)
		kcs.RcConvert[0]['a'] = uint64(3)
		kcs.RcConvert[1]['C'] = uint64(2) << n
		kcs.RcConvert[1]['c'] = uint64(2) << n
		kcs.RcConvert[1]['G'] = uint64(1) << n
		kcs.RcConvert[1]['g'] = uint64(1) << n
		kcs.RcConvert[1]['A'] = uint64(3) << n
		kcs.RcConvert[1]['a'] = uint64(3) << n
		kcs.RcConvert[2]['C'] = uint64(2) << 62
		kcs.RcConvert[2]['c'] = uint64(2) << 62
		kcs.RcConvert[2]['G'] = uint64(1) << 62
		kcs.RcConvert[2]['g'] = uint64(1) << 62
		kcs.RcConvert[2]['A'] = uint64(3) << 62
		kcs.RcConvert[2]['a'] = uint64(3) << 62
	}

	return &kcs
}

func (kcs *KCount64) Count(seq *[]byte) error {
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
	lab := make([]KLab64, nKmers)
	iLab := 0
	if kcs.Canonical {
		for iSeq := range len(seqSpl.SeqSplit) {
			// Initialize the first word
			w := KLab64{0.0, 0.0}
			wrc := KLab64{0.0, 0.0}
			// First sub word
			for i := 0; i < kcs.SubK; i++ {
				w.w1 = (w.w1 << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
				wrc.w1 = (wrc.w1 << 2) | kcs.RcConvert[0][seqSpl.SeqSplit[iSeq][kcs.K-i-1]]
			}
			// Second word
			for i := kcs.SubK; i < kcs.K; i++ {
				w.w2 = (w.w2 << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
				wrc.w2 = (wrc.w2 << 2) | kcs.RcConvert[0][seqSpl.SeqSplit[iSeq][kcs.K-i-1]]
			}

			// Store the first kmer
			lab[iLab] = MinKlab64(w, wrc)
			iLab++

			// Continue with the following kmer(s)
			for i := kcs.K; i < len(seqSpl.SeqSplit[iSeq]); i++ {
				w.w1 = (w.w1<<kcs.Fwd)>>kcs.Bwd | kcs.Convert[seqSpl.SeqSplit[iSeq][i-32]]
				w.w2 = (w.w2 << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
				wrc.w1 = (wrc.w1 >> 2) | kcs.RcConvert[1][seqSpl.SeqSplit[iSeq][i]]
				wrc.w2 = (wrc.w2 >> 2) | kcs.RcConvert[2][seqSpl.SeqSplit[iSeq][i-kcs.SubK]]
				lab[iLab] = MinKlab64(w, wrc)
				iLab++
			}
		}
	} else {
		for iSeq := range len(seqSpl.SeqSplit) {
			// Initialize the first word
			w := KLab64{0.0, 0.0}
			// First sub word
			for i := 0; i < kcs.SubK; i++ {
				w.w1 = (w.w1 << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
			}
			// Second word
			for i := kcs.SubK; i < kcs.K; i++ {
				w.w2 = (w.w2 << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
			}

			// Store the first kmer
			lab[iLab] = w
			iLab++

			// Continue with the following kmer(s)
			for i := kcs.K; i < len(seqSpl.SeqSplit[iSeq]); i++ {
				w.w1 = (w.w1<<kcs.Fwd)>>kcs.Bwd | kcs.Convert[seqSpl.SeqSplit[iSeq][i-32]]
				w.w2 = (w.w2 << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
				lab[iLab] = w
				iLab++
			}
		}
	}

	// Sort Kmers
	slices.SortFunc(lab, func(a, b KLab64) int {
		return cmp.Or(cmp.Compare(a.w1, b.w1), cmp.Compare(a.w2, b.w2))
	})

	// Count kmers
	tmpLab := make([][]uint64, 2)
	tmpLab[0] = make([]uint64, nKmers)
	tmpLab[1] = make([]uint64, nKmers)
	tmpCnt := make([]float64, nKmers)
	iKmer := 0
	for i := 0; i < nKmers; {
		tmpLab[0][iKmer] = lab[i].w1
		tmpLab[1][iKmer] = lab[i].w2
		cnt := 1
		for j := i + 1; j < nKmers && lab[i].w1 == lab[j].w1 && lab[i].w2 == lab[j].w2; {
			cnt++
			j++
		}
		tmpCnt[iKmer] = float64(cnt)
		iKmer++
		i += cnt
	}

	// Copy labels and counts
	kcs.Kmers[0] = make([]uint64, iKmer)
	kcs.Kmers[1] = make([]uint64, iKmer)
	copy(kcs.Kmers[0], tmpLab[0][0:iKmer])
	copy(kcs.Kmers[1], tmpLab[1][0:iKmer])
	kcs.Counts = *mat.NewDense(iKmer, 1, nil)
	kcs.Counts.SetCol(0, tmpCnt[0:iKmer])

	return nil
}

func (kcs *KCount64) MergeKmers(kl *[][]uint64) error {
	nIn := len(kcs.Kmers[0])
	nOut := len((*kl)[0])
	tmpLab := make([][]uint64, 2)
	tmpLab[0] = make([]uint64, nIn+nOut)
	tmpLab[1] = make([]uint64, nIn+nOut)
	tmpCnt := make([]float64, nIn+nOut)
	tmpI := 0
	i := 0
	j := 0
	for i < nIn && j < nOut {
		if kcs.Kmers[0][i] == (*kl)[0][j] {
			if kcs.Kmers[1][i] == (*kl)[1][j] {
				tmpLab[0][tmpI] = kcs.Kmers[0][i]
				tmpLab[1][tmpI] = kcs.Kmers[1][i]
				tmpCnt[tmpI] = kcs.Counts.At(i, 0)
				tmpI++
				i++
				j++
			} else if kcs.Kmers[1][i] < (*kl)[1][j] {
				tmpLab[0][tmpI] = kcs.Kmers[0][i]
				tmpLab[1][tmpI] = kcs.Kmers[1][i]
				tmpCnt[tmpI] = kcs.Counts.At(i, 0)
				tmpI++
				i++
			} else {
				tmpLab[0][tmpI] = (*kl)[0][j]
				tmpLab[1][tmpI] = (*kl)[1][j]
				tmpCnt[tmpI] = 0.0
				tmpI++
				j++
			}
		} else if kcs.Kmers[0][i] < (*kl)[0][j] {
			tmpLab[0][tmpI] = kcs.Kmers[0][i]
			tmpLab[1][tmpI] = kcs.Kmers[1][i]
			tmpCnt[tmpI] = kcs.Counts.At(i, 0)
			tmpI++
			i++
		} else {
			tmpLab[0][tmpI] = (*kl)[0][j]
			tmpLab[1][tmpI] = (*kl)[1][j]
			tmpCnt[tmpI] = 0.0
			tmpI++
			j++
		}
	}
	for i < nIn {
		tmpLab[0][tmpI] = kcs.Kmers[0][i]
		tmpLab[1][tmpI] = kcs.Kmers[1][i]
		tmpCnt[tmpI] = kcs.Counts.At(i, 0)
		tmpI++
		i++
	}
	for j < nOut {
		tmpLab[0][tmpI] = (*kl)[0][j]
		tmpLab[1][tmpI] = (*kl)[1][j]
		tmpCnt[tmpI] = 0.0
		tmpI++
		j++
	}

	// Copy new data
	kcs.Kmers[0] = make([]uint64, tmpI)
	kcs.Kmers[1] = make([]uint64, tmpI)
	copy(kcs.Kmers[0], tmpLab[0][0:tmpI])
	copy(kcs.Kmers[1], tmpLab[1][0:tmpI])
	kcs.Counts = *mat.NewDense(tmpI, 1, nil)
	kcs.Counts.SetCol(0, tmpCnt[0:tmpI])

	return nil
}

func (kcs *KCount64) GetSkippedDegeneratedBases() int {
	return kcs.SkipDeg
}

func (kcs *KCount64) GetSkippedTooShortBases() int {
	return kcs.SkipShort
}

func (kcs *KCount64) GetSkippedBases() int {
	return kcs.SkipDeg + kcs.SkipShort
}

func (kcs *KCount64) GetKmers() *[][]uint64 {
	return &kcs.Kmers
}

func (kcs *KCount64) GetCounts() *mat.Dense {
	return &kcs.Counts
}

func (kcs *KCount64) GetNKmers() int {
	return len(kcs.Kmers[0])
}

func (kcs *KCount64) IsCanonical() bool {
	return kcs.Canonical
}

func (kcs *KCount64) GetKmersToSkip() *[]uint8 {
	out := make([]uint8, kcs.GetNKmers())
	return &out
}
