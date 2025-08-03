package kmer

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

type KCountSmall struct {
	K         int
	Canonical bool
	Convert   []uint32
	RcConvert []uint32
	Fwd       uint32
	Bwd       uint32
	Kmers     [][]uint64
	Counts    mat.Dense
	SkipDeg   int
	SkipShort int
}

func NewKCountSmall(K int, c bool) *KCountSmall {
	var kcs KCountSmall

	nKmers := int(math.Pow(4.0, float64(K)))

	kcs.K = K
	kcs.Canonical = c
	kcs.Fwd = uint32((16 - K + 1) * 2)
	kcs.Bwd = uint32((16 - K) * 2)
	kcs.SkipDeg = 0
	kcs.SkipShort = 0
	kcs.Counts = *mat.NewDense(nKmers, 1, nil)
	kcs.Kmers = make([][]uint64, 1)
	kcs.Kmers[0] = make([]uint64, nKmers)
	kcs.Convert = make([]uint32, 256)

	// Set up conversion
	kcs.Convert['C'] = uint32(1)
	kcs.Convert['c'] = uint32(1)
	kcs.Convert['G'] = uint32(2)
	kcs.Convert['g'] = uint32(2)
	kcs.Convert['T'] = uint32(3)
	kcs.Convert['t'] = uint32(3)

	// Set Kmers
	for i := range nKmers {
		kcs.Kmers[0][i] = uint64(i)
	}

	// Set reverse complement conversion if necessary
	if c {
		kcs.RcConvert = make([]uint32, 256)
		n := 2*K - 2
		kcs.RcConvert['C'] = uint32(2) << n
		kcs.RcConvert['c'] = uint32(2) << n
		kcs.RcConvert['G'] = uint32(1) << n
		kcs.RcConvert['g'] = uint32(1) << n
		kcs.RcConvert['A'] = uint32(3) << n
		kcs.RcConvert['a'] = uint32(3) << n
	}

	return &kcs
}

func (kcs *KCountSmall) Count(seq *[]byte) error {
	// Split the input sequence to keep only countable words
	seqSpl := NewKSplit(kcs.K)
	err := seqSpl.SplitSeq(seq)
	if err != nil {
		return err
	}

	// Retrieve the number of skipped bases
	kcs.SkipDeg = seqSpl.NSkipped
	kcs.SkipShort = seqSpl.NTooShort

	// Init. count variable
	cnt := make([]float64, len(kcs.Kmers[0]))

	// Count words
	if kcs.Canonical {
		for iSeq := range len(seqSpl.SeqSplit) {
			// Initialize the first word
			w := uint32(0)
			wrc := uint32(0)
			for i := 0; i < kcs.K; i++ {
				w = (w << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
				wrc = (wrc >> 2) | kcs.RcConvert[seqSpl.SeqSplit[iSeq][i]]
			}

			// Count the first word
			if w < wrc {
				cnt[w]++
			} else {
				cnt[wrc]++
			}

			// Continue with the following word
			for i := kcs.K; i < len(seqSpl.SeqSplit[iSeq]); i++ {
				w = (w<<kcs.Fwd)>>kcs.Bwd | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
				wrc = (wrc >> 2) | kcs.RcConvert[seqSpl.SeqSplit[iSeq][i]]
				if w < wrc {
					cnt[w]++
				} else {
					cnt[wrc]++
				}
			}
		}
	} else {
		for iSeq := range len(seqSpl.SeqSplit) {
			// Initialize the first word
			w := uint32(0)
			for i := 0; i < kcs.K; i++ {
				w = (w << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
			}

			// Count the first word
			cnt[w]++

			// Continue with the following word
			for i := kcs.K; i < len(seqSpl.SeqSplit[iSeq]); i++ {
				w = (w<<kcs.Fwd)>>kcs.Bwd | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
				cnt[w]++
			}
		}
	}

	// Set the count variable
	kcs.Counts.SetCol(0, cnt)

	return nil
}

func (kcs *KCountSmall) MergeKmers(kl *[][]uint64) error {
	// Useless function for small Kmer
	return nil
}

func (kcs *KCountSmall) GetSkippedDegeneratedBases() int {
	return kcs.SkipDeg
}

func (kcs *KCountSmall) GetSkippedTooShortBases() int {
	return kcs.SkipShort
}

func (kcs *KCountSmall) GetSkippedBases() int {
	return kcs.SkipDeg + kcs.SkipShort
}

func (kcs *KCountSmall) GetKmers() *[][]uint64 {
	return &kcs.Kmers
}

func (kcs *KCountSmall) GetCounts() *mat.Dense {
	return &kcs.Counts
}

func (kcs *KCountSmall) GetNKmers() int {
	return len(kcs.Kmers[0])
}

func (kcs *KCountSmall) IsCanonical() bool {
	return kcs.Canonical
}
