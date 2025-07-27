package kmer

import (
	"math"
)

type KCountSmall struct {
	K         int
	Convert   []uint32
	Fwd       uint32
	Bwd       uint32
	Counts    []uint32
	SkipDeg   int
	SkipShort int
}

func NewKCountSmall(K int) *KCountSmall {
	var kcs KCountSmall

	kcs.K = K
	kcs.Fwd = uint32((16 - K + 1) * 2)
	kcs.Bwd = uint32((16 - K) * 2)
	kcs.SkipDeg = 0
	kcs.SkipShort = 0
	kcs.Counts = make([]uint32, int(math.Pow(4.0, float64(K))))
	kcs.Convert = make([]uint32, 256)

	// Set up conversion
	kcs.Convert['C'] = uint32(1)
	kcs.Convert['c'] = uint32(1)
	kcs.Convert['G'] = uint32(2)
	kcs.Convert['g'] = uint32(2)
	kcs.Convert['T'] = uint32(3)
	kcs.Convert['t'] = uint32(3)

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

	// Count words
	for iSeq := range len(seqSpl.SeqSplit) {
		// Initialize the first word
		w := uint32(0)
		for i := 0; i < kcs.K; i++ {
			w = (w << 2) | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
		}

		// Count the first word
		kcs.Counts[w]++

		// Continue with the following word
		for i := kcs.K; i < len(seqSpl.SeqSplit[iSeq]); i++ {
			w = (w<<kcs.Fwd)>>kcs.Bwd | kcs.Convert[seqSpl.SeqSplit[iSeq][i]]
			kcs.Counts[w]++
		}
	}

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
