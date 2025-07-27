package kmer

import "errors"

type KSplit struct {
	K         int
	KeptBases []int
	SeqSplit  [][]byte
	NSkipped  int
	NTooShort int
}

func NewKSplit(K int) *KSplit {
	var ks KSplit

	ks.K = K
	ks.NSkipped = 0
	ks.NTooShort = 0
	ks.KeptBases = make([]int, 256)
	ks.KeptBases['A'] = 1
	ks.KeptBases['a'] = 1
	ks.KeptBases['C'] = 1
	ks.KeptBases['c'] = 1
	ks.KeptBases['G'] = 1
	ks.KeptBases['g'] = 1
	ks.KeptBases['T'] = 1
	ks.KeptBases['t'] = 1
	ks.SeqSplit = [][]byte{}

	return &ks
}

func (ks *KSplit) SplitSeq(seq *[]byte) error {
	iLoc := 0
	sLen := len(*seq)

	// Continue
	for iLoc < sLen {
		for (iLoc < sLen) && (ks.KeptBases[(*seq)[iLoc]] == 0) {
			iLoc++
			ks.NSkipped++
		}
		from := iLoc
		for (iLoc < sLen) && (ks.KeptBases[(*seq)[iLoc]] == 1) {
			iLoc++
		}
		iLen := iLoc - from
		if iLen >= ks.K {
			// Keep this fragment
			ks.SeqSplit = append(ks.SeqSplit, (*seq)[from:iLoc])
		} else {
			ks.NTooShort += iLen
		}
	}

	if len(ks.SeqSplit) == 0 {
		return errors.New("no sequence kept (Kmer SplitSeq)")
	}
	return nil
}
