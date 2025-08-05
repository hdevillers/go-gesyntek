package kmer

import "errors"

type KLabel struct {
	K     int
	Sel   uint64
	Bases []byte
}

func NewKLabel(k int) *KLabel {
	var kl KLabel
	kl.K = k
	kl.Sel = uint64(3)
	kl.Bases = []byte{'A', 'C', 'G', 'T'}
	return &kl
}

func (kl *KLabel) ParseUint64(w uint64, l *[]byte, n int, at int) {
	for i := range n {
		(*l)[at-i] = kl.Bases[int(w&kl.Sel)]
		w = w >> 2
	}
}

func (kl *KLabel) Uint64ToBytes(ws *[][]uint64, out *[][]byte) error {
	nws := len((*ws)[0])
	nout := len(*out)
	if nws != nout {
		return errors.New("wrong number of kmer in type conversion")
	}

	if kl.K <= MaxK64Bits {
		for i := range nws {
			(*out)[i] = make([]byte, kl.K)
			kl.ParseUint64((*ws)[0][i], &(*out)[i], kl.K, kl.K-1)
		}
	} else {
		return errors.New("kmer longer than 32 bases are not supported yet")
	}
	return nil
}

// Merge two list of uint64 => update object a
func (kl *KLabel) MergeUint64(a *[][]uint64, b *[][]uint64) error {
	aLen := len((*a)[0])
	bLen := len((*b)[0])
	tmp := make([]uint64, aLen+bLen)
	ai := 0
	bi := 0
	n := 0
	for ai < aLen && bi < bLen {
		if (*a)[0][ai] == (*b)[0][bi] {
			tmp[n] = (*a)[0][ai]
			n++
			ai++
			bi++
		} else if (*a)[0][ai] < (*b)[0][bi] {
			tmp[n] = (*a)[0][ai]
			n++
			ai++
		} else {
			tmp[n] = (*b)[0][bi]
			n++
			bi++
		}
	}
	for ai < aLen {
		tmp[n] = (*a)[0][ai]
		n++
		ai++
	}
	for bi < bLen {
		tmp[n] = (*b)[0][bi]
		n++
		bi++
	}
	(*a)[0] = make([]uint64, n)
	for i := range n {
		(*a)[0][i] = tmp[i]
	}

	return nil
}
