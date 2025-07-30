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

func (kl *KLabel) Parse64Bits(w uint64, l *[]byte, n int, at int) {
	for i := range n {
		(*l)[at+i] = kl.Bases[int(w&kl.Sel)]
		w = w >> 2
	}
}

func (kl *KLabel) Uint64ToBytes(ws *[][]uint64, out *[][]byte) error {
	nws := len(*ws)
	nout := len(*out)
	if nws != nout {
		return errors.New("wrong number of kmer in type conversion")
	}

	if kl.K < 32 {
		for i := range nws {
			(*out)[i] = make([]byte, kl.K)
			kl.Parse64Bits((*ws)[0][i], &(*out)[i], kl.K, 0)
		}
	} else {
		return errors.New("kmer longer than 31 bases are not supported yet")
	}
	return nil
}
