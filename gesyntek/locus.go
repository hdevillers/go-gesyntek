package gesyntek

import (
	"fmt"

	"github.com/hdevillers/go-seq/seq"
)

/*
	Locus class: contain data for a single locus
*/

// Structure
type Locus struct {
	SeqId      string
	SeqLabel   string
	SeqStart   int
	SeqEnd     int
	SeqStrand  string
	SeqUpStr   seq.Seq
	SeqDownStr seq.Seq
	HasUpStr   bool
	HasDownStr bool
}

// Constructor
func NewLocus(i string, l string, s int, e int, p string) *Locus {
	var locus Locus

	// Init. variable
	locus.SeqId = i
	locus.SeqLabel = l
	locus.SeqStart = s
	locus.SeqEnd = e
	locus.SeqStrand = p
	locus.HasUpStr = false
	locus.HasDownStr = false

	return &locus
}

// Extract up- and down-stream sequences
func (locus *Locus) ExtractUpDownSequence(s *seq.Seq, w int) error {
	// Extract sub-sequence on the 'left'
	if locus.SeqStart > w {
		leftDNA := make([]byte, w)
		from := locus.SeqStart - w - 1
		to := locus.SeqStart - 1
		copy(leftDNA, s.Sequence[from:to])

		if locus.SeqStrand == "+" {
			// This is up-stream sequence
			locus.HasUpStr = true
			id := fmt.Sprintf("%s_upstream_w%d", locus.SeqId, w)
			locus.SeqUpStr = seq.Seq{Id: id, Sequence: leftDNA}
		} else {
			// This is down-stream and sequence has to be rev-comp
			leftDNArc := make([]byte, w) // To be included into go-seq module!
			rc := make([]byte, 85)
			rc[65] = 'T'
			rc[67] = 'G'
			rc[71] = 'C'
			rc[84] = 'A'
			for i := range w {
				leftDNArc[w-i-1] = rc[leftDNA[i]]
			}
			locus.HasDownStr = true
			id := fmt.Sprintf("%s_downstream_w%d", locus.SeqId, w)
			locus.SeqDownStr = seq.Seq{Id: id, Sequence: leftDNArc}
		}
	} // Else nothing to do (or create an empty sequence with w*N?)

	// Extract sub_sequence on the 'right'
	if locus.SeqEnd+w < s.Length() {
		rightDNA := make([]byte, w)
		from := locus.SeqEnd
		to := from + w - 1
		copy(rightDNA, s.Sequence[from:to])

		if locus.SeqStrand == "+" {
			// This is down-stream sequence
			locus.HasDownStr = true
			id := fmt.Sprintf("%s_downstream_w%d", locus.SeqId, w)
			locus.SeqDownStr = seq.Seq{Id: id, Sequence: rightDNA}
		} else {
			// This is up-stream and sequence has to be rev-comp
			rightDNArc := make([]byte, w)
			rc := make([]byte, 85)
			rc[65] = 'T'
			rc[67] = 'G'
			rc[71] = 'C'
			rc[84] = 'A'
			for i := range w {
				rightDNArc[w-i-1] = rc[rightDNA[i]]
			}
			locus.HasUpStr = true
			id := fmt.Sprintf("%s_upstream_w%d", locus.SeqId, w)
			locus.SeqUpStr = seq.Seq{Id: id, Sequence: rightDNArc}
		}
	} // Else nothing...
	return nil
}
