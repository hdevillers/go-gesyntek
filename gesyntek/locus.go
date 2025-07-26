package gesyntek

/*
	Locus class: contain data for a single locus
*/

// Structure
type Locus struct {
	SeqId     string
	SeqLabel  string
	SeqStart  string
	SeqEnd    string
	SeqStrand string
}

// Constructor
func NewLocus(i, l, s, e, p string) *Locus {
	var locus Locus

	// Init. variable
	locus.SeqId = i
	locus.SeqLabel = l
	locus.SeqStart = s
	locus.SeqEnd = e
	locus.SeqStrand = p

	return &locus
}
