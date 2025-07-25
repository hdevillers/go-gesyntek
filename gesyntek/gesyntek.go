package gesyntek

// DEFINING STRUCTURES

// GeSynteK unitary object (one locus data)
type GeSynteK struct {
	SeqId     string
	SeqLabel  string
	SeqStart  string
	SeqEnd    string
	SeqStrand string
}

// Init. GeSynteK object
func NewGeSynteK(i, l, s, e, p string) *GeSynteK {
	var gsk GeSynteK

	// Init. variable
	gsk.SeqId = i
	gsk.SeqLabel = l
	gsk.SeqStart = s
	gsk.SeqEnd = e
	gsk.SeqStrand = p

	return &gsk
}
