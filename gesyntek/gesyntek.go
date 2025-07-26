package gesyntek

import (
	"bufio"
	"errors"
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/hdevillers/go-seq/seqio"
)

/*
	GeSynteK class: mainly consists of a collection of loci
*/

// Default values
const (
	WINDOW_LEN int    = 10000
	KMER_LEN   int    = 8
	GFF_TARGET string = "gene"
	GFF_ID     string = "ID"
)

// Structure
type GeSynteK struct {
	WindowLen int
	KmerLen   int
	Loci      []Locus
	SeqIdLoci map[string][]int
	GffTarget string
	GffId     string
}

// Init. GeSynteK object
func NewGeSynteK(w int, k int, t string, i string) *GeSynteK {
	var gsk GeSynteK

	// Initialize arguments
	gsk.WindowLen = w
	gsk.KmerLen = k
	gsk.Loci = make([]Locus, 0)
	gsk.SeqIdLoci = make(map[string][]int)
	gsk.GffTarget = t
	gsk.GffId = i

	return &gsk
}

// Load data from gff files
func (gsk *GeSynteK) LoadGFF(gff string) error {
	// Create the file hendler
	fh, err := os.Open(gff)
	if err != nil {
		return err
	}
	defer fh.Close()

	// Create the buffer
	fb := bufio.NewScanner(fh)

	// Running variable
	iloc := 0

	// Initialize the regex to get the locus name
	re := regexp.MustCompile(gsk.GffId + "=([\\w\\-\\.]+)")

	// Scan each line of the GFF file
	for fb.Scan() {
		err = fb.Err()
		if err != nil {
			return err
		}

		// Read line
		line := fb.Text()
		elem := strings.Split(line, "\t")

		// Check line target
		if elem[2] == gsk.GffTarget {
			// Retrieve the name of the locus
			ln := re.FindStringSubmatch(elem[8])
			if len(ln) == 0 {
				return errors.New("failed to retrieve locus name in the GFF file")
			}

			// Convert location into integer
			start, err := strconv.Atoi(elem[3])
			if err != nil {
				return err
			}
			end, err := strconv.Atoi(elem[4])
			if err != nil {
				return err
			}

			// Create a locus and append
			gsk.Loci = append(gsk.Loci, *NewLocus(elem[0], ln[1], start, end, elem[6]))
			gsk.SeqIdLoci[elem[0]] = append(gsk.SeqIdLoci[elem[0]], iloc)
			iloc++
		} // Else, do nothing
	}

	// Return no error
	return nil
}

// Load up and down stream sequences
func (gsk *GeSynteK) LoadFasta(fasta string) error {
	// Open the fasta file
	seqIn := seqio.NewReader(fasta, "fasta", false)
	seqIn.CheckPanic()
	defer seqIn.Close()

	for seqIn.Next() {
		seqIn.CheckPanic()
		seq := seqIn.Seq()

		if inds, ok := gsk.SeqIdLoci[seq.Id]; ok {
			for i := 0; i < len(inds); i++ {
				gsk.Loci[inds[i]].ExtractUpDownSequence(&seq, gsk.WindowLen)
			}
		}
	}
	return nil
}
