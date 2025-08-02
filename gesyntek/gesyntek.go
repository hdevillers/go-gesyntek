package gesyntek

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/hdevillers/go-gesyntek/kmer"
	"github.com/hdevillers/go-seq/seqio"
	"gonum.org/v1/gonum/mat"
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
	WindowLen      int
	KmerLen        int
	Loci           []Locus
	SeqIdLoci      map[string][]int
	GffTarget      string
	GffId          string
	DistCpt        kmer.KDist
	DistMethod     string
	DistValues     [][]float64
	DistMap        [][]int
	DistDigit      int
	NeedMerge      bool
	IsStandardized bool
}

// Init. GeSynteK object
func NewGeSynteK(w int, k int, t string, i string, m string, d int) *GeSynteK {
	var gsk GeSynteK

	// Initialize arguments
	gsk.WindowLen = w
	gsk.KmerLen = k
	gsk.Loci = make([]Locus, 0)
	gsk.SeqIdLoci = make(map[string][]int)
	gsk.GffTarget = t
	gsk.GffId = i
	gsk.DistMethod = m
	gsk.DistDigit = d
	gsk.NeedMerge = false
	if k > kmer.MaxKSmall {
		gsk.NeedMerge = true
	}
	gsk.IsStandardized = false

	return &gsk
}

// Load data from gff files
func (gsk *GeSynteK) LoadGFF(gff string) error {
	// Create the file handler
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

// Count Kmers in up and down stream sequences
func (gsk *GeSynteK) CountKmers() error {
	for i := range len(gsk.Loci) {
		err := gsk.Loci[i].CountUpDownKmers(gsk.KmerLen)
		if err != nil {
			return err
		}
	}
	return nil
}

// Compute Kmers distance
func (gsk *GeSynteK) ComputeKmerDistance() error {
	// Initialize distance computing class according to method
	switch gsk.DistMethod {
	case "Euclidean":
		gsk.DistCpt = kmer.NewKDistEuclidean()
	case "Cosine":
		gsk.DistCpt = kmer.NewKDistCosine()
	case "Mash":
		gsk.DistCpt = kmer.NewKDistMash(gsk.KmerLen)
	default:
		return errors.New("unsupported kmer distance method")
	}

	// Initialize distance attributes
	nLoci := len(gsk.Loci)
	nDist := (nLoci * (nLoci - 1)) / 2
	iMax := nLoci - 1
	jMin := 1
	if gsk.DistCpt.NeedSelfComparison() {
		nDist = nDist + nLoci
		iMax = nLoci
		jMin = 0
	}
	gsk.DistValues = make([][]float64, nDist)
	gsk.DistMap = make([][]int, nDist)

	// Set up
	z := 0
	for i := 0; i < iMax; i++ {
		for j := i + jMin; j < nLoci; j++ {
			gsk.DistValues[z] = make([]float64, 2)
			gsk.DistMap[z] = make([]int, 2)
			gsk.DistMap[z][0] = i
			gsk.DistMap[z][1] = j
			// Compute up-stream kmer distance
			if gsk.Loci[i].HasUpStr && gsk.Loci[j].HasUpStr {
				err := gsk.DistCpt.Compute(gsk.Loci[i].KmerUpStr.GetCounts(), gsk.Loci[j].KmerUpStr.GetCounts())
				if err != nil {
					return err
				}
				gsk.DistValues[z][0] = gsk.DistCpt.GetDistance()
			} else {
				gsk.DistValues[z][0] = -1
			}
			// Compute down-stream kmer distance
			if gsk.Loci[i].HasDownStr && gsk.Loci[j].HasDownStr {
				err := gsk.DistCpt.Compute(gsk.Loci[i].KmerDownStr.GetCounts(), gsk.Loci[j].KmerDownStr.GetCounts())
				if err != nil {
					return err
				}
				gsk.DistValues[z][1] = gsk.DistCpt.GetDistance()
			} else {
				gsk.DistValues[z][1] = -1
			}
			z++
		}
	}

	return nil
}

// Write out up- and down-stream sequences
func (gsk *GeSynteK) WriteUpDownFasta(ob string) error {
	upOut := seqio.NewWriter(ob+"_UpStream.fasta", "fasta", false)
	doOut := seqio.NewWriter(ob+"_DownStream.fasta", "fasta", false)
	defer upOut.Close()
	defer doOut.Close()

	for i := range len(gsk.Loci) {
		if gsk.Loci[i].HasUpStr {
			upOut.Write(gsk.Loci[i].SeqUpStr)
		}
		if gsk.Loci[i].HasDownStr {
			doOut.Write(gsk.Loci[i].SeqDownStr)
		}
	}

	return nil
}

// Write out distance between each pair of Loci
func (gsk *GeSynteK) WritePairwiseDistance(ob string) error {
	f, err := os.Create(ob + "_Pairwise_" + gsk.DistMethod + ".tsv")
	if err != nil {
		return err
	}
	defer f.Close()

	fw := bufio.NewWriter(f)

	fs := "%.0" + fmt.Sprint(gsk.DistDigit) + "f"

	fw.WriteString("First.Locus\tSecond.Locus\tUpstream.Distance\tDownstream.Distance\n")
	for i := range len(gsk.DistValues) {
		c1 := "NA"
		if gsk.DistValues[i][0] > -0.5 {
			c1 = fmt.Sprintf(fs, gsk.DistValues[i][0])
		}
		c2 := "NA"
		if gsk.DistValues[i][1] > -0.5 {
			c2 = fmt.Sprintf(fs, gsk.DistValues[i][1])
		}
		fmt.Fprintf(fw, "%s\t%s\t%s\t%s\n", gsk.Loci[gsk.DistMap[i][0]].SeqLabel,
			gsk.Loci[gsk.DistMap[i][1]].SeqLabel, c1, c2)
	}
	fw.Flush()

	return nil
}

func (gsk *GeSynteK) WriteKmerCounts(ob string) error {
	// Two files, one for upstream kmers and one for downstream
	fup, err := os.Create(ob + "_UpStream_KmerCounts.tsv")
	if err != nil {
		return err
	}
	defer fup.Close()
	fdo, err := os.Create(ob + "_DownStream_KmerCounts.tsv")
	if err != nil {
		return err
	}
	defer fdo.Close()

	fupw := bufio.NewWriter(fup)
	fdow := bufio.NewWriter(fdo)

	// Create and write the header
	header := "Kmers"
	nLoci := len(gsk.Loci)
	for i := range nLoci {
		header = header + "\t" + gsk.Loci[i].SeqLabel
	}
	fupw.WriteString(header + "\n")
	fdow.WriteString(header + "\n")

	// Convert kmer numerical ids (uint64) into bytes
	kl := kmer.NewKLabel(gsk.KmerLen)
	labUpInt := gsk.Loci[0].KmerUpStr.GetKmers()
	labDownInt := gsk.Loci[0].KmerDownStr.GetKmers()
	nUpKmers := len((*labUpInt)[0])
	nDownKmers := len((*labDownInt)[0])
	labUpByte := make([][]byte, nUpKmers)
	labDownByte := make([][]byte, nDownKmers)
	err = kl.Uint64ToBytes(labUpInt, &labUpByte)
	if err != nil {
		return err
	}
	err = kl.Uint64ToBytes(labDownInt, &labDownByte)
	if err != nil {
		return err
	}

	// Extract and copy count values into an array
	var upVal [][]float64
	var doVal [][]float64
	upVal = make([][]float64, nLoci)
	doVal = make([][]float64, nLoci)
	for i := range nLoci {
		if gsk.Loci[i].HasUpStr {
			upVal[i] = make([]float64, nUpKmers)
			mat.Col(upVal[i], 0, gsk.Loci[i].KmerUpStr.GetCounts())
		} else {
			upVal[i] = make([]float64, 0)
		}

		if gsk.Loci[i].HasDownStr {
			doVal[i] = make([]float64, nDownKmers)
			mat.Col(doVal[i], 0, gsk.Loci[i].KmerDownStr.GetCounts())
		} else {
			doVal[i] = make([]float64, 0)
		}
	}

	// Write count values
	numFmt := "\t%.0f"
	if gsk.IsStandardized {
		numFmt = "\t%.04f"
	}
	for i := range nUpKmers {
		fupw.Write(labUpByte[i])
		for j := range nLoci {
			if gsk.Loci[j].HasUpStr {
				fmt.Fprintf(fupw, numFmt, upVal[j][i])
			} else {
				fupw.WriteString("\tNA")
			}
		}
		fupw.WriteByte('\n')
	}
	for i := range nDownKmers {
		fdow.Write(labDownByte[i])
		for j := range nLoci {
			if gsk.Loci[j].HasDownStr {
				fmt.Fprintf(fdow, numFmt, doVal[j][i])
			} else {
				fdow.WriteString("\tNA")
			}
		}
		fdow.WriteByte('\n')
	}
	fupw.Flush()
	fdow.Flush()

	return nil
}

// Standardize counts
func (gsk *GeSynteK) StandardizeCounts() {
	for i := range len(gsk.Loci) {
		if gsk.Loci[i].HasUpStr {
			kmer.Standardize(gsk.Loci[i].KmerUpStr.GetCounts())
		}
		if gsk.Loci[i].HasDownStr {
			kmer.Standardize(gsk.Loci[i].KmerDownStr.GetCounts())
		}
	}
	gsk.IsStandardized = true
}

// Merge Kmer label for each counts
func (gsk *GeSynteK) MergeKmers() error {
	if gsk.NeedMerge {
		// First compute union of kmer labels
		var UpLab [][]uint64
		var DoLab [][]uint64
		if gsk.KmerLen <= kmer.MaxK64Bits {
			UpLab = make([][]uint64, 1)
			UpLab[0] = make([]uint64, 0)
			DoLab = make([][]uint64, 1)
			DoLab[0] = make([]uint64, 0)
		} else {
			return errors.New("unsupported kmer value (too high) for merging")
		}
		kl := kmer.NewKLabel(gsk.KmerLen)
		for i := range len(gsk.Loci) {
			if gsk.Loci[i].HasUpStr {
				kl.MergeUint64(&UpLab, gsk.Loci[i].KmerUpStr.GetKmers())
			}
			if gsk.Loci[i].HasDownStr {
				kl.MergeUint64(&DoLab, gsk.Loci[i].KmerDownStr.GetKmers())
			}
		}

		// Scan each counts and insert missing label (with zero-count)
		for i := range len(gsk.Loci) {
			if gsk.Loci[i].HasUpStr {
				err := gsk.Loci[i].KmerUpStr.MergeKmers(&UpLab)
				if err != nil {
					return err
				}
			}
			if gsk.Loci[i].HasDownStr {
				err := gsk.Loci[i].KmerDownStr.MergeKmers(&DoLab)
				if err != nil {
					return err
				}
			}
		}
	}

	return nil
}
