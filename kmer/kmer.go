package kmer

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/hdevillers/go-seq/seqio"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/stat"
)

func Standardize(a *mat.Dense) {
	// Extract count values into a float64
	cnt := mat.Col(nil, 0, a)

	// Get mean and sd
	mean, sd := stat.MeanStdDev(cnt, nil)

	for i := range len(cnt) {
		cnt[i] = (cnt[i] - mean) / sd
	}

	//a.Reset()
	a.SetCol(0, cnt)
}

type Kmer struct {
	K       int
	Counter []KCount
	Labels  []string
	Dist    KDist
	IsStd   bool
}

func NewKmer(k int) *Kmer {
	var km Kmer
	km.K = k
	km.Counter = make([]KCount, 0)
	km.Labels = make([]string, 0)
	km.IsStd = false
	return &km
}

func (km *Kmer) LoadFasta(f string) error {
	// Open the fasta file
	seqIn := seqio.NewReader(f, "fasta", false)
	seqIn.CheckPanic()
	defer seqIn.Close()

	// Add a new counter
	ic := len(km.Counter)
	if km.K <= MaxKSmall {
		km.Counter = append(km.Counter, NewKCountSmall(km.K))
	} else if km.K <= MaxK64Bits {
		km.Counter = append(km.Counter, NewKCount31(km.K))
	} else {
		return errors.New("value of K is too high (maximal supported value is 31)")
	}

	// Count kmers
	for seqIn.Next() {
		seqIn.CheckPanic()
		seq := seqIn.Seq()
		err := km.Counter[ic].Count(&seq.Sequence)
		if err != nil {
			return err
		}
	}

	// Add a label from the fasta file
	lab := filepath.Base(f)
	lab, _ = strings.CutSuffix(lab, filepath.Ext(lab))
	km.Labels = append(km.Labels, lab)

	return nil
}

func (km *Kmer) StandardizeCounts() {
	for i := range len(km.Counter) {
		Standardize(km.Counter[i].GetCounts())
	}
	km.IsStd = true
}

// Merge kmer labels for each count
func (km *Kmer) MergeKmers() error {
	if km.K > MaxKSmall {
		if len(km.Counter) > 1 {
			// Compute union of all kmer labels
			var kunion [][]uint64
			if km.K <= MaxK64Bits {
				kunion = make([][]uint64, 1)
				kunion[0] = make([]uint64, 0)
			} else {
				return errors.New("value of K is too high (maximal supported value is 31)")
			}
			kl := NewKLabel(km.K)
			for i := range len(km.Counter) {
				kl.MergeUint64(&kunion, km.Counter[i].GetKmers())
			}

			// Scan each counter and update kmer labels
			for i := range len(km.Counter) {
				err := km.Counter[i].MergeKmers(&kunion)
				if err != nil {
					return err
				}
			}
		}

	}
	return nil
}

// Write Kmer counts
func (km *Kmer) WriteKmerCounts(ob string) error {
	f, err := os.Create(ob + "_KmerCounts.tsv")
	if err != nil {
		return err
	}
	defer f.Close()

	// Create a buffer
	fw := bufio.NewWriter(f)

	// Create and write the header
	header := "Kmers"
	nSeq := len(km.Counter)
	for i := range nSeq {
		header = header + "\t" + km.Labels[i]
	}
	fw.WriteString(header + "\n")

	// Convert kmer uint64 labels into bytes
	kl := NewKLabel(km.K)
	kNum := km.Counter[0].GetKmers()
	nKmers := len((*kNum)[0])
	kByte := make([][]byte, nKmers)
	err = kl.Uint64ToBytes(kNum, &kByte)
	if err != nil {
		return err
	}

	// Extract count values
	cnt := make([][]float64, nSeq)
	for i := range nSeq {
		cnt[i] = make([]float64, nKmers)
		mat.Col(cnt[i], 0, km.Counter[i].GetCounts())
	}

	// Write count values
	numFmt := "\t%.0f"
	if km.IsStd {
		numFmt = "\t%.04f"
	}
	for i := range nKmers {
		fw.Write(kByte[i])
		for j := range nSeq {
			fmt.Fprintf(fw, numFmt, cnt[j][i])
		}
		fw.WriteByte('\t')
	}
	fw.Flush()

	return nil
}
