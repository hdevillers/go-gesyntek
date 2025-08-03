package main

import (
	"flag"
	"fmt"

	"github.com/hdevillers/go-gesyntek/kmer"
)

type inputFlag []string

func (i *inputFlag) String() string {
	return fmt.Sprintf("%v", *i)
}

func (i *inputFlag) Set(value string) error {
	*i = append(*i, value)
	return nil
}

var inputs inputFlag

func main() {
	flag.Var(&inputs, "input", "Input sequence file(s).")
	format := flag.String("format", "fasta", "Format of the input sequence file(s).")
	outputBase := flag.String("output-base", "GeSynteK_kmer", "Output base path.")
	kmerLen := flag.Int("kmer-length", 4, "Kmer length.")
	standardize := flag.Bool("standardize", false, "Standardize kmer counts.")
	canonical := flag.Bool("canonical", false, "Count canonical kmers.")
	flag.Parse()

	if len(inputs) == 0 {
		panic("You must provide an input sequence file.")
	}

	km := kmer.NewKmer(*kmerLen, *canonical)

	// Load sequence
	for i := range len(inputs) {
		err := km.LoadSequences(inputs[i], *format)
		if err != nil {
			panic(err)
		}
	}

	// Standardize counts if required
	if *standardize {
		km.StandardizeCounts()
	}

	// Merge kmer
	err := km.MergeKmers()
	if err != nil {
		panic(err)
	}

	// Write out kmer counts
	err = km.WriteKmerCounts(*outputBase)
	if err != nil {
		panic(err)
	}

}
