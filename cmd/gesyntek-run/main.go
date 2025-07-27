package main

import (
	"flag"
	"fmt"

	"github.com/hdevillers/go-gesyntek/gesyntek"
)

func main() {
	gff := flag.String("gff", "", "GFF input file (loci description).")
	gffTarget := flag.String("gff-target", gesyntek.GFF_TARGET, "GFF feature type to target.")
	gffID := flag.String("gff-id", gesyntek.GFF_ID, "GFF description flag to define the locus ID.")
	fasta := flag.String("fasta", "", "Input Fasta file(s).")
	kmerLen := flag.Int("kmer-length", gesyntek.KMER_LEN, "Kmer length to consider.")
	windowLen := flag.Int("window-length", gesyntek.WINDOW_LEN, "Window length around loci.")

	flag.Parse()

	if *gff == "" {
		panic("You must provide an input GFF file.")
	}

	if *fasta == "" {
		panic("You must provide an input Fasta file.")
	}

	// Initialize the structure
	gsk := gesyntek.NewGeSynteK(*windowLen, *kmerLen, *gffTarget, *gffID)

	// Load the GFF file
	err := gsk.LoadGFF(*gff)
	if err != nil {
		panic(err)
	}

	// Extract the up and down stream sequences of each locus
	err = gsk.LoadFasta(*fasta)
	if err != nil {
		panic(err)
	}

	// Count Kmers in up and down stream sequences
	err = gsk.CountKmers()
	if err != nil {
		panic(err)
	}

	tmp := gsk.Loci[0].KmerUpStr.GetCounts()
	for i := range 64 {
		fmt.Println((*tmp)[i])
	}
}
