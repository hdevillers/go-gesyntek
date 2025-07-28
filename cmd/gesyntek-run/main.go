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
	distMethod := flag.String("dist-method", "Euclidean", "Kmer distance method.")

	flag.Parse()

	if *gff == "" {
		panic("You must provide an input GFF file.")
	}

	if *fasta == "" {
		panic("You must provide an input Fasta file.")
	}

	// Initialize the structure
	gsk := gesyntek.NewGeSynteK(*windowLen, *kmerLen, *gffTarget, *gffID, *distMethod)

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

	// Compute all pairwise distances
	err = gsk.ComputeKmerDistance()
	if err != nil {
		panic(err)
	}

	for i := range len(gsk.DistValues) {
		fmt.Printf("Comparison of %s and %s:\n", gsk.Loci[gsk.DistMap[i][0]].SeqLabel, gsk.Loci[gsk.DistMap[i][1]].SeqLabel)
		fmt.Printf("    Up-stream distance: %.04f\n", gsk.DistValues[i][0])
		fmt.Printf("    Down-stream distance: %.04f\n", gsk.DistValues[i][1])
	}

	/*for i := range 64 {
		fmt.Printf("Comptage Up G2 G5: %d\t%d\n", (*gsk.Loci[1].KmerUpStr.GetCounts())[i], (*gsk.Loci[4].KmerUpStr.GetCounts())[i])
		fmt.Printf("Comptage Down G2 G5: %d\t%d\n", (*gsk.Loci[1].KmerDownStr.GetCounts())[i], (*gsk.Loci[4].KmerDownStr.GetCounts())[i])
	}*/
}
