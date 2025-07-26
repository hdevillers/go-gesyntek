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
	kmerLen := flag.Int("kmer-length", gesyntek.KMER_LEN, "Kmer length to consider.")
	windowLen := flag.Int("window-length", gesyntek.WINDOW_LEN, "Window length around loci.")

	flag.Parse()

	if *gff == "" {
		panic("You must provide an input GFF file.")
	}

	// Initialize the structure
	gsk := gesyntek.NewGeSynteK(*windowLen, *kmerLen, *gffTarget, *gffID)

	// Load the GFF file
	err := gsk.LoadGFF(*gff)
	if err != nil {
		panic(err)
	}

	for i := 0; i < len(gsk.Loci); i++ {
		fmt.Println(gsk.Loci[i].SeqId)
		fmt.Println(gsk.Loci[i].SeqStart)
		fmt.Println(gsk.Loci[i].SeqEnd)
		fmt.Println(gsk.Loci[i].SeqStrand)
		fmt.Println(gsk.Loci[i].SeqLabel)
		fmt.Println(gsk.SeqIdLoci[gsk.Loci[i].SeqId])
	}
}
