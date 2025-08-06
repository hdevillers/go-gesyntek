# go-GeSynteK

`go-GeSynteK` for Gene Synteny with Kmers is a tool that evaluate synteny around orthologous copies of a gene between different strains/species. It only requires a `Fasta` file with the chromosomes of the compares species and a `GFF3` files that provide the coordinate of the gene of interest on these chromosomes. It is written in `Go`, it can count *canonical* or *regular* kmers up to 64 nucleotides and it is quite fast.

## Install

### From sources

To build `go-GeSynteK` from source, you will need to install `Go`. See instructions [here](https://go.dev/doc/install).

Get sources from github, compile, test and install:

```{bash}
# Clone
git clone https://github.com/hdevillers/go-gesyntek.git

# Enter the cloned directory
cd go-gesyntek

# Compile
make

# Test the code
make test

# Install (by default in /usr/local/bin => required root mode)
sudo make install
```

### From `conda`

A `conda` will be available soon.

## Quick start

`go-GeSynteK`needs only two input files:

* A `Fasta` file with the different chromosomes/scaffolds you want to investigate.
* A `GFF3` file with the coordinates of your favorite gene on these different chromosomes.

To run the analysis with kmers of 9 nucleotides and using a Mash-like distance to compare kmers frequencies:

```{bash}
gesyntek-run -gff data.gff -fasta data.fasta -kmer-lengh 9 \
    -dist-method Mash -output-base out
```

The tool comes with a utility to draw an heatmap from the computed distances:

```{bash}
gesyntek-heatmap -input out_Pairwise_Mash.tsv -output heatmap.pdf
```

It will generate a `PDF` file containing the heatmap.

To get the full list of options/arguments of these two commands, you can display the help section:

```{bash}
gesyntek-run -h
gesyntek-heatmap -h
```

Last, `go-GeSynteK` comes with an additional command that simply compute kmer frequencies from `Fasta` or `Fastq` file(s).

```{bash}
# Only one input
kmer-count -input genome1.fasta -output-base out -kmer-length 25

# Multiple inputs
kmer-count -input genome1.fasta -input genome2.fasta \
    -input genome3.fasta -intput genome4.fasta       \ 
    -output-base out -kmer-length 25

# Count canonical kmers
kmer-count -input genome1.fasta -output-base out     \
    -canonical -kmer-length 25
```

Maximal kmer length that our tool can consider is **64** nucleotides.
