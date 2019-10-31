# UST 

UST is a bioinformatics tool for constructing maximal string set representation from sequencing data.


## Usage

Read the instructions below to compile, then:

    ./ust -i [unitigs.fa] -k [kmer_size] 
  
e.g.

    ./ust -i list_reads.unitigs.fa -k 21 

Importants parameters are:

    -k [int] : The k-mer size, i.e. length of the nodes of the de Bruijn graph.
    
    -i [input-file] : Unitigs file produced by become in FASTA format.    
	
    -a [0 or 1] : Default is 0. A value of 1 tells UST to preserve abundance. Use this when bcalm2 output is generated with -all-abundance counts.
    
    
## Pre-requisites:

GCC >= 4.8 or a very recent C++11 capable compiler

## Installation

Compile from source as follows:

    git clone https://github.com/amatur/UST 
    cd UST
    make

## Input formats

File input format can only be BCALM2 unitigs file in fasta format. 
   
## Output

UST outputs the file "stitchedUnitigs.fa" in the working folder, which is the set of strings in maximal path cover of the de Bruijn graph in FASTA format.

