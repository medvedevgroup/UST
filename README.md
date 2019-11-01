# UST 

UST is a bioinformatics tool for constructing minimum weight string set representation from sets of k-mers.    
    
## Requirements

GCC >= 4.8 or a C++11 capable compiler

## Installation

Compile from source:

    git clone https://github.com/medvedevgroup/UST 
    cd UST
    make
    
## Usage

After compiling, use

    ./ust -i [unitigs.fa] -k [kmer_size] 
  
e.g.

    ./ust -i examples/k11.unitigs.fa -k 11 

Importants parameters are:

    -k [int] : The k-mer size, i.e. length of the nodes of the de Bruijn graph. You need to use the exact same k value used in bcalm 2.
    
    -i [input-file] : Unitigs file produced by become in FASTA format.    
	
    -a [0 or 1] : Default is 0. A value of 1 tells UST to preserve abundance. Use this when bcalm2 output is generated with -all-abundance counts.


## Input formats

File input format can only be BCALM2 unitigs file in fasta format. 
   
## Output

UST outputs the file "stitchedUnitigs.fa" in the working folder, which is the set of strings in maximal path cover of the de Bruijn graph in FASTA format.

