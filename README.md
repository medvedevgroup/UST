# UST 

UST is a bioinformatics tool for constructing a spectrum-preserving string set (SPSS) representation from sets of k-mers.
    
## Requirements

GCC >= 4.8 or a C++11 capable compiler

    
## Quick start

To install, compile from source:

    git clone https://github.com/medvedevgroup/UST 
    cd UST
    make

After compiling, use

    ./ust -i [unitigs.fa] -k [kmer_size] 
  
e.g.

    ./ust -i examples/k11.unitigs.fa -k 11 


The important parameters are:

    -k [int] : The k-mer size that was used to generate the input, i.e. the length of the nodes of the node-centric de Bruijn graph. 
    
    -i [input-file] : Unitigs file produced by [BCALM2 in FASTA format](https://github.com/GATB/bcalm#output). 
	
    -a [0 or 1] : Default is 0. A value of 1 tells UST to preserve abundance. Use this option when the input file was generated with the  `-all-abundance counts` option of BCALM2.

The output is a FASTA file "stitchedUnitigs.fa" in the working folder, which is the SPSS representaiton of the input.


## Detailed Usage

In order to build a SPSS representation for your k-mer set, you must first run [BCALM2](https://github.com/GATB/bcalm) on your set of k-mers. BCALM2 will construct a set of unitigs. Those unitigs are then fed as input to `ust`, which outputs a FASTA file with the SPSS representation. Note that the k parameter to `ust` must match the `-kmer-size` used when running BCALM2.

If you would like to store the data on disk in compressed form (like UST-Compress in our paper), you can then run [MFCompress](http://bioinformatics.ua.pt/software/mfcompress/) on the output of UST as follows: `MFCompressC mykmers.ust.fa`

If you would like to build a membership data structure based on UST, then 
- Install [bwtdisk](http://people.unipmn.it/manzini/bwtdisk/) and [dbgfm](https://github.com/jts/dbgfm). 
- Make sure the path to both tools is in your environment PATH variable.
- Run `ust-fm.sh`

## Citation

If using UST in your research, please cite 
* Amatur Rahman and Paul Medvedev, "Representation of k-mer sets using spectrum-preserving string sets", in submission.

   
