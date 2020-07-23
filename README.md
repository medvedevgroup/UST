# UST

UST is a bioinformatics tool for constructing a spectrum-preserving string set (SPSS) representation from sets of k-mers.

__Note__: This software has been subsumed by [ESSCompress](https://github.com/medvedevgroup/ESSCompress/). To use UST, download ESSCompress and follow the UST instructions in the [README](https://github.com/medvedevgroup/ESSCompress/blob/master/README.md#Running-in-UST-mode).

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

*  `k [int]` : The k-mer size that was used to generate the input, i.e. the length of the nodes of the node-centric de Bruijn graph.
*  `i [input-file]` : Unitigs file produced by [BCALM2 in FASTA format](https://github.com/GATB/bcalm#output).
*  `a [0 or 1]` : Default is 0. A value of 1 tells UST to preserve abundance. Use this option when the input file was generated with the  `-all-abundance counts` option of BCALM2.

The output is a FASTA file with extenstion "ust.fa" in the working folder, which is the SPSS representaiton of the input. If the program is run with the option -a 1, an additional count file with extension "ust.counts" will also be generated.


## Detailed Usage

In order to build a SPSS representation for your k-mer set, you must first run [BCALM2](https://github.com/GATB/bcalm) on your set of k-mers. BCALM2 will construct a set of unitigs. Those unitigs are then fed as input to `ust`, which outputs a FASTA file with the SPSS representation. Note that the k parameter to `ust` must match the `-kmer-size` used when running BCALM2.

If you would like to store the data on disk in compressed form (like UST-Compress in our paper), you can then install and run [MFCompress](http://bioinformatics.ua.pt/software/mfcompress/) on the output of UST as follows: `MFCompressC mykmers.ust.fa`

If you would like to build a membership data structure based on UST, then
- Install [bwtdisk](http://people.unipmn.it/manzini/bwtdisk/) and [dbgfm](https://github.com/jts/dbgfm).
- Change the two variables "DBGFM_DIRECTORY" and "BWTDISK_DIRECTORY" in the script `ust-fm.sh` to point to the locations where dbgfm and bwtdisk are installed. Alternatively, you can add the path to both tools in your environment PATH variable and then modify the script accordingly.
- Run `ust-fm.sh` as follows: `ust-fm.sh mykmers.ust.fa`

## Citation

If using UST in your research, please cite
* Amatur Rahman and Paul Medvedev, [Representation of k-mer sets using spectrum-preserving string sets](http://doi.org/10.1007/978-3-030-45257-5\_10), RECOMB 2020.
* Here is the bibtex entry:

```
@inproceedings{RahmanMedvedevRECOMB20,
  author    = {Amatur Rahman and Paul Medvedev},
  title     = {Representation of $k$-mer sets using spectrum-preserving string sets},
  booktitle = {Research in Computational Molecular Biology - 24th Annual International Conference, {RECOMB} 2020, Padua, Italy, May 10-13, 2020, Proceedings},
  series    = {Lecture Notes in Computer Science},
  volume    = {12074},
  pages     = {152--168},
  publisher = {Springer},
  year      = {2020
}
```

Note that the general notion of an SPSS was independently introduced under the name of simplitigs. Therefore, if citing this general notion, please also cite:
* Brinda K, Baym M, and Kucherov G, [Simplitigs as an efficient and scalable representation of de Bruijn graphs](https://doi.org/10.1101/2020.01.12.903443), bioRxiv 2020.




