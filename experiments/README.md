# Experiments 

In this page, we describe all the necessary steps to reproduce the results in the paper.


## File Descriptions

cobs.txt: contains the ID of the files for the cortex data from the paper [COBS](https://arxiv.org/abs/1905.09624).

datasets.txt: contains the commands for filtering SRA data.

## Commands for other tools

#### SRA fastq-dump
For human gut metagenome data used in paper:

`fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR341725`

This will generate two "fastq.gz" files.

#### BCALM 2


##### Version to download
- If you want support for counts, you need to download the latest commit (as of October 31, 2019), or commit with hash f4e0012e8056c56a04c7b00a927c260d5dbd2636.  
- Otherwise, release v2.2.1 or earlier releases will also work.

##### Prepare Input 
Make a file (i.e. "list_reads") with the path of the two "fastq.gz" files. 

##### Run bcalm 2 without counts
`bcalm -in list_reads -kmer-size 31 -abundance-min 2 -max-memory 1000`

This creates a file name "list_reads.unitigs.fa" in the working directory. The output header is:

	><id> LN:i:<length> KC:i:<abundance> KM:f:<abundance> L:<+/->:<other id>:<+/-> [..]
	
	
##### Run bcalm 2 with counts
`bcalm -in list_reads -kmer-size 31 -abundance-min 2 -max-memory 1000 -all-abundance-counts`

The output header looks like this:

	>0 LN:i:80 ab:Z:2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
	GAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTG
	>1 LN:i:32 ab:Z:2 2
	GACACATGCAGCTCCCGGAGACGGTCACAGCT
	
#### dbgfm
- Install [bwtdisk](http://people.unipmn.it/manzini/bwtdisk/) and [dbgfm](https://github.com/jts/dbgfm). 
- Run `ust-fm.sh`. This assumes you have the path to both tools added in your environment.

#### MFCompress
- Download [MFCompress](http://bioinformatics.ua.pt/software/mfcompress/)
- To compress: `MFCompressC stitchedUnitigs.fa`
- To decompress: `MFCompressD stitchedUnitigs.fa.mfc`

#### KMC 3
`kmc -k31 -ci2  -v  @list_reads /home/user/kmc31.2_outdir /home/user/kmc31.2_outdir/tmp`

#### Cosmo/VARI
- Follow the instructions to install [VARI](https://github.com/cosmo-team/cosmo/tree/VARI). Last checked with commit d35bc3dd2d6ba7861232c49274dc6c63320cedc1.
- Run kmc. Make a file `list_kmc` with the path to the kmc outputs (i.e. `/home/user/kmc31.2_outdir/kmc31.2`) .
- Run `cosmo-build -d list_kmc`
- This outputs the binary succinct DBG file `list_kmc.dbg`.

#### Squeakr-exact
`squeakr count -e -k 31 -c 2 -s 2000 -t 1  -o output_filename SRR341725_pass_1.fastq SRR341725_pass_2.fastq`

#### DSK
`dsk  -file list_reads -kmer-size 31 -abundance-min 2`

#### McCortex
- Compile McCortex with k = 31.
- To get FASTA files with unitigs `mccortex31 unitigs DRR000001.ctx > DRR000001.fa`   

###### To generate the compressed UST format:
- run bcalm2 on DRR000001.fa
- run ust on DRR000001.unitgs.fa
- run MFCompressC on stitchedUnitigs.fa	
	- Final compressed output is stitchedUnitigs.fa.mfc


## Command to measure time and memory usage
`/usr/bin/time`
