# RScaffolder: Multi-Reference Assisted Bacterial Genome Scaffolder 

## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Output Files](#output-files)
- [Citation](#citation)

## Introduction
Whole genome sequencing is nowadays a routine tool for biological and clinical research. As a consequence, the number of bacterial genomes in public databases has seen exponential growth over the past two decades. The vast majority of these, however, remain in draft status due to time and/or financial constraints. Although analyses on gene and nucleotide level are generally possible at this stage, the low contiguity of the draft genomes limit their use in high resolution functional and structural comparative genomics. 
Reference-assisted scaffolding harnesses the information of available closely-related reference genomes to improve contiguity subsequent to a de novo assembly step. Unfortunately, choosing a suitable reference is not straight-forward as genome organization can change rapidly during evolution. To circumvent this problem, recent approaches exploit multiple reference genomes to guide the scaffolding process. These programs however, are either are too optimistic and thereby produce many false contig joints or perform poorly with larger numbers of references. Therefore, we developed a new algorithm which efficiently scales with the increasing number of available genomes. Our approach starts by ordering and orienting the contigs along each reference genome. More specifically, flanking regions of each contig are mapped onto all available sequences to derive a set of signed permutations of contigs (with repetitions). We further define two consecutive high quality mappings within a permutation to be “connected” only if their distance is less than the size of the smallest contig length. These subsets of “connected mappings” are then used to resolve ambiguities in the assembly graph. Yet, with regards to possible rearrangement events, we promote scaffolding only in cases where the connection is unambiguous across the reference genomes. In contrast to many existing solutions, our algorithm is not limited by the number or contiguity of the reference genomes and is robust with regards to false positive joints when run with references of greater evolutionary distance. This is achieved not only by the required unambiguity of observed “connections” but also by employing contig coverage information to control multiplicity of the contigs in the resulting scaffolds.
We benchmark our algorihtm with 6 Acinetobacter baumannii genomes assembled both from simulated and real Illumina HiSeq 2500 paired-end read sets together with a set of references comprising 46 publically available finished Acinetobacter genomes. Applying our program, we find that misassembly-corrected contiguity, NGA50, has risen approximately 2 to 20 fold reducing the fragments count by 49% to 76%. Because the work flow is capable of handling large numbers of reference genomes, we further demonstrate how it can be utilised to examine the genomic structure of the assembly by identifying conserved and variable genomic sites across the references taxonomic range.

## Installation
- There is no installation required. Simply clone the repository (see below). Incompatibilities with your system and Paths have to be changed manually in the code. I plan to enhance usability in the next weeks.
- Please check the dependencies section, as RScaffolder requires a couple of external tools.
- Please be aware, that currently this tool

### Github
Choose location to put the repo in, for example in your home directory (no root access required):
```bash
% cd $HOME
```
Clone the latest version of the repository:
```bash
% git clone https://github.com/BIONF/RScaffolder.git
% ls RScaffolder
```

## Dependencies

### Python 3.x modules
- Biopython
- Nucmer
- Networkx
- Graphviz

- MUMmer3.23
- 


## Usage

usage: rscaffolder.py [-h] -i </path/to/contigfile> [-rg <path>]
                      [-spdir <path>] [-o </path/to/out/>]
                      [-DEBUGexclude <str> [<str> ...]] [-f]

Multi-Reference-Genome-Assisted Scaffolding of a target assembly

optional arguments:
  -h, --help            show this help message and exit
  -i </path/to/contigfile>, --contigs_filepath </path/to/contigfile>
                        path to an contigs fasta file
  -rg <path>, --reference-genomes-directory <path>
                        specify a directory containing the reference genomes
                        (multi-) fasta nucleotide sequence files
  -spdir <path>, --assembly-graph-dir <path>
                        specify the spades output directory to parse the
                        assembly graph file in GFA 1.0 format
  -o </path/to/out/>, --outdir </path/to/out/>
                        path to output file prefix
  -DEBUGexclude <str> [<str> ...]
                        References with sequence headers containing any of the
                        specified string will get deleted.
  -f, --force           if set, the files in the specified output directory
                        will be overwritten

## Output Files

| Extension | Description |
| --------- | ----------- |
| .gff | This is the master annotation in GFF3 format, containing both sequences and annotations. It can be viewed directly in Artemis or IGV. |
| .gbk | This is a standard Genbank file derived from the master .gff. If the input to prokka was a multi-FASTA, then this will be a multi-Genbank, with one record for each sequence. |
| .fna | Nucleotide FASTA file of the input contig sequences. |
| .faa | Protein FASTA file of the translated CDS sequences. |
| .ffn | Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA) |
| .sqn | An ASN1 format "Sequin" file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc. |
| .fsa | Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines. |
| .tbl | Feature Table file, used by "tbl2asn" to create the .sqn file. |
| .err | Unacceptable annotations - the NCBI discrepancy report. |
| .log | Contains all the output that Prokka produced during its run. This is a record of what settings you used, even if the --quiet option was enabled. |
| .txt | Statistics relating to the annotated features found. |
| .tsv | Tab-separated file of all features: locus_tag,ftype,gene,EC_number,product |

## Citation

Not yet published.
