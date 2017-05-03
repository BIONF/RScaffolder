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


## Installation
- There is no installation required. Simply clone the repository (see below). Incompatibilities with your system i.e. paths etc. have to be changed manually in the code. I plan to enhance usability in the next weeks.
- Please check the dependencies section, as RScaffolder requires a couple of external tools.

#### Github
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

#### Python 3.x modules
- Biopython
- Nucmer
- Networkx
- Graphviz

#### MUMmer3.23


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

| Filename | Description |
| --------- | ----------- |
| ./mapped | created in the output directory containing the result files of the mapping with nucmer against all reference genomes  |
| contigs_filtered.fasta | This file contains all contigs after filtering (size or name if specified) in FASTA format |
| contigs_filtered_HT.fasta | Nucleotide FASTA file of the input subsequences (the head and tail) of each contig. |
| gfa_contig_connections.txt | Connected contig ends based on assembly graph information |
| gfa_contigs.paths | Links in assembly graph |
| ht_mapping_distances.tsv | Contais a list of mapped contig end positions for each reference and their distances |
| ht_connectivity_graph.svg | Vector graphic illustration of the connection graph (only non-repetitve contigs) |
| ht_connectivity_graph2.svg | Second vector graphic illustration now incorporating the repetitive contigs |
| contigs_stats.tsv | Statistics with relevant contig information |
| scaffolds.fasta | Final output, contigs are now merged into scaffolds with 100 Ns (standard) as spacer |

## Citation

Not yet published.
