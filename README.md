# Xmuse
miXture Models for nUcleic acid SequencE analysis. Our paper about using mixture models to analyze sequence features is now [on bioXriv]()!

## Prerequisite
- python >=3.7
- sci-kit learn >=0.22
- seaborn >=0.12.2
- numpy>=1.17.4
- pandas>=1.3.5


## Getting Started
1. Clone or download this repository.
### Using Xmuse as a package
2. Run jupyter notebook
3. Open Analyze human intron sequence features.ipynb
### Using Xmuse as a command line tool
Check out the `example_info.md` in the `example` folder.
3. `cd example`
4. `python3 ../xmuse.py -if intron_vs_exon_ASC.fa -it intron_vs_exon_label.txt -m ASC -k 4 -t 4 -sw 4 -o ASC_result/`
5. More information can be found [here](example/example_info.md)

## Overview
Xmuse is a flexible and robust tool to study biological sequence features using mixture models/latent Dirichlet allocation. Xmuse is adaptable for multiple types of sequence feature searching tasks:
- Identify features associated with positions using groups of aligned sequences through evaluating the k-mer compositions across positions
- Summarize repetitive features in sequences and identify subgroups from sequences of various lengths
- Summarize positional features in sequences and identify subgroups from aligned sequences of the same length

## Get Xmuse
Clone or download this repository

## Usage
Xmuse can be used as either a command line tool or a python package.

### Command line options
python3 xmuse.py -if input.fa -t tags.txt -m ASC -k 4 -t 4 -o output_folder

#### Required parameters
- `-if`: The path to the input fasta file(s). If multiple fasta files need to be provided, use ',' to separate the paths
- `it`: Two options are allowed here. (1) Provide the path to a txt file in which each line corresponds to each sequence in the input fasta file. This option should be used when single fasta file is provided by `if`. (2) Provide a string as the tag for a bunch of sequences. If multiple fasta files are provided by `if` and each one corresponds to one tag, use ',' to separate the tags.
- `o`: The path to the output directory.

#### Optional parameters

- `m`: The mode to run LDA. Provide one of the following options. (1) ASC, aligned sequence characterization. This mode will summarize k-mers at positions and model the features across positions. Sequences provided must have the same length. (2) AWSC, aligned whole sequence characterization. This model will summarize positional k-mers in sequences or groups of sequences. Sequences provided must have the same length. (3) UWSC, unaligned whole sequence feature characterization. This model will summarize k-mers in sequences or groups of sequences. Default: UWSC
- `k`: The integer as the size of k-mers to be used in sequence feature modeling. Default: 4
- `t`: The interger as the number of topics. Default: 4
- `sw`: The integer as the size of the sliding window to be allied in k-mer counting. Specifically used in ASC mode. Default: 0
- `sg`: The integer as the size of subgroup size. Sequences will be randomly put into subgroups to increase the counts of k-mers in each sample, which may increase the sequence modeling performance. Specifically used in UWSC and AWSC modes. Default: 1
- `s`: The integer as the position of the first nucleotide in a sequence. Specifically used in ASC and AWSC modes. Default: 0
- `fm`: The path to a pickle file of a muse object. Use this parameter when the task is to transform new sequences into topic distributions. If not provided, new model will be calculated.

### Xmuse library
Check out the jupyter notebook in the example directory!