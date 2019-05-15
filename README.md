# Read Mapping for High-Throughput Library Screening
## Introduction
When using CRISPR library for genetic screening, guide RNA sequences integrated into the pool of cells can be amplified through PCR and then be sequenced. The designed guide sequences usually takes only part of the read. This script is designed to extract parts of read sequence and map them to designed libraries. 

## Guidance
The main script is [get_library_index.py](scripts/get_library_index.py)
### Prerequisite
The program is developed in **Python 3.7**. Dependent packages include:
* argparse
* pandas
* numpy
* re

### Input files
#### Backbone sequence
Suppose the sequenced region looks like this:
`const1 --- alt1 --- const2 --- alt2 --- const3`, just split backbone into rows containing parts of constant sequences, like this:
```
const1
const2
const3
```
Note that the direction of backbone sequence should be same as the direction of library file. 

#### Library files
Librirary files must be given by `csv` format, with columns separated by a comma (,). The first colmn is the target gene ID, the sdecond column is the sgRNA ID. From the third column to the last, each column is a pool of designed sequence corresponding to the alternative sequences within reads.
If there is only one library file containing multiple sequence columns, these sequence are combined according to the index (each row appear in the same read). Otherwise, if there are more than one library files, these sequences are randomly combined. 

#### Reads files
One or two fasta/fastq files are provided according to whether the sequencing is single-end or paired-end. 

### Output files
Output file is a single `csv` files with two columns, ID and counts. As previously described, if the sequences are not randomly combined, this ID is just the same as provided in the library file. If multiple library files are provided, the ID is a combination of IDs from different libraries linked by `_`.

The script will also output status of mapping and running time to standard output. 

### Basic usage
For single-end sequencing, 
```
$ python3 get_library_index.py -i inputfile -l libraries -b backbone
```

For paired-end sequencing,
```
$ python3 get_library_index.py -i inputfile1 -i inputfile2 -l libraries -b backbone -p
```

For more parameters, use:
```
$ python3 get_library_index.py -h


usage: get_library_index.py [-h] [--input INPUT] [--library LIBRARY]
                            [--backbone BACKBONE]
                            [--outputprefix OUTPUTPREFIX] [--paired]
                            [--undirectional] [--seed SEED]
                            [--matchnum MATCHNUM]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        input files
  --library LIBRARY, -l LIBRARY
                        library files csv
  --backbone BACKBONE, -b BACKBONE
                        backbone sequence file
  --outputprefix OUTPUTPREFIX, -o OUTPUTPREFIX
                        output file prefix
  --paired, -p          paired-end reads
  --undirectional, -ud  reads are not of same direction
  --seed SEED, -sd SEED
                        seed length for constructing backbone and mapping
  --matchnum MATCHNUM, -mn MATCHNUM
                        minimum match number
```
`-sd` and `-mn` need to be adjusted according to sequence quality and purpose. The default value is 4 and 3 respectively. See supplementary file [mapping rate](supplementary/S1.csv) for reference. 

### Example
Sample input and output files are in [data](data) directory.
```
$ python3 get_library_index.py -i single_2_ibar.fasta -l ibar_lib_100.csv -b ibar_back.txt
```
```
$ python3 get_library_index.py -i single_geckov2.fasta -l geckov2_100.csv -b backbone.txt
```

## Other
[random_seq.py](scripts/random_seq.py) helps with generating random read sequences from given library. Number of mutations in backbone can be specified. 
```
$ python3 random_seq.py -h

usage: random_seq.py [-h] [--library LIBRARY] [--backbone BACKBONE] [--n N]
                     [--length LENGTH] [--output OUTPUT] [--lenrange LENRANGE]
                     [--indel INDEL] [--mut MUT]

optional arguments:
  -h, --help            show this help message and exit
  --library LIBRARY, -l LIBRARY
                        library files csv
  --backbone BACKBONE, -b BACKBONE
                        backbone sequence file
  --n N, -n N           amounts of sequences
  --length LENGTH, -len LENGTH
                        length of backbone sequence
  --output OUTPUT, -o OUTPUT
                        output file: .fa, .fasta, .fq, .fastq
  --lenrange LENRANGE, -rg LENRANGE
                        length range of backbone sequence: len +- rg
  --indel INDEL, -indel INDEL
                        indel
  --mut MUT, -mut MUT   single nucleotide mutation
```

One example:
```
$ python3 random_seq.py -l geckov2_100.csv -b backbone.txt -o single_geckov2.fasta -len 16 -rg 2 -n 1000 -indel 0 -mut 0
```

