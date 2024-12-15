# Phylogenetic Tree Constructor (Neighbour-Joining)
A Simple Model Program for Constructing DNA Sequence-base Phylogenetic Tree using Neighbour-Joining Algorithm

> by Naufal Muhammad Alif - 10422039

> supplementary of Tugas Makalah Matematika Diskrit IF1202 2024-2025

## Description
+ This is a program (__NJ_v2.py__) that constructs a phylogenetic tree of specified organsims based on their DNA sequences. 
+ The program <ins>receives a multiple sequence alignment (MSA) file in FASTA format of target organisms' DNA sequences.</ins>
+ Then, the program calculate pairwise distance between the sequences using __Kimura two-distance parameter model__ or __Jukes-Cantor model__, specified by the user upon input.
+ A phylogenetic tree will be constructed using __neighbour-joining algorithm__ and converted into a "Newick tree" format.
+ The Newick tree will be visualized using `Bio.Phylo` and `matplotlib.pyplot` packages

## Samples
+ Two MSA sample files are available (sample_1.fas and sample_2.fas)
+ Partial sequences of 16S rRNA gene retrieved from NCBI GenBank & RefSeq database are used
+ MSA was done through MEGA11 using MUSCLE algorithm
+ sample_1.fas includes:
  - `>NR 042778.1 Streptococcus thermophilus strain ATCC 19258`
  - `>NR 042776.1 Streptococcus salivarius strain ATCC 7073`
  - `>NR 117496.1 Streptococcus pneumoniae strain ATCC 33400`
  - `>NR 117503.1 Streptococcus agalactiae ATCC 13813`
  - `>NR 112088.1 Streptococcus pyogenes strain JCM 5674`
  - `>NC 000913.3:4035531-4037072 Escherichia coli str. K-12 substr.`
+ sample_2.fas includes:
  - `>NR 043424.1 Pseudomonas putida strain IAM 1236`
  - `>NR 114749.1 Pseudomonas protegens strain CHA0`
  - `>FJ972536.1 Pseudomonas fluorescens NO7`
  - `>NR 026078.1 Pseudomonas aeruginosa strain DSM 50071`
  - `>NR 041742.1 Legionella pneumophila subsp. pascullei strain U8W`
  - `>NR 024570.1 Escherichia coli strain U 5/41`
  - `>NR 102794.2 Enterobacter cloacae strain ATCC 13047`
  - `>NR 145647.1 Enterobacter asburiae strain JM-458`

### notes
+ not equipped with input validation nor error-handling capability
+ only accept multiple sequence alignment (MSA) file in FASTA format
+ the MSA must be done separately using other program such as MEGA or ClustalW.
+ FASTA entry identifier must follow the standard NCBI nomenclature, for example
`>NR 043424.1 Pseudomonas putida strain IAM 1236 16S ribosomal RNA partial sequence`
