# bsponge_chloroplast

Supplement scripts and data files for the manuscript entitled "The signs of adaptive mutations identified in the chloroplast genome of the algae endosymbiont of Baikal sponge."

## Description

The scripts on python intended to generage figures in the manuscript use several data files, provided with the scripts

python version 2.7 was used.

The data files consists of:

- Choloroplast genome fragment, discussed in the manuscript, as a sequence file
- Sequencing reads from five samples used in the study, aligned to the sequence of the genome fragment, in .bam format
- Complete chloroplast genomes of several algae species in genbank format
- Phylogenetic trees in newick format

The python libraries required to run the scripts:

- biopython (1.66)
- pysam (0.14.1)
- matplotlib (2.2.2)

For the pie chart with counts of polymorphic sites for each gene (fig. 3 in the manuscript), the log file is provided with numeric values of the distributions presented in the chart (polymorph_distr.out) 

mafft, fastme and readseq software is required to re-calculate phylogenetic trees

In addition, the scripts which could be used to assemble genome fragment from raw sequencing reads and a sequence of reference genome are provided.
But the raw sequencing archives are not deposited yet, this data is used in several other projects and the submission is postponed a while. 

The software was used on Ubuntu Linux x86_64 server.




