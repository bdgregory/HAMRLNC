# HAMRLINC: High-throughput Analysis of Modified Ribonucleotides and Long Intergenic Non-Coding RNAs
![HAMRbox_Workflow_v5](https://github.com/harrlol/HAMRLINC/assets/87460010/aa2b83f4-889b-42d1-9594-1f6715326ee6)


## Overview
- HAMRLINC is a multipurpose toolbox that expedites the analysis pipeline of two algorithms: [HAMR](https://github.com/GregoryLab/HAMR) and [Evolinc](https://github.com/Evolinc/Evolinc-I/tree/master). The former was developed by [Paul Ryvkin et al](https://rnajournal.cshlp.org/content/19/12/1684), and the latter by [Andrew D.L. Nelson et al](https://www.frontiersin.org/articles/10.3389/fgene.2017.00052/full). HAMRLINC aims to make the original methods more accessible by automating the tedious pre-processing steps, allowing users to analyze RNA-seq data at an experiment scale. 
- HAMRLINC is high-throughput and performs RNA-modification analysis and long intergenic non-coding RNAs(lincRNA) annotation at a bioproject scale. HAMRLINC performs constitutive trimming of acquired reads using Trim-Galore, and makes use of STAR as the aligning tool which reduces runtime, while allowing the users to use Tophat as an option. 

## Command Line Arguments and Description

| Command | Description |
| :---: | :---: |
| Required |
| -o | \<project directory\> <br> where you want your entire hamr project to be |
| -t | \<SRA accession list.txt\> or \<folder of raw fastq files\> <br> a txt file of all srr accession code to your desired reads or a path containing them |
| -c | \<filenames for each fastq.csv\> <br> a csv file that corresponds each srr code (or name of fastq file) to your desired nomenclature for each read |
| -g | \<reference genome.fa> <br> a fasta file of the genome of the model organism |
| -i | \<reference genome annotation.gff3> <br> a gff3 file of the genome of the model organism, note we require gff3 instead of gtf |
| -l | \<read length\> <br> an integer, the read length of this sequencing experiment, if non-unanimous use the shortest length |
| -s | \<genome size in bp\> <br> an integer, the number of base pairs of the genome of this model organism |
| -e | \<genome annotation generator code\> <br> see below for abbreviation code, one code per organism/cultivar |
| Optional |
| -n | \[number of threads\] <br> default=4 |
| -a | \[use Tophat2 instead of STAR\] <br> default uses STAR |
| -b | \[Tophat2 library choice: fr-unstranded, fr-firststrand, fr-secondstrand\] <br> default=fr-firststrand |
| -f | \[filter\] <br> default=filter_SAM_number_hits.pl |
| -p | \[suppress pamlinc\] |
| -u | \[suppress featurecount\] |
| -v | \[evolinc option: M or MO\] <br> default=M |
| -Q | \[HAMR: minimum qualuty score\] <br> default=30 |
| -C | \[HAMR: minimum coverage\] <br> default=50 |
| -E | \[HAMR: sequencing error\] <br> default=0.01 |
| -P | \[HAMR: maximum p-value\] <br> default=1 |
| -F | \[HAMR: maximum FDR\] <br> default=0.05 |
| -T | \</path/to/transposable_elements_file\> <br> Only required under evolinc MO option |
| -G | \</path/to/CAGE_RNA_file\> <br> Only required under evolinc MO option |
| -D | \</path/to/known_lincRNA_file\> <br> Only required under evolinc MO option |
| -m | \[HAMR model\] <br> default=euk_trna_mods.Rdata |
| -h | \[help message\]|

## Annotation Generator Code
| Abbreviation Code | Organism |
| --- | --- |
| AT | Arabidopsis thaliana |
| BD | Brachypodium distachyon |
| ZM | Zea mays |
| OSJ | Oryza sativa japonica |
| OSI | Oryza sativa indica |
| OSIR64 | Oryza sativa IR64 |
