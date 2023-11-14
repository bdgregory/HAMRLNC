# HAMRLINC: High-throughput Analysis of Modified Ribonucleotides and Long Intergenic Non-Coding RNAs
![HAMRLINC_Workflow_v4](https://github.com/harrlol/HAMRLINC/assets/87460010/d34fc1f2-2c5d-4b41-98ea-7a97f35fe4c3)



## Overview
- HAMRLINC is a multipurpose toolbox that expedites the analysis pipeline of two algorithms: [HAMR](https://github.com/GregoryLab/HAMR) and [Evolinc](https://github.com/Evolinc/Evolinc-I/tree/master). The former was developed by [Paul Ryvkin et al](https://rnajournal.cshlp.org/content/19/12/1684), and the latter by [Andrew D.L. Nelson et al](https://www.frontiersin.org/articles/10.3389/fgene.2017.00052/full). HAMRLINC aims to make the original methods more accessible by automating the tedious pre-processing steps, and expand on their functionalities with its built-in post-processing steps, allowing users to perform RNA modification prediction with intuitive output formats.
- HAMRLINC is high-throughput and performs RNA-modification analysis and long intergenic non-coding RNAs(lincRNA) annotation at a bioproject scale. HAMRLINC performs constitutive trimming of acquired reads using Trim-Galore, and makes use of STAR (Tophat option available) as the default aligning tool; mapped reads are pre-processed using selected methods from GATK and samtools.
- HAMRLINC is optimized for partial parallel processing, and modularization. Specifying a larger thread count where hardware permits will greatly increase the efficiency of a run. If only partial functionality is needed (e.g. Only analyzing modified ribonucleotides), users can implement flags to activate  function modules desired. See below for more details. 

## Command Line Arguments and Description

| Command | Description |
| :---: | :---: |
| Required |
| -o | \<project directory\> <br> where you want your entire hamr project to be |
| -c | \<filenames for each fastq.csv\> <br> a csv file that corresponds each srr code (or name of fastq file) to your desired nomenclature for each read |
| -g | \<reference genome.fa> <br> a fasta file of the genome of the model organism |
| -i | \<reference genome annotation.gff3> <br> a gff3 file of the genome of the model organism, note we require gff3 instead of gtf |
| -l | \<read length\> <br> an integer, the read length of this sequencing experiment, if non-unanimous use the shortest length |
| -s | \<genome size in bp\> <br> an integer, the number of base pairs of the genome of this model organism |
| Optional |
| -n | \[number of threads\] <br> default=4 |
| -d | \[raw fastq folder\] <br> a path to a folder containing raw fastq files if needed, in such case, -c csv should have each fastq file as key
| -a | \[use Tophat2 instead of STAR\] <br> default uses STAR |
| -b | \[Tophat2 library choice: fr-unstranded, fr-firststrand, fr-secondstrand\] <br> default=fr-firststrand |
| -f | \[filter\] <br> default=filter_SAM_number_hits.pl |
| -k | \[activate modification analysis (left arm)\] |
| -p | \[activate lincRNA identification (inner right arm)\] |
| -u | \[activate regular featurecount (outer right arm)\] |
| -v | \[evolinc option: M or MO\] <br> default=M |
| -Q | \[HAMR minimum qualuty score: 0-40\] <br> default=30 |
| -C | \[HAMR minimum coverage: 0-∞\] <br> default=50 |
| -E | \[HAMR sequencing error: 0-1\] <br> default=0.01 |
| -P | \[HAMR maximum p-value: 0-1\] <br> default=1 |
| -F | \[HAMR maximum FDR: 0-1\] <br> default=0.05 |
| -O | \[Panther [organism taxon ID](http://pantherdb.org/services/oai/pantherdb/supportedgenomes)\] <br> default="3702" |
| -A | \[Panther [annotation dataset](http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets)\] <br> default="GO:0008150" |
| -Y | \[Panther test type: FISHER or BINOMIAL\] <br> default="FISHER" |
| -R | \[Panther correction type: FDR, BONFERRONI, or NONE\] <br> default="FDR" |
| -T | \</path/to/transposable_elements_file\> <br> Only required under evolinc MO option |
| -G | \</path/to/CAGE_RNA_file\> <br> Only required under evolinc MO option |
| -D | \</path/to/known_lincRNA_file\> <br> Only required under evolinc MO option |
| -m | \[HAMR model\] <br> default=euk_trna_mods.Rdata |
| -h | \[help message\]|
