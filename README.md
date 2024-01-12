# HAMRLINC: High-throughput Analysis of Modified Ribonucleotides and Long Intergenic Non-Coding RNAs
![HAMRLINC_Workflow_v4](https://github.com/harrlol/HAMRLINC/assets/87460010/d34fc1f2-2c5d-4b41-98ea-7a97f35fe4c3)



## Overview
- HAMRLINC is a multipurpose toolbox that expedites the analysis pipeline of two algorithms: [HAMR](https://github.com/GregoryLab/HAMR) and [Evolinc](https://github.com/Evolinc/Evolinc-I/tree/master). The former was developed by [Paul Ryvkin et al](https://rnajournal.cshlp.org/content/19/12/1684), and the latter by [Andrew D.L. Nelson et al](https://www.frontiersin.org/articles/10.3389/fgene.2017.00052/full). HAMRLINC aims to make the original methods more accessible by automating the tedious pre-processing steps and expanding on their functionalities with its built-in post-processing steps, allowing users to perform RNA modification prediction with intuitive output formats.
- HAMRLINC is high-throughput and performs RNA-modification analysis and long intergenic non-coding RNAs(lincRNA) annotation at a bioproject scale. HAMRLINC performs constitutive trimming of acquired reads using Trim-Galore, and makes use of STAR (Tophat option available) as the default aligning tool; mapped reads are pre-processed using selected methods from GATK and samtools.
- HAMRLINC is optimized for partial parallel processing and modularization. Specifying a larger thread count where hardware permits will greatly increase the efficiency of a run. If only partial functionality is needed (e.g. Only analyzing modified ribonucleotides), users can implement flags to activate  function modules desired. See below for more details. 

## Command Line Arguments and Description

Read the [doc](https://chosenobih.github.io/hamrlinc_docs/Tutorial/) for detailed descriptions on selected flags.

| Command | Description |
| :---: | :---: |
| Required |
| -o | \<pipeline output directory\> <br> name of the directory where you would like your hamrlinc run to be |
| -c | \<filenames for each fastq.csv\> <br> a csv file that corresponds each srr code (or name of fastq file) to your desired nomenclature for each read |
| -g | \<reference genome.fa> <br> a fasta file of the genome of the model organism |
| -i | \<reference genome annotation.gff3> <br> a gff3 file of the genome of the model organism, note we require gff3 instead of gtf |
| -l | \<read length\> <br> an integer, the read length of this sequencing experiment, if non-unanimous use the shortest length |
| -s | \<genome size in bp\> <br> an integer, the number of base pairs of the genome of this model organism |
| Optional |
| -n | \[number of threads\] <br> default=4 |
| -d | \[raw fastq folder\] <br> a path to a folder containing raw fastq files if needed
| -a | \[use Tophat2 instead of STAR\] <br> default uses STAR |
| -b | \[Tophat2 library choice: fr-unstranded, fr-firststrand, fr-secondstrand\] <br> default=fr-firststrand |
| -x | \[Genome index directory for tophat2 by user input\] <br> default=None|
| -f | \[filter\] <br> default=filter_SAM_number_hits.pl |
| -k | \[activate modification analysis (left arm)\] |
| -p | \[activate lincRNA identification (inner right arm)\] |
| -u | \[activate regular featurecount (outer right arm)\] |
| -v | \[evolinc option: M or MO\] <br> default=M |
| -Q | \[HAMR minimum quality score: 0-40\] <br> default=30 |
| -C | \[HAMR minimum coverage: 0-∞\] <br> default=10 |
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
| -S | \[path to hamr.py for internal testing and server usage\] <br> default=/HAMR/hamr.py |
| -H | \[path to pthr_go_annots.py for internal testing and server usage\] <br> default=/HAMR/hamr.py |

Running HAMRLINC
-----------------------

### Required dependencies
1. Linux-based computer, server, or cluster.
2. [Docker](https://docs.docker.com/engine/install/)
3. Minimum memory of 32 GB and minimum disk space of 120 GB. 

```
#pull HAMRLINC docker image:  
docker pull chosenobih/hamrlinc:v0.3
```
```
#clone HAMRLINC repo
git clone https://github.com/chosenobih/HAMRLINC.git
cd HAMRLINC
```
```
#download the genome file for Arabidopsis thaliana from ENSEMBL
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```
```
#download the annotation file for Arabidopsis thaliana from ENSEMBL
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.57.gff3.gz
gunzip Arabidopsis_thaliana.TAIR10.57.gff3.gz
```
```
#run HAMRLINC in SE mode with SRA IDs
docker run --rm -v $(pwd):/working-dir -w /working-dir chosenobih/hamrlinc:v0.3 -o hamrlinc_test -c /demo/PRJNA478205.csv -g Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -i Arabidopsis_thaliana.TAIR10.57.gff3 -l 50 -s 135000000 -n 8 -k
```

Running HAMRLINC as an application on CyVerse's Discovery Environment
---------------------------------------------------------------------
HAMRLINC has been integrated as an app on [CyVerse's Discovery Environment (DE)](https://de.cyverse.org/), and it is available for use by researchers. Search for “HAMRLINC" and then select the 1.0.0 version. A short tutorial on how to run the app is available at this [CyVerse wiki](https://cyverse.atlassian.net/wiki/spaces/DEapps/pages/1819639809/HAMRLINC+v1.0). CyVerse's DE provides an easy-to-use graphic user interphase for running several Life Sciences computational pipelines.

Step-by-step walkthrough
------------------------
For more detailed documentation and step-by-step tutorial for running HAMRLINC using docker, please visit [HAMRLINC Documentation page](https://chosenobih.github.io/hamrlinc_docs/)

Issues
------
If you encounter any issues while running HAMRLINC, please open an issue on this GitHub repo, and we'll attend to it as soon as possible.

Copyright
---------
```
Copyright (c) 2023 HAMRLINC Team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE 
```
