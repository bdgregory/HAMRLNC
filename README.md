# HAMRLNC: High-throughput Annotation of Modified Ribonucleotides and Long Non-Coding RNAs
![HAMRLINC_workflow](https://github.com/harrlol/HAMRLINC/assets/87460010/cb994f27-f48a-4210-b0e7-e3c270923eab)

## Overview
- HAMRLNC is a multipurpose toolbox that expedites the analysis pipeline for [HAMR](https://github.com/GregoryLab/HAMR) and [Evolinc](https://github.com/Evolinc/Evolinc-I/tree/master). The former was developed by [Paul Ryvkin et al](https://rnajournal.cshlp.org/content/19/12/1684), and the latter by [Andrew D.L. Nelson et al](https://www.frontiersin.org/articles/10.3389/fgene.2017.00052/full). HAMRLNC aims to make the original methods more accessible by automating the tedious pre-processing steps and expanding on their functionalities with its built-in post-processing steps, allowing users to visualize epitranscriptomic analysis with experimental condition contexts.
- HAMRLNC is high-throughput and performs RNA-modification annotation and long intergenic non-coding RNAs(lincRNA) annotation at a bioproject scale. HAMRLNC performs constitutive trimming of acquired reads using Trim-Galore, and makes use of STAR as the default aligning tool; mapped reads are pre-processed using selected methods from [GATK](https://gatk.broadinstitute.org/hc/en-us), [gffread](https://ccb.jhu.edu/software/stringtie/gff.shtml#gffread), [CPC2](https://cpc2.gao-lab.org), [infernal](http://eddylab.org/infernal/), [samtools](http://www.htslib.org/doc/samtools.html), etc. Users can also opt to quantify transcripts alongside these steps. 
- HAMRLNC is optimized for partial parallel processing and modularization. Specifying a larger thread count where hardware permits will greatly increase the speed of a single run. If only partial functionality is needed (e.g. only analyzing modified ribonucleotides), users can implement flags to activate the function modules desired. See below for more details. 

## Command Line Arguments and Description

Read the [doc](https://chosenobih.github.io/hamrlinc_docs/Tutorial/) for detailed descriptions on selected flags.

| Command | Description |
| :---: | :---: |
| Required |
| -o | \<pipeline output directory\> <br> name of the directory where you would like your hamrlnc run to be |
| -c | \<filenames for each fastq.csv\> <br> a csv file that corresponds each srr code (or name of fastq file) to your desired nomenclature for each read |
| -g | \<reference genome.fa> <br> a fasta file of the genome of the model organism |
| -i | \<reference genome annotation.gff3> <br> a gff3 file of the genome of the model organism, note we require gff3 instead of gtf |
| `-l` | \<read length\> <br> an integer, the read length of this sequencing experiment, if non-unanimous use the shortest length |
| Optional |
| -n | \[number of threads\] <br> default=4 |
| -d | \[raw fastq folder\] <br> default=NA |
| -t | \[trim raw fastq\] <br> default=false |
| -D | \[raw bam folder\] <br> default=NA |
| -b | \[sort raw bam\] <br> default=false |
| `-I` | \[STAR genome index folder\] <br> default=NA |
| -k | \[activate modification annotation workflow\] <br> default=false |
| -p | \[activate lncRNA annotation workflow\] <br> default=false |
| -u | \[activate featurecount workflow\] <br> default=false |
| -f | \[HAMR filter\] <br> default=filter_SAM_number_hits.pl |
| -m | \[HAMR model\] <br> default=euk_trna_mods.Rdata |
| -Q | \[HAMR minimum quality score: 0-40\] <br> default=30 |
| -C | \[HAMR minimum coverage: 0-∞\] <br> default=10 |
| -E | \[HAMR sequencing error: 0-1\] <br> default=0.01 |
| -P | \[HAMR maximum p-value: 0-1\] <br> default=1 |
| -F | \[HAMR maximum FDR: 0-1\] <br> default=0.05 |
| -O | \[Panther [organism taxon ID](http://pantherdb.org/services/oai/pantherdb/supportedgenomes)\] <br> default="3702" |
| -A | \[Panther [annotation dataset](http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets)\] <br> default="GO:0008150" |
| -Y | \[Panther test type: FISHER or BINOMIAL\] <br> default="FISHER" |
| -R | \[Panther correction type: FDR, BONFERRONI, or NONE\] <br> default="FDR" |
| -y | \[keep intermediate bam files\] <br> default=false |
| -q | \[halt program upon completion of checkpoint 2\] <br> default=false |
| -G | \[attribute used for featurecount\] <br> default="gene_id" |
| -x | \[max intron length for lncRNA-annotation-unique STAR mapping\] <br> default=NA |
| -H | \[SERVER alt path for panther\] |
| -U | \[SERVER alt path for HAMRLNC\] |
| -W | \[SERVER alt path for GATK\]
| -S | \[SERVER alt path for HAMR\]
| -J | \[SERVER alt path for CPC2\]
| -M | \[SERVER alt path for Rfam\]
| -h | \[help message\]|

Running HAMRLNC
-----------------------

### Required dependencies
1. Linux-based computer, server, or cluster.
2. [Docker](https://docs.docker.com/engine/install/)
3. Minimum memory of 32 GB and minimum disk space of 120 GB, could require higher specs for organisms with larger genomes like human.

```
# pull HAMRLNC docker image:  
docker pull harrlol/hamrlnc:v0.1.0a
```
```
# clone HAMRLNC repo
git clone https://github.com/harrlol/HAMRLNC.git
cd HAMRLNC
```
```
# download the genome file for Arabidopsis thaliana from ENSEMBL
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
```
```
# download the annotation file for Arabidopsis thaliana from ENSEMBL
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gff3.gz
gunzip Arabidopsis_thaliana.TAIR10.59.gff3.gz
```
```
# make sure your fa and gff3 files are in your working directory, and enter that directory
cd /your/working/directory
# run HAMRLNC with SRA IDs with all three arms activated
docker run \
  --rm -v $(pwd):/working-dir \
  -w /working-dir harrlol/hamrlnc:v0.1.0a \
  -o test_run \
  -c /demo/demo_filenames.csv \
  -g Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  -i Arabidopsis_thaliana.TAIR10.59.gff3 \
  -l 50 -n 8 -k -p -u
```

Running HAMRLNC as an application on CyVerse's Discovery Environment
---------------------------------------------------------------------
HAMRLNC has been integrated as an app on [CyVerse's Discovery Environment (DE)](https://de.cyverse.org/), and it is available for use by researchers. Search for “HAMRLINC" and then select the 1.0.0 version. A short tutorial on how to run the app is available at this [CyVerse wiki](https://cyverse.atlassian.net/wiki/spaces/DEapps/pages/1819639809/HAMRLINC+v1.0). CyVerse's DE provides an easy-to-use graphic user interphase for running several Life Sciences computational pipelines.

Step-by-step walkthrough
------------------------
For more detailed documentation and step-by-step tutorial for running HAMRLNC using docker, please visit [HAMRLINC Documentation page](https://chosenobih.github.io/hamrlinc_docs/)

Issues
------
If you encounter any issues while running HAMRLNC, please open an issue on this GitHub repo, and we'll attend to it as soon as possible.

Copyright
---------
```
Copyright (c) 2024 HAMRLNC Team

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
