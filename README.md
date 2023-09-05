# HAMRLINC: High-throughput Analysis of Modified Ribonucleotides and Long Intergenic Non-Coding RNAs
![HAMRbox_Workflow_v4](https://github.com/harrlol/HAMRLINC/assets/87460010/bb27dd18-dd9d-45e7-be6e-f295473121c5)


## Overview
- HAMRLINC is a multipurpose toolbox that expedites the analysis pipeline of two algorithms: [HAMR](https://github.com/GregoryLab/HAMR) and [Evolinc](https://github.com/Evolinc/Evolinc-I/tree/master). The former was developed by [Paul Ryvkin et al](https://rnajournal.cshlp.org/content/19/12/1684), and the latter by [Andrew D.L. Nelson et al](https://www.frontiersin.org/articles/10.3389/fgene.2017.00052/full). HAMRLINC aims to make the original methods more accessible by automating the tedious pre-processing steps, allowing users to analyze RNA-seq data at an experiment scale. 
- HAMRLINC is high-throughput and performs RNA-modification analysis and long intergenic non-coding RNAs(lincRNA) annotation at a bioproject scale. HAMRLINC performs constitutive trimming of acquired reads using Trim-Galore, and makes use of STAR as the aligning tool which reduces runtime, while allowing the users to use Tophat as an option. 


