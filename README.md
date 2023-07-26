# MitoR

## Description
MitoR is an efficient R package designed to detect single nucleotide polymorphisms (SNPs), insertions and deletions (indels), and copy number variations (CNVs) in the mitochondrial genome. 

This software offers a user-friendly tool that simplifies the execution of alignment, mapping, sorting, indexing, quality controls, detection, and visualization of SNPs, indels, and CNVs. 

Additionally, MitoR generates comprehensive reports that not only indicate the detected variations but also associates them with widely recognized databases such as HmtVar, VarSome, Franklin, and dbSNP. Moreover, it allows users to save all the obtained results, facilitating the creation of a local database for population analysis.

The best part is that users do not need advanced computer knowledge to effectively use MitoR.

------------------------------------------------------------------------------------------------------
MitoR is a bioinformatic tool designed for detecting mitochondrial genome rearrangements. It offers a user-friendly analysis of single nucleotide polymorphisms (SNPs), insertions and deletions (indels), and copy number variations (CNVs) starting sequencing data. 

The software seamlessly integrates information from reputable databases like HMTVAR, VarSome, Franklin, and dbSNP, providing reliable and context-rich results. 

Additionally, MitoR presents visual representations for easy interpretation and generates a local database, enhancing the accuracy and accessibility of variant annotation. 

With its comprehensive features, MitoR empowers geneticists to gain valuable insights into the intricacies of the mitochondrial genome.
And the best part is that users do not need advanced computer knowledge to effectively use MitoR.

## Requirements
1. You must work on a Linux environment.
   
3. Download the following Linux packages: "lsb_release", "bzip2", "libncurses5-dev"(for Ubuntu/Kali) or "ncurses-devel" (for Redhat/Fedora) , "unzip", "rpm2cpio", "cpio", "gzip", "tar", "wget", "java".
   
5. Download the R package: devtools.
   
7. Organize each FASTA/FASTQ file of the patient in one unique folder which name must be the id of the correspondant patient.


## Installation
To install this package execute the following command: 

`remotes::install_github(".../MitoR")`


## Authors
- Engr. Juan Cruz Becerra
- Engr. Daniela Orschanski
- Dr. Engr. Elmer Andrés Fernández.
