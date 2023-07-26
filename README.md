# MitoR

## Description
MitoR is a R package designed to detect single nucleotide polymorphisms (SNPs), insertions and deletions (indels), and copy number variations (CNVs) in the mitochondrial genome. 
This software provides a user-friendly tool that simplifies the execution of alignment, mapping, sorting, indexing, quality controls, detection, and visualization of SNPs, indels, and CNVs. Additionally, MitoR generates comprehensive reports that not only indicate the detected variations but also associates them with widely recognized databases such as HmtVar, VarSome, Franklin, and dbSNP. Moreover, it allows users to save all the obtained results, facilitating the creation of a local database for population analysis.
Users do not require advanced computer knowledge to utilize MitoR effectively.

## Instalation
To install this package execute the following command: 

`remotes::install_github("DaniOrschanski/MitoR")`

## Requirements
1. You must work on a Linux environment.
2. Download the Linux packages: "lsb_release", "bzip2", "libncurses5-dev"(for Ubuntu/Kali) or "ncurses-devel" (for Redhat/Fedora) , "unzip", "rpm2cpio", "cpio", "gzip", "tar", "wget", "java".
3. Organize each FASTA/FASTQ file of the patient in one unique folder which name must be the id of the correspondant patient.
