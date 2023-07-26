# MitoR

## Description
MitoR is a bioinformatic tool designed for detecting mitochondrial genome rearrangements. It offers a user-friendly analysis of single nucleotide polymorphisms (SNPs), insertions and deletions (indels), and copy number variations (CNVs) starting sequencing data. 

MitoR generates comprehensive reports that not only indicate the detected variations but also associates them with widely recognized databases such as HmtVar, VarSome, Franklin, and dbSNP, providing reliable and context-rich results. These reports also include visual representations for easy interpretation, ensuring that users can readily explore the findings.

Moreover, it allows users to save all the obtained results, facilitating the creation of a local database for population analysis, enhancing the accuracy and accessibility of variant annotation.
 
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
