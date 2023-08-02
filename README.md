# MitoR: 
## The complete Whole Genome Mitochondrial Analysis Toolkit in R.

## Description
MitoR is a bioinformatic tool designed for detecting mitochondrial genome rearrangements. It offers a user-friendly analysis of single nucleotide polymorphisms (SNPs), insertions and deletions (indels), and copy number variations (CNVs) from sequenced mithocondrial DNA data. 

MitoR generates complete reports that not only indicate the detected variations but also associates them with widely recognized databases such as HmtVar, VarSome, Franklin, and dbSNP, providing reliable and context-rich results. These reports also include visual representations for easy interpretation, ensuring that users can readily explore the findings.

Moreover, it allows users to save all the obtained results, facilitating the creation of a local database for population analysis, enhancing the accuracy and accessibility of variant annotation from local population, thus allowing to filter out highly common variants in the analyzed population.
 
With its comprehensive features, MitoR empowers geneticists to gain valuable insights into the intricacies of the mitochondrial genome.
And the best part is that users do not need advanced computer knowledge to effectively use MitoR.

## Requirements
1. You must work on a Linux environment.
   
3. Download the following Linux packages: "lsb_release", "bzip2", "libncurses5-dev"(for Ubuntu/Kali) or "ncurses-devel" (for Redhat/Fedora) , "unzip", "rpm2cpio", "cpio", "gzip", "tar", "wget", "java".
   
4. Organize each FASTA/FASTQ file of the patient in one unique folder which name must be the id of the correspondant patient.


## Installation
To install this package execute the following command: 

```
install.packages("devtools")
library(devtools)
install_github("DaniOrschanski/MitoR")
```
## Authors
* **Elmer Andrés Fernández** - *Initial work & Idea* - [Profile](https://www.researchgate.net/profile/Elmer_Fernandez) - [CIDIE]- [CONICET](http://www.conicet.gov.ar) - [Fundación para el Progreso de la Medicina - FPM](https://fpmlab.org.ar/)
* **Ing. Biom. Juan Becerra** -*Developer & Maintener
* **Ing. Biom. Daniela Orchansky** -*Developer & Maintener
## Collaborators
* **Horacio Martinetto** - *Genetic Advisor application [Profile](https://www.fleni.org.ar/profesionales/martinetto-horacio/), [FLENI](https://www.fleni.org.ar/)
* **Agata Fernandez** - *Genetic Advisor application [FLENI](https://www.fleni.org.ar/)
* **Juan Carloz Vazquez** - *Developer advisor

