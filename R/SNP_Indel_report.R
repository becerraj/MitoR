#' @title Analyze a DNA sequence to look for SNPs & INDELs
#' @description It generates a XLSX file with the SNPs & INDELs found, with information about them from HMTVAR. It also creates hyperlinks to Franklin, Varsome, dbSNP for each variant.
#' It will also generate a VCF file which will be saved in the same directory as the input path is set. The VCF is filtered with the setted parameters.
#' The analysis will be automatically loaded in MitoR DataBase.
#' @param path Either the BAM file or the directory containing the fasta or fastq files
#' @param ... Filters applied to the analysis of SNPs or INDELs to keep or remove a detected variant. Default numbers are the recommended by GATK
#' @return XLSX file containing several sheets are returned:
#' - SNP and INDEL report: for EACH variant you will get hyperlinks to Franklin, Varsome and dbSNP databases, information from HMTVAR database, filtering information from GATK and PICARD, allele depth and coverage information.
#' - Plot report: for EACH mutated gene you will get a single plot showing the read coverage. And for each mutation you will also get a single plot with the same information, showing the sequence with the centered variant and the next and previous nucleotides (range of 5 nucleotides)
#' - Softwares report: Information regarding the versions, dates of analysis and data from the used softwares (BWA, GATK, PICARD, SamTools)
#' @export
#' @examples
#' SNP_Indel_report("../MitoR/patient1.vcf")
#' generateXLSX("../MitoR", QD_SNPS = < 2.5, FS_SNPS = > 45.0)

SNP_Indel_Analyze <- function(path, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {
  wbs <- SNP_Indel_report(path, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0")

  numPaciente <- strsplit(path, "/") %>% unlist()
  numPaciente <- numPaciente[length(numPaciente) - 3]

  if (substr(path, nchar(path) - 3, nchar(path)) == ".bam") {
    report_path <- sprintf("%s%s_SNP_Indels_Report", substr(path, start = 0, stop = (nchar(path)-nchar(basename(path)))), numPaciente)
  } else {
    report_path <- sprintf("%s/MitoR/%s_SNP_Indels_Report", path, numPaciente)
  }

  wb_SNP_Indel_patient <- wbs[[1]]
  wb_SNP_Indel_soft <- wbs[[2]]
  wb_SNP_Indel_plot <- wbs[[3]]

  # Merge all the WB into one to create the XLSX file
  workbooks <- list(wb_SNP_Indel_patient, wb_SNP_Indel_plot, wb_SNP_Indel_soft)
  sheetNames <- list("SNP-INDEL", "Plot SNP-INDEL", "Softwares")
  XLSX_file <- combine_workbooks(workbooks, sheetNames)

  openxlsx::saveWorkbook(XLSX_file, report_path,  overwrite = TRUE)
}

#' @title Analyze a DNA sequence to look for SNPs & INDELs
#' @description It detects SNPs & INDELs from a BAM or fasta/fastq file, and adds information about them from HMTVAR. It also creates hyperlinks to Franklin, Varsome, dbSNP for each variant.
#' It will also generate a VCF file which will be saved in the same directory as the input path is set. The VCF is filtered with the setted parameters.
#' It creates plots for each mutated gene and each variant.
#' The analysis will be automatically loaded in MitoR DataBase.
#' @param path Either the BAM file or the directory containing the fasta or fastq files
#' @param ... Filters applied to the analysis of SNPs or INDELs to keep or remove a detected variant. Default numbers are the recommended by GATK
#' @return Several R workbooks (from openxlsx package) are returned:
#' SNP and INDEL report: for EACH variant you will get hyperlinks to Franklin, Varsome and dbSNP databases, information from HMTVAR database, filtering information from GATK and PICARD, allele depth and coverage information.
#' Softwares report: Information regarding the versions, dates of analysis and data from the used softwares (BWA, GATK, PICARD, SamTools)
#' Plot report: for EACH mutated gene you will get a single plot showing the read coverage. And for each mutation you will also get a single plot with the same information, showing the sequence with the centered variant and the next and previous nucleotides (range of 5 nucleotides)
#' @examples
#' SNP_Indel_report("../MitoR/patient1.vcf")
#' generateXLSX("../MitoR", QD_SNPS = < 2.5, FS_SNPS = > 45.0)

SNP_Indel_report <- function(path, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {

  #Check if there are BAM or BAI files already generated or generates them:
  path <- getBAMBAI(path)

  # Modifies the parameters to the exact way GATK needs to understand them
  params <- fix_filter_values(QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS,
                                  QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  QD_SNPS <- params[1]
  FS_SNPS <- params[2]
  MQ_SNPS <- params[3]
  MQRankSum_SNPS <- params[4]
  ReadPosRankSum_SNPS <- params[5]
  QD_INDELS <- params[6]
  FS_INDELS <- params[7]
  ReadPosRankSum_INDELS <- params[8]

  VCF_file <- bam_to_vcf(path, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  #VCF_file <- "~/Paciente/46608/MitoR/46608_filtered_MitoR.vcf"
  message("VCF file was successfully generated")

  ######################################################################################
  #wb_SNP_Indel <- xlsx_SNP_Indels_report
  wb_SNP_Indel <- generateXLSX(VCF_file)
  wb_SNP_Indel_patient <- wb_SNP_Indel[[1]]
  wb_SNP_Indel_soft <- wb_SNP_Indel[[2]]

  message("generate XLSX was successfully completed")

  # Look for the BAM file to do the pileupPlot function
  mitor_db <- sprintf("%s/mitorDB/DB", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/RDS_DB.rds", mitor_db))
  numPaciente <- names(RDS_DB)[length(names(RDS_DB))]

  BAM_file <- sprintf("%s%s_dedup_reads.bam", substr(VCF_file, start = 0, stop = (nchar(VCF_file)-nchar(basename(VCF_file)))), numPaciente)

  plots_to_xlsx <- pileupPlot(BAM_file, plot_mutation = TRUE, range = 5)
  message("The visualization for SNPs and indels was successfully generated")

  return(list(wb_SNP_Indel_patient, wb_SNP_Indel_soft, plots_to_xlsx))
}

#' @title Creates a VCF file from a BAM file
#' @description Selecting the path where the BAM file is saved, bam_to_vcf creates a VCF (Variant Call Format) file.
#' It uses a basic pipeline to create it, using the softwares: TrimGalore, BWA, PICARD, GATK, SamTools.
#' The default filters used for keeping or removing variants are:
#' For SNPs: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
#' For indels: 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'
#' @return Returns the path of the created VCF file named "<patient>_filtered_MitoR.vcf".
#' Several files are created and saved in a new directory called "MitoR" at the same directory as the R1 and R2 reads were located:
#' 1 The VCF file
#' 2 The VCF index file. It is saved
#' @import magrittr
#' @examples
#' 'bam_to_vcf("~/Github/MitochondriaAnalysis/MitoR/Patients/37019.bam")

bam_to_vcf <- function(sorted_BAM, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {

  path_to_MitoR <- substr(sorted_BAM, start=0, stop=(nchar(sorted_BAM)-nchar(basename(sorted_BAM))))
  patient_ID <- strsplit(path_to_MitoR, "/") %>% unlist()
  patient_ID <- patient_ID[length(patient_ID) - 1]

  mitor_sof <- sprintf("%s/mitorDB/Softwares", Sys.getenv('R_LIBS_USER'))
  SAMTOOLS <- sprintf("%s/Samtools/samtools-1.16.1/samtools", mitor_sof)
  BWA <- sprintf("%s/BWA/usr/bin/bwa", mitor_sof)
  PICARD <- sprintf("%s/picard-2.27.5/picard.jar", mitor_sof)
  GATK <- sprintf("%s/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_sof)
  reference <- sprintf("%s/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", mitor_sof)

  # STEP 5 MarkDuplicates
  system2("java", sprintf("-jar %s MarkDuplicates -I %s --VALIDATION_STRINGENCY SILENT --CREATE_INDEX True --ASSUME_SORTED True -M %s%s_marked_dup_metrics.txt -O %s%s_dedup_reads.bam", PICARD, sorted_BAM, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID), stdout = TRUE, wait = TRUE)

  # STEP 6 HaplotypeCaller
  system2("java", sprintf("-jar %s HaplotypeCaller -I %s%s_dedup_reads.bam -O %s%s_raw_variants.vcf -ip 100 -R %s", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference), stdout = TRUE, wait = TRUE)

  # STEP 7 SelectVariants
  system2("java", sprintf("-jar %s SelectVariants -V %s%s_raw_variants.vcf -O %s%s_raw_snps.vcf -R %s -select-type SNP", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference), stdout = TRUE, wait = TRUE)

  # STEP 8 SelectVariants
  system2("java", sprintf("-jar %s SelectVariants -V %s%s_raw_variants.vcf -O %s%s_raw_indels.vcf -R %s -select-type INDEL", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference), stdout = TRUE, wait = TRUE)

  # STEP 9
  system2("java", sprintf("-jar %s VariantFiltration -V %s%s_raw_snps.vcf -O %s%s_filtered_snps.vcf  -R %s --filter-expression 'QD %s || FS %s || MQ %s || MQRankSum %s || ReadPosRankSum %s' --filter-name 'mitor_indel_filter'", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS), stdout = TRUE, wait = TRUE)

  # STEP 10
  system2("java", sprintf("-jar %s VariantFiltration -V %s%s_raw_indels.vcf -O %s%s_filtered_indels.vcf  -R %s --filter-expression 'QD %s || FS %s || ReadPosRankSum %s' --filter-name 'mitor_snp_filter'", GATK, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, reference, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS), stdout = TRUE, wait = TRUE)

  # STEP 11
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s%s_filtered_indels.vcf -O %s%s_filtered_indels_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID), stdout = TRUE, wait = TRUE)
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s%s_filtered_snps.vcf -O %s%s_filtered_snps_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID), stdout = TRUE, wait = TRUE)

  # STEP 12
  system2("java", sprintf("-jar %s MergeVcfs -I %s%s_filtered_snps_tomerge.vcf -I %s%s_filtered_indels_tomerge.vcf -O %s%s_filtered_MitoR.vcf", PICARD, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID, path_to_MitoR, patient_ID), stdout = TRUE, wait = TRUE)

  # Deletes the unnecessary files
  file.remove(sprintf("%s%s_R1R2.sam", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_R1R2.bam", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_variants.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_snps.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_indels.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_indels.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_snps.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_indels_tomerge.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_snps_tomerge.vcf", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_variants.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_snps.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_raw_indels.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_indels.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_snps.vcf.idx", path_to_MitoR, patient_ID))
  file.remove(sprintf("%s%s_filtered_snps_tomerge.vcf.idx", path_to_MitoR, patient_ID))

  return(sprintf("%s%s_filtered_MitoR.vcf", path_to_MitoR, patient_ID))
}




#' @title Fix the filter parameters for GATK
#' @description Modifies the parameters to the exact way GATK needs to understand them
#' @return List of the same parameters it had as input

fix_filter_values <- function(QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS,
                              QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS) {
  params <- c(QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  for (i in seq_along(params)) {
    params[i] <- stringr::str_replace(params[i], ",", ".")
    if (!grepl("\\.", params[i])) {
      params[i] <- paste0(params[i], ".0")
    }
  }
  return(params)
}

