#' @title Creates a VCF file from reads R1 & R2 (fasta or fastq files)
#' @description Selecting the directory where the fasta or fastq files are saved, generate_VCF creates a VCF (Variant Call Format) file.
#' It uses a basic pipeline to create it, using the softwares: TrimGalore, BWA, PICARD, GATK, SamTools.
#' The filters used for keeping or removing mutation data are:
#' For SNPs: 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0'
#' For indels: 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'
#' @param keep_BAM=FALSE While creating the VCF file, a BAM file is created. If you want to keep it as a part of the analysis, it will be saved in the same directory as the VCF file.
#' @return Returns a the path of the created VCF file named "<patient>_filtered_MitoR.vcf".
#' Several files are created and saved in a new directory called "MitoR" at the same directory as the R1 and R2 reads were located:
#' 1 The VCF file
#' 2 The VCF index file. It is saved
#' 3 The BAM file (Optional. Default = FALSE)
#' In case you want to keep the BAM file, it will be placed at a new directory named of "BAMFile" at the same location as the VCF file.
#' @import magrittr
#' @export
#' @examples
#' R1_R2_directory_path = "~/Github/MitochondriaAnalysis/MitoR/Patients/37019"
#' generate_VCF(R1_R2_directory_path);
#' generate_VCF(R1_R2_directory_path, keep_BAM = TRUE)

generate_VCF <- function(path, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {

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

  mitor_sof <- sprintf("%s/mitorDB/Softwares", Sys.getenv('R_LIBS_USER'))
  SAMTOOLS <- sprintf("%s/Samtools/samtools-1.16.1/samtools", mitor_sof)
  BWA <- sprintf("%s/BWA/usr/bin/bwa", mitor_sof)
  PICARD <- sprintf("%s/picard-2.27.5/picard.jar", mitor_sof)
  GATK <- sprintf("%s/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_sof)
  reference <- sprintf("%s/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", mitor_sof)

  # Search for the R1 and R2
  patients <- searchPatient(path)
  patientR1 <- sprintf("%s/%s", path, patients[[1]])
  patientR2 <- sprintf("%s/%s", path, patients[[2]])
  id <- patients[[3]]

  # Generates the folder MitoR inside path directory
  setwd(path)
  dir.create(sprintf('%s/MitoR', path))
  setwd(sprintf('%s/MitoR', path))

  # Obtain number of threads of the PC
  cant_thread  <- (as.integer(system2("nproc", "--all", stdout = TRUE))-1)

  # Header-Contains information about the entire file, such as sample name, sample length, and alignment method.
  # Alignments in the alignments section are associated with specific information in the header section.
  header <- sprintf('@RG\\tID:%s\\tLB:MitoR\\tSM:bar', id)

  # BWA mem - Maps the alignment of R1 and R2 with the reference:
  system2(BWA, sprintf("mem -t %s -o %s/MitoR/%s_R1R2.sam %s %s %s -R '%s'", cant_thread, path, id, reference, patientR1, patientR2, header), stdout = TRUE, wait = TRUE)

  # STEP 3 SAM to BAM
  system2(SAMTOOLS, sprintf("view -bS %s/MitoR/%s_R1R2.sam > %s/MitoR/%s_R1R2.bam", path, id, path, id), stdout = TRUE, wait = TRUE)

  # STEP 4 SortSam
  system2("java", sprintf("-jar %s SortSam -I %s/MitoR/%s_R1R2.bam -O %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate", PICARD, path, id, path, id), stdout = TRUE, wait = TRUE)

  # STEP 5 MarkDuplicates
  system2("java", sprintf("-jar %s MarkDuplicates -I %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --CREATE_INDEX True --ASSUME_SORTED True -M %s/MitoR/%s_marked_dup_metrics.txt -O %s/MitoR/BAMfile/%s_dedup_reads.bam", PICARD, path, id, path, id, path, id), stdout = TRUE, wait = TRUE)

  # STEP 6 HaplotypeCaller
  system2("java", sprintf("-jar %s HaplotypeCaller -I %s/MitoR/%s_dedup_reads.bam -O %s/MitoR/%s_raw_variants.vcf -ip 100 -R %s", GATK, path, id, path, id, reference), stdout = TRUE, wait = TRUE)

  # STEP 7 SelectVariants
  system2("java", sprintf("-jar %s SelectVariants -V %s/MitoR/%s_raw_variants.vcf -O %s/MitoR/%s_raw_snps.vcf -R %s -select-type SNP", GATK, path, id, path, id, reference), stdout = TRUE, wait = TRUE)

  # STEP 8 SelectVariants
  system2("java", sprintf("-jar %s SelectVariants -V %s/MitoR/%s_raw_variants.vcf -O %s/MitoR/%s_raw_indels.vcf -R %s -select-type INDEL", GATK, path, id, path, id, reference), stdout = TRUE, wait = TRUE)

  # STEP 9
  system2("java", sprintf("-jar %s VariantFiltration -V %s/MitoR/%s_raw_snps.vcf -O %s/MitoR/%s_filtered_snps.vcf  -R %s --filter-expression 'QD %s || FS %s || MQ %s || MQRankSum %s || ReadPosRankSum %s' --filter-name 'mitor_indel_filter'", GATK, path, id, path, id, reference, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS), stdout = TRUE, wait = TRUE)

  # STEP 10
  system2("java", sprintf("-jar %s VariantFiltration -V %s/MitoR/%s_raw_indels.vcf -O %s/MitoR/%s_filtered_indels.vcf  -R %s --filter-expression 'QD %s || FS %s || ReadPosRankSum %s' --filter-name 'mitor_snp_filter'", GATK, path, id, path, id, reference, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS), stdout = TRUE, wait = TRUE)

  # STEP 11
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s/MitoR/%s_filtered_indels.vcf -O %s/MitoR/%s_filtered_indels_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path, id, path, id), stdout = TRUE, wait = TRUE)
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s/MitoR/%s_filtered_snps.vcf -O %s/MitoR/%s_filtered_snps_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path, id, path, id), stdout = TRUE, wait = TRUE)

  # STEP 12
  system2("java", sprintf("-jar %s MergeVcfs -I %s/MitoR/%s_filtered_snps_tomerge.vcf -I %s/MitoR/%s_filtered_indels_tomerge.vcf -O %s/MitoR/%s_filtered_MitoR.vcf", PICARD, path, id, path, id, path, id), stdout = TRUE, wait = TRUE)

  # Remove files:
  file.remove(sprintf("%s/MitoR/%s_R1R2.sam", path, id))
  file.remove(sprintf("%s/MitoR/%s_R1R2.bam", path, id))
  #file.remove(sprintf("%s/MitoR/%s_dedup_reads.bam", path, id))
  file.remove(sprintf("%s/MitoR/%s_raw_variants.vcf", path, id))
  file.remove(sprintf("%s/MitoR/%s_raw_snps.vcf", path, id))
  file.remove(sprintf("%s/MitoR/%s_raw_indels.vcf", path, id))
  file.remove(sprintf("%s/MitoR/%s_filtered_indels.vcf", path, id))
  file.remove(sprintf("%s/MitoR/%s_filtered_snps.vcf", path, id))
  file.remove(sprintf("%s/MitoR/%s_filtered_indels_tomerge.vcf", path, id))
  file.remove(sprintf("%s/MitoR/%s_filtered_snps_tomerge.vcf", path, id))
  file.remove(sprintf("%s/MitoR/%s_raw_variants.vcf.idx", path, id))
  file.remove(sprintf("%s/MitoR/%s_raw_snps.vcf.idx", path, id))
  file.remove(sprintf("%s/MitoR/%s_raw_indels.vcf.idx", path, id))
  file.remove(sprintf("%s/MitoR/%s_filtered_indels.vcf.idx", path, id))
  file.remove(sprintf("%s/MitoR/%s_filtered_snps.vcf.idx", path, id))
  file.remove(sprintf("%s/MitoR/%s_filtered_snps_tomerge.vcf.idx", path, id))

  return(sprintf("%s/MitoR/%s_filtered_MitoR.vcf", path, id))
}

#' @title searchPatient
#' @description search for R1 and R2 fasta/fastq files from the path of the patient's folder.
#' @param path of the directory of the patient. It must contain fasta/q files.
#' @return List:
#' - Path of fasta/q file - R1
#' - Path of fasta/q  file - R2
#' - Id of the patient.
#' @examples
#' searchPatient(~/MitoR/patients/39337)

searchPatient <- function(path) {
  files <- unlist(list.files(path))

  patientR1 <- files[stringr::str_detect(files, "R1")]
  if (length(patientR1) > 1) {
    accepted_files <- c(".fasta", ".fastq")
    patientR1 <- patientR1[stringr::str_detect(patientR1, accepted_files)]
    if (length(patientR1) > 1) {
      patientR1 <- patientR1[stringr::str_detect(patientR1, ".fastq")]
      if (length(patientR1) > 1) {
        stop("In the path you provided there are several .fasta or .fastq files containing the expression 'R1'. Please keep just one of them so that we take the correct one.")
      }
    } else if (is.null(patientR1)) {
      stop("Not able to find R1 fasta or fastq file. Please make sure that the file name contains:
           - The expression 'R1'
           - It is a .fasta or .fastq file")
    }
  }
  patientR2 <- files[stringr::str_detect(files, "R2")]
  if (length(patientR2) > 1) {
    accepted_files <- c(".fasta", ".fastq")
    patientR2 <- patientR2[stringr::str_detect(patientR2, accepted_files)]
    if (length(patientR2) > 1) {
      patientR2 <- patientR2[stringr::str_detect(patientR2, ".fastq")]
      if (length(patientR2) > 1) {
        stop ("In the path you provided there are several .fasta or .fastq files containing the expression 'R2'. Please keep just one of them so that we take the correct one.")
      }
    } else if (is.null(patientR2)) {
      stop ("Not able to find R2 fasta or fastq file. Please make sure that the file name contains:
           - The expression 'R2'
           - It is a .fasta or .fastq file")
    }
  }

  patient <- unlist(stringr::str_split(substr(stringr::str_extract(patientR1, pattern = "^(.+?)R1"), 1, nchar(stringr::str_extract(patientR1, pattern = "^(.+?)R1"))-2), ""))

  if (is.null(patient)) {
    patient <- basename(path)

  } else if (!(grepl('^[A-Za-z0-9]+$', patient[length(patient)], ignore.case = TRUE))) {
    patient <- paste(patient[-(length(patient))], sep = "", collapse = "")
  }
  return(list(patientR1, patientR2, patient))
}




