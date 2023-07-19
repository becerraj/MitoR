#' @title Generation of BAM file from FASTA or FASTQ files.
#' @description Given the path of a directory, this function will identify the R1 and R2 files
#' and will generate a SAM file, then a BAM and a sorted BAM file by applying commands from BWA,
#' Picard, Samtools and GATK and using the mitochondrial human genome reference.
#' All this softwares and reference are downloaded automathically with the installation of MitoR.
#' @param path_dir  path of the directory of the patient that will be analyzed.
#' @return path of the sorted BAM file.

fasta_to_bam <- function(path_dir) {
  mitor_files <- path.expand("~")
  #mitor_files <- sprintf("%s/MitoR", Sys.getenv('R_LIBS_USER'))

  SAMTOOLS <- sprintf("%s/MitoRSoftware/Samtools/samtools-1.16.1/samtools", mitor_files)
  BWA <- sprintf("%s/MitoRSoftware/BWA/usr/bin/bwa", mitor_files)
  PICARD <- sprintf("%s/MitoRSoftware/picard-2.27.5/picard.jar", mitor_files)
  GATK <- sprintf("%s/MitoRSoftware/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_files)
  reference <- sprintf("%s/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", mitor_files)

  id <- basename(path_dir)
  file_list <- list.files(path_dir)
  patientR1 <- file_list[grep("R1", file_list)]
  patientR2 <- file_list[grep("R2", file_list)]
  patientR1 <- sprintf("%s/%s", path_dir, patientR1)
  patientR2 <- sprintf("%s/%s", path_dir,  patientR2)

  # Generate a folder inside the patient's folder where the files will be saved
  dir.create(sprintf("%s/MitoR", path_dir))

  # Threads of the PC
  cant_thread  <- (as.integer(system2("nproc", "--all", stdout = TRUE))-1)

  # Header-Contains information about the entire file, such as sample name,
  # sample length, and alignment method. Alignments in the alignments section
  # are associated with specific information in the header section.
  header <- sprintf('@RG\\tID:%s\\tLB:MitoR\\tSM:bar', id)

  # BWA mem - Alignment of R1 and R2 from the patient with the reference
  system2(BWA, sprintf("mem -t %s -o %s/MitoR/%s_R1R2.sam %s %s %s -R '%s'", cant_thread, path_dir, id, reference, patientR1, patientR2, header), stdout = TRUE, wait = TRUE)

  # SAM to BAM
  system2(SAMTOOLS, sprintf("view -bS %s/MitoR/%s_R1R2.sam > %s/MitoR/%s_R1R2.bam", path_dir, id, path_dir, id), stdout = TRUE, wait = TRUE)

  # SortSam (SortBAM)
  system2("java", sprintf("-jar %s SortSam -I %s/MitoR/%s_R1R2.bam -O %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate", PICARD, path_dir, id, path_dir, id), stdout = TRUE, wait = TRUE)

  #Remove big files
  file.remove(sprintf("%s/MitoR/%s_R1R2.sam", path_dir, id))
  file.remove(sprintf("%s/MitoR/%s_R1R2.bam", path_dir, id))

  return(sprintf("%s/MitoR/%s_sortedR1R2.bam", path_dir, id))
}
