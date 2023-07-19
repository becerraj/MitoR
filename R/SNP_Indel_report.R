# MitoR conjunto
# SNP_Indel_report

# Funcion solo para analisis de SNP Indels. El path tiene que necesariamente ser el directorio donde estan los FASTA
# El directorio tiene que estar con el nombre del ID del paciente. Los R1 y R2 tienen que tener el nombre de: IDpacienteR1
SNP_Indel_Analyze <- function(path, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {
  wbs <- SNP_Indel_report(path, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0")

  numPaciente <- strsplit(path, "/") %>% unlist()
  numPaciente <- numPaciente[length(numPaciente) - 3]

  report_path <- sprintf("%s%s_SNP_Indels_Report", substr(path, start = 0, stop = (nchar(path)-nchar(basename(path)))), numPaciente)

  wb_SNP_Indel_report <- wbs[[1]]
  wb_SNP_Indel_soft <- wbs[[2]]
  wb_SNP_Indel_plot <- wbs[[3]]

  workbooks <- list(wb_SNP_Indel_report, wb_SNP_Indel_soft, wb_SNP_Indel_plot)

  XLSX_file <- combinarWorkbooks(workbooks, report_path)

  openxlsx::saveWorkbook(XLSX_file, report_path,  overwrite = TRUE)
}

# Funcion para usar en el analisis particular o general. El input puede ser el BAM o puede ser el directorio.
SNP_Indel_report <- function(path, QD_SNPS = "< 2.0", FS_SNPS = "> 60.0", MQ_SNPS = "< 40.0", MQRankSum_SNPS = "< -12.5", ReadPosRankSum_SNPS = "< -8.0", QD_INDELS = "< 2.0", FS_INDELS = "> 200.0", ReadPosRankSum_INDELS = "< -20.0") {

  file_type <- strsplit(path, "\\.") %>% unlist()
  file_type <- file_type[length(file_type)]

  if (!file_type == "bam") {
    path <- fasta_to_bam(path)
  }

  VCF_file <- bam_to_vcf(path, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  wb2 <- generateXLSX(VCF_file)

  wb_SNP_Indel_report <- wb_SNP_Indel[1]
  wb_SNP_Indel_soft <- wb_SNP_Indel[2]

  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))
  numPaciente <- names(RDS_DB)[length(names(RDS_DB))]
  wb_SNP_Indel_plot <- pileupPlot(numPaciente)

  return(list(wb_SNP_Indel_report, wb_SNP_Indel_soft, wb_SNP_Indel_plot))
}


bam_to_vcf <- function(sorted_BAM) {
  path_to_MitoR <- substr(sorted_BAM, start=0, stop=(nchar(sorted_BAM)-nchar(basename(sorted_BAM))))
  titulo <- strsplit(path_to_MitoR, "/") %>% unlist()
  titulo <- titulo[length(titulo) - 1]
  # STEP 5 MarkDuplicates
  system2("java", sprintf("-jar %s MarkDuplicates -I %s/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --CREATE_INDEX True --ASSUME_SORTED True -M %s/%s_marked_dup_metrics.txt -O %s/BAMfile/%s_dedup_reads.bam", PICARD, sorted_BAM, path_to_MitoR, titulo, path_to_MitoR, titulo), stdout = TRUE, wait = TRUE)

  # STEP 6 HaplotypeCaller
  system2("java", sprintf("-jar %s HaplotypeCaller -I %s/%s_dedup_reads.bam -O %s/%s_raw_variants.vcf -ip 100 -R %s", GATK, path_to_MitoR, titulo, path_to_MitoR, titulo, referencia), stdout = TRUE, wait = TRUE)

  # STEP 7 SelectVariants
  system2("java", sprintf("-jar %s SelectVariants -V %s/%s_raw_variants.vcf -O %s/%s_raw_snps.vcf -R %s -select-type SNP", GATK, path_to_MitoR, titulo, path_to_MitoR, titulo, referencia), stdout = TRUE, wait = TRUE)

  # STEP 8 SelectVariants
  system2("java", sprintf("-jar %s SelectVariants -V %s/%s_raw_variants.vcf -O %s/%s_raw_indels.vcf -R %s -select-type INDEL", GATK, path_to_MitoR, titulo, path_to_MitoR, titulo, referencia), stdout = TRUE, wait = TRUE)

  # STEP 9
  system2("java", sprintf("-jar %s VariantFiltration -V %s/%s_raw_snps.vcf -O %s/%s_filtered_snps.vcf  -R %s --filter-expression 'QD %s || FS %s || MQ %s || MQRankSum %s || ReadPosRankSum %s' --filter-name 'mitor_indel_filter'", GATK, path_to_MitoR, titulo, path_to_MitoR, titulo, referencia, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS), stdout = TRUE, wait = TRUE)

  # STEP 10
  system2("java", sprintf("-jar %s VariantFiltration -V %s/%s_raw_indels.vcf -O %s/%s_filtered_indels.vcf  -R %s --filter-expression 'QD %s || FS %s || ReadPosRankSum %s' --filter-name 'mitor_snp_filter'", GATK, path_to_MitoR, titulo, path_to_MitoR, titulo, referencia, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS), stdout = TRUE, wait = TRUE)

  # STEP 11
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s/%s_filtered_indels.vcf -O %s/%s_filtered_indels_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path_to_MitoR, titulo, path_to_MitoR, titulo), stdout = TRUE, wait = TRUE)
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s/%s_filtered_snps.vcf -O %s/%s_filtered_snps_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path_to_MitoR, titulo, path_to_MitoR, titulo), stdout = TRUE, wait = TRUE)

  # STEP 12
  system2("java", sprintf("-jar %s MergeVcfs -I %s/%s_filtered_snps_tomerge.vcf -I %s/%s_filtered_indels_tomerge.vcf -O %s/%s_filtered_MitoR.vcf", PICARD, path_to_MitoR, titulo, path_to_MitoR, titulo, path_to_MitoR, titulo), stdout = TRUE, wait = TRUE)

  # Eliminamos los archivos generados que no son de utilidad
  file.remove(sorted_BAM)
  file.remove(sprintf("%s/%s_R1R2.sam", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_R1R2.bam", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_raw_variants.vcf", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_raw_snps.vcf", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_raw_indels.vcf", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_filtered_indels.vcf", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_filtered_snps.vcf", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_filtered_indels_tomerge.vcf", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_filtered_snps_tomerge.vcf", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_raw_variants.vcf.idx", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_raw_snps.vcf.idx", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_raw_indels.vcf.idx", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_filtered_indels.vcf.idx", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_filtered_snps.vcf.idx", path_to_MitoR, titulo))
  file.remove(sprintf("%s/%s_filtered_snps_tomerge.vcf.idx", path_to_MitoR, titulo))

  return(sprintf("%s/%s_filtered_MitoR.vcf", path_to_MitoR, titulo))
}


#' @title bam_from_fastq
#' @description Al darle un path de un directorio, va a buscar los archivos R1 y R2 que ademas sean .fasta o .fastq
#' @param path es el directorio donde tiene que buscar los dos archivos
#' @return Nos devuelve tres elementos:
#' - Archivo .fasta o .fastq R1
#' - Archivo .fasta o .fastq R2
#' - Nombre del paciente. Es decir, la identificacion que poseen los archivos R1 y R2 en comun
fasta_to_bam <- function(path_dir) {
  titulo <- basename(path_dir)
  file_list <- list.files(path_dir)
  pacienteR1 <- file_list[grep("R1", file_list)]
  pacienteR2 <- file_list[grep("R2", file_list)]
  pacienteR1 <- sprintf("%s/%s", path_dir, pacienteR1)
  pacienteR2 <- sprintf("%s/%s", path_dir,  pacienteR2)

  # Generate a folder inside the patient's folder where the files will be saved
  dir.create(sprintf("%s/MitoR", path_dir))

  mitor_files <- sprintf("%s/MitoR", Sys.getenv('R_LIBS_USER'))
  referencia <- paste(mitor_files, "/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", sep="")

  # Threads of the PC
  cant_thread  <- (as.integer(system2("nproc", "--all", stdout = TRUE))-1)

  # Header-Contains information about the entire file, such as sample name,
  # sample length, and alignment method. Alignments in the alignments section
  # are associated with specific information in the header section.
  header <- sprintf('@RG\\tID:%s\\tLB:MitoR\\tSM:bar', titulo)

  # BWA mem - Alignment of R1 and R2 from the patient with the reference
  system2(BWA, sprintf("mem -t %s -o %s/MitoR/%s_R1R2.sam %s %s %s -R '%s'", cant_thread, dirPrinc, titulo, referencia, pacienteR1, pacienteR2, header), stdout = TRUE, wait = TRUE)

  # SAM to BAM
  system2(Samtools, sprintf("view -bS %s/MitoR/%s_R1R2.sam > %s/MitoR/%s_R1R2.bam", dirPrinc, titulo, dirPrinc, titulo), stdout = TRUE, wait = TRUE)

  # SortSam (SortBAM)
  system2("java", sprintf("-jar %s SortSam -I %s/MitoR/%s_R1R2.bam -O %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate", Picard, dirPrinc, titulo, dirPrinc, titulo), stdout = TRUE, wait = TRUE)

  #Remove big files
  file.remove(sprintf("%s/MitoR/%s_R1R2.sam", path_file, titulo))
  file.remove(sprintf("%s/MitoR/%s_R1R2.bam", path_file, titulo))

  return(sprintf("%s/MitoR/%s_sortedR1R2.bam", dirPrinc, titulo))
}

# Función para combinar varios workbooks en uno
combinarWorkbooks <- function(workbooks, nombreArchivo) {
  # Crear un nuevo workbook de destino
  wb_SNP_Indels <- createWorkbook()

  # Iterar sobre los workbooks individuales
  for (wb_actual in workbooks) {
    # Leer las hojas del workbook actual
    hojas <- getSheetNames(wb_actual)

    # Iterar sobre las hojas del workbook actual
    for (hoja in hojas) {
      # Leer los datos de la hoja actual
      datos <- read.xlsx(wb_actual, sheet = hoja)

      # Agregar una nueva hoja al workbook de destino
      addWorksheet(wb_SNP_Indels, sheetName = hoja)

      # Escribir los datos en la hoja del workbook de destino
      writeData(wb_SNP_Indels, sheet = hoja, x = datos)
    }
  }

  return(wb_SNP_Indels)
}


