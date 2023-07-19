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

  parametros <- fix_filter_values(QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS,
                                  QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  QD_SNPS <- parametros[1]
  FS_SNPS <- parametros[2]
  MQ_SNPS <- parametros[3]
  MQRankSum_SNPS <- parametros[4]
  ReadPosRankSum_SNPS <- parametros[5]
  QD_INDELS <- parametros[6]
  FS_INDELS <- parametros[7]
  ReadPosRankSum_INDELS <- parametros[8]

  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  SAMTOOLS <- sprintf("%s/MitoRSoftware/Samtools/samtools-1.16.1/samtools", mitor_files)
  BWA <- sprintf("%s/MitoRSoftware/BWA/usr/bin/bwa", mitor_files)
  PICARD <- sprintf("%s/MitoRSoftware/picard-2.27.5/picard.jar", mitor_files)
  GATK <- sprintf("%s/MitoRSoftware/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", mitor_files)
  referencia <- sprintf("%s/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", mitor_files)

  # Pide la carpeta donde se encuentran las lecturas R1 y R2
  pacientes <- searchPatient(path)
  pacienteR1 <- sprintf("%s/%s", path, pacientes[[1]])
  pacienteR2 <- sprintf("%s/%s", path, pacientes[[2]])
  titulo <- pacientes[[3]]

  # Armamos las carpetas dentro del paciente: MitoR
  setwd(path)
  dir.create(sprintf('%s/MitoR', path))
  #system2("mkdir", sprintf('%s/MitoR', path), wait = TRUE)
  setwd(sprintf('%s/MitoR', path))

  # Buscar la cantidad de threads posibles para la maquina
  cant_thread  <- (as.integer(system2("nproc", "--all", stdout = TRUE))-1)

  # Header-Contains information about the entire file, such as sample name, sample length, and alignment method.
  # Alignments in the alignments section are associated with specific information in the header section.
  header <- sprintf('@RG\\tID:%s\\tLB:MitoR\\tSM:bar', titulo)

  # BWA mem - Mapea la alineacion del R1 y R2 con la referencia
  system2(BWA, sprintf("mem -t %s -o %s/MitoR/%s_R1R2.sam %s %s %s -R '%s'", cant_thread, path, titulo, referencia, pacienteR1, pacienteR2, header), stdout = TRUE, wait = TRUE)

  # STEP 3 SAM to BAM
  system2(SAMTOOLS, sprintf("view -bS %s/MitoR/%s_R1R2.sam > %s/MitoR/%s_R1R2.bam", path, titulo, path, titulo), stdout = TRUE, wait = TRUE)

  # STEP 4 SortSam
  system2("java", sprintf("-jar %s SortSam -I %s/MitoR/%s_R1R2.bam -O %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate", PICARD, path, titulo, path, titulo), stdout = TRUE, wait = TRUE)

  # STEP 5 MarkDuplicates
  system2("java", sprintf("-jar %s MarkDuplicates -I %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --CREATE_INDEX True --ASSUME_SORTED True -M %s/MitoR/%s_marked_dup_metrics.txt -O %s/MitoR/BAMfile/%s_dedup_reads.bam", PICARD, path, titulo, path, titulo, path, titulo), stdout = TRUE, wait = TRUE)

  # STEP 6 HaplotypeCaller
  system2("java", sprintf("-jar %s HaplotypeCaller -I %s/MitoR/%s_dedup_reads.bam -O %s/MitoR/%s_raw_variants.vcf -ip 100 -R %s", GATK, path, titulo, path, titulo, referencia), stdout = TRUE, wait = TRUE)

  # STEP 7 SelectVariants
  system2("java", sprintf("-jar %s SelectVariants -V %s/MitoR/%s_raw_variants.vcf -O %s/MitoR/%s_raw_snps.vcf -R %s -select-type SNP", GATK, path, titulo, path, titulo, referencia), stdout = TRUE, wait = TRUE)

  # STEP 8 SelectVariants
  system2("java", sprintf("-jar %s SelectVariants -V %s/MitoR/%s_raw_variants.vcf -O %s/MitoR/%s_raw_indels.vcf -R %s -select-type INDEL", GATK, path, titulo, path, titulo, referencia), stdout = TRUE, wait = TRUE)

  # STEP 9
  system2("java", sprintf("-jar %s VariantFiltration -V %s/MitoR/%s_raw_snps.vcf -O %s/MitoR/%s_filtered_snps.vcf  -R %s --filter-expression 'QD %s || FS %s || MQ %s || MQRankSum %s || ReadPosRankSum %s' --filter-name 'mitor_indel_filter'", GATK, path, titulo, path, titulo, referencia, QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS), stdout = TRUE, wait = TRUE)

  # STEP 10
  system2("java", sprintf("-jar %s VariantFiltration -V %s/MitoR/%s_raw_indels.vcf -O %s/MitoR/%s_filtered_indels.vcf  -R %s --filter-expression 'QD %s || FS %s || ReadPosRankSum %s' --filter-name 'mitor_snp_filter'", GATK, path, titulo, path, titulo, referencia, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS), stdout = TRUE, wait = TRUE)

  # STEP 11
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s/MitoR/%s_filtered_indels.vcf -O %s/MitoR/%s_filtered_indels_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path, titulo, path, titulo), stdout = TRUE, wait = TRUE)
  system2("java", sprintf("-jar %s RenameSampleInVcf -I %s/MitoR/%s_filtered_snps.vcf -O %s/MitoR/%s_filtered_snps_tomerge.vcf --NEW_SAMPLE_NAME bar", PICARD, path, titulo, path, titulo), stdout = TRUE, wait = TRUE)

  # STEP 12
  system2("java", sprintf("-jar %s MergeVcfs -I %s/MitoR/%s_filtered_snps_tomerge.vcf -I %s/MitoR/%s_filtered_indels_tomerge.vcf -O %s/MitoR/%s_filtered_MitoR.vcf", PICARD, path, titulo, path, titulo, path, titulo), stdout = TRUE, wait = TRUE)

  # Eliminamos los archivos generados que no son de utilidad
  file.remove(sprintf("%s/MitoR/%s_R1R2.sam", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_R1R2.bam", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_sortedR1R2.bam", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_raw_variants.vcf", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_raw_snps.vcf", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_raw_indels.vcf", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_filtered_indels.vcf", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_filtered_snps.vcf", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_filtered_indels_tomerge.vcf", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_filtered_snps_tomerge.vcf", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_raw_variants.vcf.idx", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_raw_snps.vcf.idx", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_raw_indels.vcf.idx", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_filtered_indels.vcf.idx", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_filtered_snps.vcf.idx", path, titulo))
  file.remove(sprintf("%s/MitoR/%s_filtered_snps_tomerge.vcf.idx", path, titulo))

  return(sprintf("%s/MitoR/%s_filtered_MitoR.vcf", path, titulo))
}

#' @title searchPatient
#' @description Al darle un path de un directorio, va a buscar los archivos R1 y R2 que ademas sean .fasta o .fastq
#' @param path es el directorio donde tiene que buscar los dos archivos
#' @return Nos devuelve tres elementos:
#' - Archivo .fasta o .fastq R1
#' - Archivo .fasta o .fastq R2
#' - Nombre del paciente. Es decir, la identificacion que poseen los archivos R1 y R2 en comun
#' @examples
#' searchPatient(~/MitoR/Pacientes/39337)
#' @import stringr

searchPatient <- function(path) {
  files <- unlist(list.files(path))
  # Del path provisto, tomamos los que tengan en su nombre R1 y R2, y ademas sean archivo fasta o fastq
  pacienteR1 <- files[stringr::str_detect(files, "R1")]
  if (length(pacienteR1) > 1) {
    accepted_files <- c(".fasta", ".fastq")
    pacienteR1 <- pacienteR1[stringr::str_detect(pacienteR1, accepted_files)]
    if (length(pacienteR1) > 1) {
      pacienteR1 <- pacienteR1[stringr::str_detect(pacienteR1, ".fastq")]
      if (length(pacienteR1) > 1) {
        stop("In the path you provided there are several .fasta or .fastq files containing the expression 'R1'. Please keep just one of them so that we take the correct one.")
      }
    } else if (is.null(pacienteR1)) {
      stop("Not able to find R1 fasta or fastq file. Please make sure that the file name contains:
           - The expression 'R1'
           - It is a .fasta or .fastq file")
    }
  }
  pacienteR2 <- files[stringr::str_detect(files, "R2")]
  if (length(pacienteR2) > 1) {
    accepted_files <- c(".fasta", ".fastq")
    pacienteR2 <- pacienteR2[stringr::str_detect(pacienteR2, accepted_files)]
    if (length(pacienteR2) > 1) {
      pacienteR2 <- pacienteR2[stringr::str_detect(pacienteR2, ".fastq")]
      if (length(pacienteR2) > 1) {
        stop ("In the path you provided there are several .fasta or .fastq files containing the expression 'R2'. Please keep just one of them so that we take the correct one.")
      }
    } else if (is.null(pacienteR2)) {
      stop ("Not able to find R2 fasta or fastq file. Please make sure that the file name contains:
           - The expression 'R2'
           - It is a .fasta or .fastq file")
    }
  }
  # El titulo del paciente sera aquello que no fuera R1 o R2 del nombre anterior
  paciente <- unlist(stringr::str_split(substr(stringr::str_extract(pacienteR1, pattern = "^(.+?)R1"), 1, nchar(stringr::str_extract(pacienteR1, pattern = "^(.+?)R1"))-2), ""))
  # Llegado el caso que el archivo se llame solamente R1 o R2
  if (is.null(paciente)) {
    paciente <- basename(path)
    # Habiendo quitado el R1, si hay un caracter especial como ultimo valor, tambien lo sacamos
  } else if (!(grepl('^[A-Za-z0-9]+$', paciente[length(paciente)], ignore.case = TRUE))) {
    paciente <- paste(paciente[-(length(paciente))], sep = "", collapse = "")
  }
  return(list(pacienteR1, pacienteR2, paciente))
}

#' @title Fix the filter values
#' @description blabla
#' @param muchos blabla
#' @return parametros

fix_filter_values <- function(QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS,
                              QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS) {
  parametros <- c(QD_SNPS, FS_SNPS, MQ_SNPS, MQRankSum_SNPS, ReadPosRankSum_SNPS, QD_INDELS, FS_INDELS, ReadPosRankSum_INDELS)
  for (i in seq_along(parametros)) {
    parametros[i] <- stringr::str_replace(parametros[i], ",", ".")
    if (!grepl("\\.", parametros[i])) {
      parametros[i] <- paste0(parametros[i], ".0")
    }
  }
  return(parametros)
}
