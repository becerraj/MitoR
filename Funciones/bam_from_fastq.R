#' @title Get BAM file from the FASTQ
#' @description Generation of BAM file and it's index from the FastQ format.
#' @param path_dir path od the directory which name is the id of the patient and contains the FastQ files (R1 and R2).
#' @param keepFastQ indicates if the FastQ files are removed (if FALSE) or not (TRUE) after the conversion to BAM
#' @return path of the BAM file
#' @examples
#' BAM_path <- bam_from_fastq(path_dir = "~/CNVAnalisisDefinitivo/PacientesMito/46608")
#' @export
#' @import Rsamtools


bam_from_fastq <- function(path_dir, keepFastQ = TRUE){

  titulo <- basename(path_dir)
  file_list <- list.files(path_dir)
  pacienteR1 <- file_list[grep("R1", file_list)]
  pacienteR2 <- file_list[grep("R2", file_list)]
  pacienteR1 <- sprintf("%s/%s", path_dir, pacienteR1)
  pacienteR2 <- sprintf("%s/%s", path_dir,  pacienteR2)

  # Generate a folder inside the patient's folder where the files will be saved
  dir.create(sprintf("%s/MitoR", path_dir))

  #dir_home <- path.expand("~")
  #mitor_files <- sprintf("%s/MitoR", Sys.getenv('R_LIBS_USER'))
  #referencia <- paste(mitor_files, "/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", sep="")

  # Threads of the PC
  cant_thread  <- (as.integer(system2("nproc", "--all", stdout = TRUE))-1)


  # Header-Contains information about the entire file, such as sample name,
  # sample length, and alignment method. Alignments in the alignments section
  # are associated with specific information in the header section.
  header <- sprintf('@RG\\tID:%s\\tLB:MitoR\\tSM:bar', titulo)

  # BWA mem - Alignment of R1 and R2 from the patient with the reference
  system2(BWA, sprintf("mem -t %s -o %s/MitoR/%s_R1R2.sam %s %s %s -R '%s'", cant_thread, path_dir, titulo, referencia, pacienteR1, pacienteR2, header), stdout = TRUE, wait = TRUE)
  #bwa mem -t 3 -o /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/MitoR/37339_R1R2.sam /home/daniela/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/37339_R1.fastq.gz /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/37339_R2.fastq.gz -R @RG\\tID:37339\\tLB:MitoR\\tSM:bar

  # SAM to BAM
  system2(SAMTOOLS, sprintf("view -bS %s/MitoR/%s_R1R2.sam > %s/MitoR/%s_R1R2.bam", path_dir, titulo, path_dir, titulo), stdout = TRUE, wait = TRUE)
  #/home/daniela/MitoRSoftware/Samtools/samtools-1.16.1/samtools view -bS /home/daniela/Paciente/44562/MitoR/44562_R1R2.sam > /home/daniela/Paciente/44562/MitoR/44562_R1R2.bam


  # SortSam (SortBAM)
  system2("java", sprintf("-jar %s SortSam -I %s/MitoR/%s_R1R2.bam -O %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate", PICARD, path_dir, titulo, path_dir, titulo), stdout = TRUE, wait = TRUE)

  # index sorted Bam (creates .bam.bai)
  indexBam(sprintf("%s/MitoR/%s_sortedR1R2.bam", path_dir, titulo))

  #Remove big files
  file.remove(sprintf("%s/MitoR/%s_R1R2.sam", path_dir, id))
  file.remove(sprintf("%s/MitoR/%s_R1R2.bam", path_dir, id))

  return(sprintf("%s/MitoR/%s_sortedR1R2.bam", path_dir, titulo))
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
  #bwa mem -t 3 -o /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/MitoR/37339_R1R2.sam /home/daniela/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/37339_R1.fastq.gz /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/37339_R2.fastq.gz -R @RG\\tID:37339\\tLB:MitoR\\tSM:bar

  # SAM to BAM
  system2(Samtools, sprintf("view -bS %s/MitoR/%s_R1R2.sam > %s/MitoR/%s_R1R2.bam", dirPrinc, titulo, dirPrinc, titulo), stdout = TRUE, wait = TRUE)
  #/home/daniela/MitoRSoftware/Samtools/samtools-1.16.1/samtools view -bS /home/daniela/Paciente/44562/MitoR/44562_R1R2.sam > /home/daniela/Paciente/44562/MitoR/44562_R1R2.bam

  # SortSam (SortBAM)
  system2("java", sprintf("-jar %s SortSam -I %s/MitoR/%s_R1R2.bam -O %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate", Picard, dirPrinc, titulo, dirPrinc, titulo), stdout = TRUE, wait = TRUE)

  #Remove big files
  file.remove(sprintf("%s/MitoR/%s_R1R2.sam", path_file, titulo))
  file.remove(sprintf("%s/MitoR/%s_R1R2.bam", path_file, titulo))

  return(sprintf("%s/MitoR/%s_sortedR1R2.bam", dirPrinc, titulo))
}
