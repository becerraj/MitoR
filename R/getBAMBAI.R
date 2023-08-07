#' @title Get the BAM file from the patient selected by the path.
#' @description Check if there are BAM or BAI files already generated on the directory. If there are not, generates them.
#' @param path_dir path of the directory of the patient that will be analyzed. It must contain R1 and R2 (fasta or fastq files).
#' @return bam.control path of the sorted BAM file

getBAMBAI <- function(path_dir) {
  if (dir.exists(paste(path_dir, "/MitoR", sep = ""))) {
    file_list <- list.files(paste(path_dir, "/MitoR", sep=""))
    bam_file <- file_list[endsWith(file_list, "sortedR1R2.bam")]
    if (!(length(nchar(bam_file)) == 0)) {
      bam.control <- sprintf("%s/MitoR/%s", path_dir, bam_file)
    } else {
      bam.control <- fasta_to_bam(path_dir)
    }
    bai_file <- file_list[endsWith(file_list, "sortedR1R2.bam.bai")]
    if (length(nchar(bai_file)) == 0) { indexBam(bam.control) }

  } else {
    bam.control <- fasta_to_bam(path_dir)
    indexBam(bam.control)
  }

  return(bam.control)
}
