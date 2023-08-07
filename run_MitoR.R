##
library(MitoR)
library(stringr)
library(magrittr)
library(Rsamtools)

select <- rstudioapi::selectDirectory(
  caption = "Select Directory",
  label = "Select",
  path = "/home/biomolecular/DATA/NGS/MITOCONDRIAL/Pacientes/"
)


if(length(select)>0){

  dir.create(file.path("/mnt/data/Mitocondrial/Pacientes",basename(select)))
  lfo <- list.files(select, full.names = T, pattern = ".fastq.gz")

  ok <- file.copy(from = lfo, to = file.path("/mnt/data/Mitocondrial/Pacientes",basename(select)))
  # print(lfo)
  # subjPath <- file.path("/mnt/data/RNAseq/Patients",basename(select))
  # # subjPath<-select
  # lf <- list.files(subjPath, full.names = T, pattern = "1.fastq")
  #
  print(ok)

  mitor_analysis(file.path("/mnt/data/Mitocondrial/Pacientes",basename(select)))

} else {
  cat("Not selected")
}
# FusionPlot(sorted.bam.file, savePlot = T)

