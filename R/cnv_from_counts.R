#' @import ExomeDepth
#' @title Get CNV calls from exome counts
#' @description Using the ExomeDepth algorithm to detect deletions or duplications
#' of regions within the genes of the patient by comparing it against the normal
#' patients (my.ref.samples).
#' @param id number or character that indicates which patient will be analyzed.
#' It has to be the same as its id in the data base.
#' @param bed bedfile that has information about the genes of the mitochondrial
#' genome. It's first 3 columns must be: chromosome, start, end. This file is 
#' given on this package as bedfileMito.RDS in the Data folder.
#' @param ControlDB database in dataframe format that contains the counts of the 
#' genes of normal patients that are used as a reference for making the analysis.
#' @param PatientsDB database in dataframe format that contains the counts of
#' the genes of any patient (normal or not) that has been analyzed.
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' @param transit transition probability: Transition probability of the hidden Markov Chain 
#' from the normal copy number state to either a deletion or a duplication. 
#' The default (0.0001) expect approximately 20 CNVs genome-wide.
#' @return a list: CNV calls, Annot (annotation), cts (exome counts)
#' @examples
#' u=CNVfromCounts(102, DataBase)
#' u=CNVfromCounts(Exons, "~/Patients/Patient102.bam", minoverlap=0.001, 102)
#' @export

cnv_from_counts <- function(id, bed, ControlDB, PatientsDB, minoverlap = 0.0001, transit = 1) {
  
  indice <- which(colnames(PatientsDB) == id)
  my.test <- as.numeric(PatientsDB[, indice, drop = TRUE])
  
  my.ref.samples <- colnames(ControlDB)
  my.reference.selected <- apply(X = ControlDB[, my.ref.samples, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)
  
  #Creates a ExomeDepth object with test and reference:
  all.exons <- new('ExomeDepth',
                   test = my.test ,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')
  
  all.exons <- CallCNVs(x = all.exons, #ExomeDepth object
                        transition.probability = transit,
                        chromosome = bed$chromosome,
                        start = bed$start,
                        end = bed$end,
                        name = bed$name)
  
  exons.GRanges <- GenomicRanges::GRanges(seqnames = bed$chromosome,
                                          IRanges::IRanges(start = bed$start, end = bed$end),
                                          names = bed$name)
  
  #Add names of genes within each variation:
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.GRanges ,
                             min.overlap = minoverlap, #mientras menor es, mas genes
                             column.name = 'name.gene')
  
  CNV_calls <- all.exons@CNV.calls
  
  #Order by size of variation:
  CNV_calls$size <- abs(CNV_calls$start - CNV_calls$end)
  CNV_calls <- CNV_calls[order(CNV_calls$size, decreasing = TRUE),]
  
  Annot <- all.exons@annotations
  return(list(CNV_calls, all.exons, Annot))
  
}