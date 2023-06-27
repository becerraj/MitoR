#Instalaciones----
Instalations<-function(){
  install.packages("BiocManager")
  install.packages("MASS")
  install.packages("Matrix")
  install.packages("spatial")
  #sudo apt -y install libcurl4-openssl-dev
  #Rscript -e "install.packages('RCurl')"

  install.packages("RCurl")
  BiocManager::install("GenomeInfoDb")
  BiocManager::install("Biostrings")
  BiocManager::install("Rsamtools")
  BiocManager::install("GenomicRanges")
  BiocManager::install("GenomicAlignments")
  #BiocManager::install("ExomeDepth")
  install.packages("ExomeDepth", dependencies = TRUE, INSTALL_opts = '--no-lock')
  install.packages("writexl")
  install.packages("openxlsx", dependencies = TRUE, INSTALL_opts = '--no-lock')
  install.packages("readODS")
  install.packages("stringr")

}


Libraries<-function(){
  library(ExomeDepth)
  library(writexl)
  library(Rsamtools)
  library(dplyr)
  library(BiocParallel)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(openxlsx)
  library(ggplot2)
  library(xlsx)
  library(readxl)
  library(readODS)
  library(tibble)
  library(stringr)
  library(dplyr)
  library(ape)
  library(rjson)
  library(httr)
  library(data.table)
  library(Rsamtools)

}


#' @title Initial Data Bases
#' @description Generation of 3 empty databases: Control(normal references),
#'  Patients (analyzed) and Exons (information about the exons detected with variations)
#' @export InitialDB

InitialDB<-function(c, p, e){
  dir.create("~/DataBases")
  dir.create("~/Results")

  #wb1 <- openxlsx::createWorkbook()
  #openxlsx::saveWorkbook(wb1, paste("~/DataBases/",c,".xlsx", sep=""),  overwrite = TRUE)
  wb1<-data.frame()
  saveRDS(wb1, paste("~/DataBases/",c,".RDS", sep=""))

  #wb2 <- openxlsx::createWorkbook()
  #openxlsx::saveWorkbook(wb2, paste("~/DataBases/", p,".xlsx", sep=""),  overwrite = TRUE)
  wb2<-data.frame()
  saveRDS(wb2, paste("~/DataBases/",p,".RDS", sep=""))

  #wb3 <- openxlsx::createWorkbook()
  #openxlsx::saveWorkbook(wb3, paste("~/DataBases/",e,".xlsx"),  overwrite = TRUE)
  wb3<-data.frame()
  saveRDS(wb3, paste("~/DataBases/",e,".RDS", sep=""))

}



#------------------------------------------------------
    #ADD PATIENT TO CONTROL DB
#------------------------------------------------------

#' @title Add patient to Control DB
#' @description Saves the counts of the lectures of each exon in a data base
#' which columns are the ids of the normal patients. It also changes the excel file
#' with the data base saved in the PC
#' @param ControlDB dataframe which contains the counts of normal patients that will be used as reference.
#' @param PatientsDB dataframe which contains the counts of patients that has been analyzed.
#' @return Control DB updated
#' @export AddControl
#' @examples
#' ControlDB<-AddControl(ControlDB, PatientsDB)

AddControl<-function(ControlDB, PatientsDB, id, bed, c, include.chr){
  #Ask for the id of patient

  if (!(id %in% colnames(ControlDB))){ #Make sure the patient is not already in controlDB

    if(id %in% colnames(PatientsDB)){ #If the patients is in PatientsDB-----------------------------
      message("This patient is already in the data base. It will be copied to the control data base")
      indice<-which(colnames(PatientsDB)==id)
      Pcounts <- PatientsDB[,indice]

    }else{ #If the patient is not on the DB------------------------------------------------------------
      path_file<-file.choose()
      if(substr(path_file,nchar(path_file)-3, nchar(path_file)) == ".bam"){ #If its .BAM
        cts <- getBamCounts(bed.frame = bed,
                            bam.files = path_file,
                            include.chr = include.chr,
                            referenceFasta = NULL)
        Pcounts<-data.frame(Pcounts=cts[,ncol(cts)])

      }
      else{ #lo mando a hacer el counts desde FASTQ
        bam.control<-BAMfromFastQ()
        cts <- getBamCounts(bed.frame = bed, #Exons_chrM
                            bam.files = bam.control ,
                            include.chr = include.chr, 
                            referenceFasta = NULL)
        
        Pcounts<-data.frame(Pcounts=cts[,ncol(cts)])
        
      }
    }#---------------------------------------------------------------------------
    
    #update the database
    
    if(max(dim(ControlDB)) == 0){ #this will be the first control patient
      ControlDB<-data.frame(Pcounts)
      names(ControlDB)[names(ControlDB) == "Pcounts"] <- id
      message("Your first control patient has been saved!")

    }else{ #not the first control patient
      ControlDB<-add_column(ControlDB, Pcounts)
      names(ControlDB)[names(ControlDB) == "Pcounts"] <- id
      message("Another control patient has been saved!")
    }

    saveRDS(ControlDB, paste("~/DataBases/",c,".RDS", sep=""))


  }else{ message("This patients is already on your Control Data Base") }

  return(ControlDB)

}


#------------------------------------------------------
    #ADD PATIENT TO PATIENTS DB
#------------------------------------------------------
#' @title Add patient to Patients DB
#' @description Saves the counts of the lectures of each exon in a data base
#' which columns are the ids of the patients.It also changes the excel file
#' with the data base saved in the PC.
#' @param ControlDB dataframe which contains the counts of normal patients that will be used as reference.
#' @param PatientsDB dataframe which contains the counts of patients that has been analyzed.
#' @return Patients DB updated
#' @export AddPatient
#' @examples
#' PatientsDB<-AddPatient(ControlDB, PatientsDB)

AddPatient<-function(ControlDB, PatientsDB, id, bed, p, include.chr){
  #id<-readline(prompt = "Patient ID? ")
  if (id %in% colnames(PatientsDB)){ #If the patients is in Patients DB
    message("This patients is already in your database")
    stop()

  }else if (id %in% colnames(ControlDB)){
    indice<-which(colnames(ControlDB)==id)
    Pcounts <- ControlDB[,indice]
    message("This patient is already in your control database")

  }else{ #If the patient is not on any DB------------------------------------------------------------
    path_file<-file.choose()
    if(substr(path_file,nchar(path_file)-3, nchar(path_file)) == ".bam"){ #If its .BAM
      cts <- getBamCounts(bed.frame = Exons,
                          bam.files = path_file,
                          include.chr = include.chr,
                          referenceFasta = NULL)
      Pcounts<-data.frame(Pcounts=cts[,ncol(cts)])

    }
    else{ #lo mando a hacer el counts desde FASTQ
      bam.control<-BAMfromFastQ()
      cts <- getBamCounts(bed.frame = bed, #Exons_chrM
                          bam.files = bam.control ,
                          include.chr = include.chr, 
                          referenceFasta = NULL)
      
      Pcounts<-data.frame(Pcounts=cts[,ncol(cts)])
    }
  }#---------------------------------------------------------------------------

  #Save counts on PatientsDB

  if(max(dim(PatientsDB)) == 0){ #this will be the first patient in DB
    PatientsDB<-data.frame(Pcounts)
    names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
    message("Your first patient has been saved!")

  }else{ #not the first control patient
    PatientsDB<-add_column(PatientsDB, Pcounts)
    names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
    message("Another patient has been saved!")
  }

  #Save the updated Patients DB
  saveRDS(PatientsDB, paste("~/DataBases/", p, ".RDS", sep=""))

  return(PatientsDB)
}


RemovePatient<-function(DB, id){
  indice<-which(colnames(DB)==id)
  DB <- subset(DB, select=-indice) 
  return(DB)
}


#' @title Get CNV calls from a .BAM file
#' @description Using the ExomeDepth algorithm to detect deletions or duplications
#' of regions within the exomes of the patient by comparing it against the normal
#' patients (my.ref.samples).
#' @param bed.frame dataframe which contains the information about 196,647 exons
#' based on human genome hg19 (chromosome, end, start, gene name, number of exon,
#'  etc.).
#' @param bam_file file in format .bam which contains the previously aligned
#' sequenced data of the patient that will be analyzed.
#' @param DB Data base which contains the exome counts of 14 patients among which
#' there are 2 normal patients which are used as a reference.
#' @param id number or character that indicates which patient is. It is used to
#' detect if this patient has been already analyzed or not. And it is used for
#' the names of output files.
#' @return a list: CNV calls, Annot (annotation), cts (exome counts)
#' @export CNVfromBAM
#' @examples
#' u=CNVfromBAM(Exons, "~/Patients/Patient102.bam", 102)
#' u=CNVfromBAM(Exons, "~/Patients/Patient102.bam", minoverlap=0.001, 102)

CNVfromBAM<-function(bed, bam_file, ControlDB, minoverlap=0.0001, id, include.chr){

  #Exon count from BAM
  cts <- getBamCounts(bed.frame = bed,
                      bam.files = bam_file,
                      include.chr = include.chr, #TENIA QUE SER TRUE para que coincidieran los seqnames con el seqleves(x)
                      referenceFasta = NULL)

  my.test<-cts[,ncol(cts)]

  #Referencies from normal patients:
  #my.ref.samples <- colnames(DB[c(13,14)])


  my.ref.samples <- colnames(ControlDB)
  my.reference.selected <- apply(X = ControlDB[, my.ref.samples, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)


  #CNVs calling:
  all.exons <- new('ExomeDepth',
                   test = my.test ,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')

  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = bed.frame$chromosome,
                        start = bed.frame$start,
                        end = bed.frame$end,
                        name = bed.frame$name)

  exons.GRanges <- GenomicRanges::GRanges(seqnames = bed.frame$chromosome,
                                          IRanges::IRanges(start=bed.frame$start,end=bed.frame$end),
                                          names = bed.frame$name)
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.GRanges ,
                             min.overlap = minoverlap, #mientras menor es, mas genes
                             column.name = 'name.gene')


  CNV_calls<-all.exons@CNV.calls

  #Order by size of variation:
  CNV_calls$size<- abs(CNV_calls$start - CNV_calls$end)
  CNV_calls<-CNV_calls[order(CNV_calls$size, decreasing = TRUE),]

  Annot<-all.exons@annotations
  return(list(CNV_calls, all.exons, Annot, cts))

}


#' @title Get CNV calls from exome counts
#' @description Using the ExomeDepth algorithm to detect deletions or duplications
#' of regions within the exomes of the patient by comparing it against the normal
#' patients (my.ref.samples).
#' @param id number or character that indicates which patient will be analyzed.
#' It has to be the same as its id in the data base. Also it is used for
#' the names of output files.
#' @param bam_file file in format .bam which contains the previously aligned
#' sequenced data of the patient that will be analyzed.
#' @param DB dataframe which contains the exome counts of 14 patients among which
#' there are 2 normal patients which are used as a reference.
#' @param minoverlap is a decimal number which indicates the needed percentage
#' of overlap from which it can be declare that a region has suffered a variation.
#' @return a list: CNV calls, Annot (annotation), cts (exome counts)
#' @export CNVfromCounts
#' @examples
#' u=CNVfromCounts(102, DataBase)
#' u=CNVfromCounts(Exons, "~/Patients/Patient102.bam", minoverlap=0.001, 102)

CNVfromCounts<-function(id, ControlDB, PatientsDB, minoverlap=0.0001){
  indice<-which(colnames(PatientsDB)==id)
  my.test <- PatientsDB[,indice]

  #Reference of normal patients:
  #my.ref.samples <- colnames(DB[c(13,14)])

  my.ref.samples <- colnames(ControlDB)
  my.reference.selected <- apply(X = ControlDB[, my.ref.samples, drop = FALSE],
                                 MAR = 1,
                                 FUN = sum)

  #CNVs calling:
  all.exons <- new('ExomeDepth',
                   test = my.test ,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')

  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = Exons$chromosome,
                        start = Exons$start,
                        end = Exons$end,
                        name = Exons$name)

  exons.GRanges <- GenomicRanges::GRanges(seqnames = Exons$chromosome,
                                          IRanges::IRanges(start=Exons$start,end=Exons$end),
                                          names = Exons$name)
  all.exons <- AnnotateExtra(x = all.exons,
                             reference.annotation = exons.GRanges ,
                             min.overlap = minoverlap, #mientras menor es, mas genes
                             column.name = 'name.gene')

  CNV_calls<-all.exons@CNV.calls

  #Order by size of variation:
  CNV_calls$size<- abs(CNV_calls$start - CNV_calls$end)
  CNV_calls<-CNV_calls[order(CNV_calls$size, decreasing = TRUE),]

  Annot<-all.exons@annotations
  return(list(CNV_calls, all.exons, Annot))

}



#-------------------------------------------------
    #ANALYZE PATIENT
#--------------------------------------------


#' @title Analyze the CNVs of the patient
#' @description Apply CNVsfromBAM or CNVsfromCounts depending from where the patient is taken from.
#' @param ControlDB dataframe which contains the counts of normal patients that will be used as reference.
#' @param PatientsDB dataframe which contains the counts of patients that has been analyzed.
#' @return a list: Patients DB (updated if the patient that has been analyzed is new),
#' CNV calls (dataframe with information about CNVs),
#' all.exons (large exomedepth)
#' @export AnalyzePatient
#' @examples
#' outAnalisis<-AnalyzePatient(ControlDB, PatientsDB)
#' PatientsDB<-outAnalisis[[1]]
#' CNV_calls<-outAnalisis[[2]]
#' all.exons<-outAnalisis[[3]]


AnalyzePatient<-function(ControlDB, PatientsDB, id, bed, p, include.chr){
  #id<-readline(prompt = "Patient ID? ")

  if(id %in% colnames(ControlDB)){ # I cant use it as a reference for itself
    indice<-which(colnames(ControlDB)==id)
    Pcounts<-ControlDB[,indice]


    #Delete the patient from ControlDB
    ControlDB <- subset(ControlDB, select = -indice)
    if(max(dim(ControlDB))==0 | ncol(ControlDB)==0){ stop("You need to have at least one control patient") }

    #Add the patient to PatientsDB if its not there----
    if(!(id %in% colnames(PatientsDB))){
      if(max(dim(PatientsDB)) == 0){ #this will be the first patient in DB
        PatientsDB<-data.frame(Pcounts)
        names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
        message("Your first patient has been saved!")

      }else{ #not the first control patient
        PatientsDB<-add_column(PatientsDB, Pcounts)
        names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
        message("Another patient has been saved!")
      }
    }
  }

  print(" ")

  if (id %in% colnames(PatientsDB)){ #If the patients is in Patients DB---------------
    message("This patients is already on your database.The process will take less than a minute")
    #CNVfromCounts
    out_counts<-CNVfromCounts(id, ControlDB, PatientsDB)
    CNV_calls<-out_counts[[1]]
    all.exons<-out_counts[[2]]
    #Annot<- out_counts[[3]]


  }else{ #If its a new patient---------------------------------------------------------

    print("nuevo paciente")
    message("This patients is not already on your database.The process may take a few minutes")

    path_file<-file.choose()
    if(substr(path_file,nchar(path_file)-3, nchar(path_file)) == ".bam"){
      bam_file<-path_file
      out_counts<-CNVfromBAM(bed, bam_file, ControlDB, minoverlap=0.0001, id, include.chr)
      cts<-out_counts[[4]]
      Pcounts<-data.frame(Pcounts=cts[,ncol(cts)])
      CNV_calls<-out_counts[[1]]
      all.exons<-out_counts[[2]]

      }else{ #lo mando a hacer el counts desde FASTQ ESTO ADELANTE
        bam.control<-BAMfromFastQ()
        cts <- getBamCounts(bed.frame = bed, #Exons_chrM
                            bam.files = bam.control ,
                            include.chr = include.chr, 
                            referenceFasta = NULL)
        
        Pcounts<-data.frame(Pcounts=cts[,ncol(cts)])

      }

    #I will add the new patient to the PatientsDB
    if(max(dim(PatientsDB)) == 0){ #this will be the first patient in DB
      PatientsDB<-data.frame(Pcounts)
      names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
      message("Your first patient has been saved!")

      }else{ #not the first control patient ESTO ADELANTE
      PatientsDB<-add_column(PatientsDB, Pcounts)
      names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id
      message("Another patient has been saved!")
      }

    saveRDS(PatientsDB, paste("~/DataBases/", p, ".RDS", sep=""))


    }#--------------------------------------------------------------

  print("aca final")
  #CNV_calls<-out_counts[[1]]
  #all.exons<-out_counts[[2]]

  return(list(PatientsDB,  CNV_calls, all.exons))
}


#' @title Order name.gene column of CNV_calls
#' @description present the exons in a more eficient way to read them.
#' Example: GENE1_2, GENE1_5 will be transformed to GENE1(2,5)
#' @param CNV_file dataframe with cnv calls information.
#' @return CNV file: is the cnv calls ordered.
#' @export OrderCNV
#' @examples
#' cnv_call=OrderCNV(cnv_calls)

OrderCNV <- function(CNV_file){

  CNV_file$genes<- strsplit(CNV_file$name.gene, split=",")

  atras<-c()
  outs<-c()
  yala<-c()
  count<-1

  for (fila in CNV_file$genes){
    y<-strsplit(unlist(fila), split="_")
    for (k in y){
      atras<-k[2] #save the number of the exon
      if (!(k[1] %in% yala)){ #if the name of the gene is not already counted
        for (j in y){
          if ( k[1] == j[[1]] ){
            atras<-append(atras, j[2])
          }
        }
        yala<-append(yala,k[1])
        atras<-list(unique(atras))
        out<-paste(k[1], "(", sapply(atras, paste, collapse=","), ")")
        outs<-append(outs, out)
        atras<-c()
        out<-c()
      }
    }
    if(length(outs)>1){
      outs<-paste(unlist(outs), collapse=" ")
    }
    #print(outs)
    CNV_file$genes[count]<-outs

     #count the number of genes
    CNV_file$nGenes[count]<-length(yala)

    count<-count+1
    #print(outs)
    outs<-c()
    yala<-c()

  }

  #mando al final la info menos importante
  CNV_file <- subset(CNV_file, select=c("chromosome", "type", "size", "id", "genes", "nGenes", "nexons","BF", "reads.expected", "reads.observed", "reads.ratio"))
  colnames(CNV_file)[1]<-"chr"
  return(CNV_file)

}




#' GetCNVsAnnotation
#' @param chr (chromosome) any of chr1..chr22
#' @param x a CNV object from ExonDepth (all.exons)
#' @param countThreshold (default = 10). The minimum expected exon count
#' @return an annotation data frame for the specified chromosome to be used for plotting
#' @export GetCNVsAnnotation

GetCNVsAnnotation <- function(chr, x, countThreshold =10){
  anno <- x@annotations
  selected <- which(anno$chromosome == chr & (x@test + x@reference) *
                      x@expected > countThreshold)

  if (length(selected) == 0) {
    warning("No exon seem to be located in the requested region. It could also be that the read count is too low for these exons? In any case no graph will be plotted.")
    return(NULL)
  }
  anno <- anno[selected, ]
  anno$expected <- x@expected[selected]
  anno$freq <- x@test[selected]/(x@reference[selected] +
                                   x@test[selected])
  anno$middle <- 0.5 * (anno$start + anno$end)
  anno$ratio <- anno$freq/anno$expected
  anno$test <- x@test[selected]
  anno$reference <- x@reference[selected]
  anno$total.counts <- anno$test + anno$reference

  if (length(x@phi) == 1)
    anno$phi <- x@phi
  else anno$phi <- x@phi[selected]
  anno <- cbind(anno,plyr::ldply(1:nrow(anno), function(i){
    c(my.min.norm.prop=qbetabinom(p = 0.025, size = anno$total.counts[i],
                                  phi = anno$phi[i], prob = anno$expected[i]),
      my.max.norm.prop=qbetabinom(p = 0.975, size = anno$total.counts[i],
                                  phi = anno$phi[i], prob = anno$expected[i]) )
  }))

  anno$my.min.norm.prop <- anno$my.min.norm/anno$total.counts
  anno$my.max.norm.prop <- anno$my.max.norm/anno$total.counts

  CNVsret <- subset(x@CNV.calls, chromosome == chr)
  attr(CNVsret,"BamHeader") <- attr(x@CNV.calls,"BamHeader")
  return(list(Exons=anno, CNVs=CNVsret))
}



#' @title Data needed to plot 1 CHR
#' @description Generates the data needed to plot 1 chromosome
#' @param annot output from GetCNVsAnnotation function (exons and cnvs)
#' @return puntos: exons of the hole chromosome
#' lp: exons of the left arm of the chromosome
#' lq: exons of the rigth arm
#' CNVs (variations of the chromosome)
#' @export .PlotCNVchromosome
#' @examples

.PlotCNVchromosome<-function(annot,thLength=500){
  #load("~/CNVAnalisis-Package/CNVAnalisis/Data-CNVs/Centromere.RData") #indica donde estan los centros de cada brazo del cromosoma

  CNVs <- annot$CNVs #cnvs de 1 chr
  CNVs$width <- (CNVs$end-CNVs$start)
  chr <- stringr::str_replace_all(annot$Exons$chromosome[1],"chr","Chr")
  CNVs <- subset(CNVs,width > thLength) #region con variacion >500

  annot <- annot$Exons #exons de 1 chr
  ordm <- order(annot$middle)
    #middle <- 0.5 * (anno$start + anno$end) = centro de cada exon

  puntos <- data.frame(x=annot$middle[ordm], y=annot$ratio[ordm])
    #x=  centro cada exon ordenado segun middle
    #y= frq/expected ordenado
  puntos$colscale <- puntos$y
  puntos$colscale[puntos$colscale>2] <- 2
    #normalizo el maximo=2


  chr <- unique(annot$chromosome)
  print(chr)
  if(nrow(subset(puntos,x < Centromere$hg19[Centromere$hg19$chr==chr,]$left))>0){
    #Centromere[[Genome]] = left and rigth from each chr
    #if: si hay algun puntos$x < left del chr selec: osea si hay algun punto del brazo izq
    lp <- loess(y~x,subset(puntos,x < Centromere$hg19[Centromere$hg19$chr==chr,]$left), span = 0.5)
      #loess:makes the curve
      #span= grade of smothing
      #y,x predictors
      #data: dots of the left arm
  }else{
    lp<-data.frame(x=puntos$x[1],y=puntos$y[1],fitted=puntos$x[1])
    #si no hay puntos en el brazo izq, tomo el primero que aparezca. CREO que esto està mal. No puedo no ponerlo porque no me puede quedar vacio en el return
  }

  #agrego el if
  if(nrow(subset(puntos,x > Centromere$hg19[Centromere$hg19$chr==chr,]$right))>0){
    lq <- loess(y~x,subset(puntos,x > Centromere$hg19[Centromere$hg19$chr==chr,]$right), span = 0.5)
    #polinomio de  los puntos del brazo derecho
  }else{
    lq<-data.frame(x=puntos$x[nrow(puntos)],y=puntos$y[nrow(puntos)],fitted=puntos$x[nrow(puntos)])
  }

  return(list(puntos, lp, lq, CNVs))
}



#' @title Data needed for complete visualization
#' @description Generates the data needed to plot all chromosomes of a patient
#' Here are included the previous functions GetCNVsAnnotation and .PlotCNV
#' @param all.exons large ExomeDepth file
#' @return Todos_puntos: all exons from all the chromosomes ordered
#' Todos_lp: exons of the left arm from all the chromosomes ordered
#' Todos_lq: exons of the right arm from all the chromosomes ordered
#' Todos_cnvs: all the variations of the patient
#' div_chr: division of chromosomes
#' cent: middle of arms in every chromosome
#'
#' all the outputs have the x axis data ordered by chromosome
#' @export DataforVisualization
#' @examples
#'

DataforVisualization<-function(all.exons){

  #Genero puntos y brazos:-----
  list_puntos<-c()
  list_lp<-c()
  list_lq<-c()
  list_cnvs<-c()

  for (i in 1:22){
    obj <- GetCNVsAnnotation(chr = i  , x=all.exons)
    out<-.PlotCNVchromosome(obj, thLength=500)
    puntos<-out[[1]]
    puntos$chr<-i
    list_puntos[[i]]<-puntos

    lp<-out[[2]]
    lp$chr<-i
    lp<-data.frame(x=lp$x, fitted=lp$fitted, chr=lp$chr)
    list_lp[[i]]<-lp

    lq<-out[[3]]
    lq$chr<-i
    lq<-data.frame(x=lq$x, fitted=lq$fitted, chr=lq$chr)
    list_lq[[i]]<-lq

    cnvs<-out[[4]]
    cnvs<-data.frame(start= cnvs$start, end=cnvs$end, reads.ratio=cnvs$reads.ratio, chr=cnvs$chromosome)
    list_cnvs[[i]]<-cnvs

  }

  #--------

  #Unifico todos los puntos y brazos:---
  Todos_puntos<-data.frame()
  Todos_lp<-data.frame()
  Todos_lq<-data.frame()
  Todos_cnvs<-data.frame()

  for(i in 1:length(list_puntos)){
    Todos_puntos<-rbind(list_puntos[[i]], Todos_puntos)
    Todos_lp<-rbind(list_lp[[i]], Todos_lp)
    Todos_lq<-rbind(list_lq[[i]], Todos_lq)
    Todos_cnvs<-rbind(list_cnvs[[i]], Todos_cnvs)
  }

  #--------

  #Para que queden los cromosomas en order en Todos_puntos: ----
  div_chr<-c()
  final<-c()
  Todos_puntos$x<-as.numeric(Todos_puntos$x)
  Todos_cnvs$start<-as.numeric(Todos_cnvs$start)
  Todos_cnvs$end<-as.numeric(Todos_cnvs$end)
  ch<-sort(unique(Todos_puntos$chr))

  cent<-Centromere$hg19[1:22,] ##

  for (i in 1:length(ch)){
    print(ch[i])
    final<- max(Todos_puntos$x[Todos_puntos$chr==ch[i]])
    finall<-final + 1000000
    div_chr<-append(div_chr, finall)

    if (i!=length(ch)){
      Todos_puntos$x[Todos_puntos$chr == ch[i+1]]<- Todos_puntos$x[Todos_puntos$chr == ch[i+1]] + (final+100000000)
      Todos_lp$x[Todos_lp$chr == ch[i+1]]<- Todos_lp$x[Todos_lp$chr == ch[i+1]] + (final+100000000)
      Todos_lq$x[Todos_lq$chr == ch[i+1]]<- Todos_lq$x[Todos_lq$chr == ch[i+1]] + (final+100000000)
      Todos_cnvs$start[Todos_cnvs$chr==ch[i+1]]<-Todos_cnvs$start[Todos_cnvs$chr == ch[i+1]] + (final+100000000)
      Todos_cnvs$end[Todos_cnvs$chr==ch[i+1]]<-Todos_cnvs$end[Todos_cnvs$chr == ch[i+1]] + (final+100000000)

      cent$left[cent$chr==ch[i+1]]<-cent$left[cent$chr==ch[i+1]] + (final+100000000)##
      cent$right[cent$chr==ch[i+1]]<-cent$right[cent$chr==ch[i+1]] + (final+100000000)##

    }
  }
  return(list(Todos_puntos,Todos_lp, Todos_lq, Todos_cnvs, div_chr, cent))
}



#' @title Generation of complete graph
#' @description Generates the plot with all chromosomes of a patient.
#' @param Todos_puntos
#' @param Todos_lp
#' @param Todos_lq
#' @param div_chr
#' @param cent
#' all parameters are outputs from the function DataforVisualizatio
#' @return p the plot
#' @export GenerateGraph
#' @examples
#'

GenerateGraph<-function(Todos_puntos,Todos_lp, Todos_lq, div_chr, cent){
  #1. Formaciòn del fondo de exones:
      #pinta 3 colores segun valores de y (breaks):
  p <- ggplot(Todos_puntos, aes(x,y,colour=colscale)) + geom_point(size=0.5) +
    scale_color_stepsn(colours = c("red", "white","blue"),breaks = c(  0.9, 1.1 ) )+
    geom_vline(xintercept=div_chr, linetype=3)+

      #titulos
    labs(x="base position",y="Ratio") +
      theme(legend.position = "none",plot.title = element_text(margin = margin(t=40,b=-30)))


  #2. 2 curvas de los brazos:
  p+ geom_line(aes(x=x,y=y),data=data.frame(x=Todos_lp$x,y=Todos_lp$fitted),colour="red")+
    geom_line(aes(x=x,y=y),data=data.frame(x=Todos_lq$x,y=Todos_lq$fitted),colour="pink")+

    theme(axis.text.x = element_text(angle=45, vjust=.5, hjust=1, size=0.5)) +
    ylim(-.1,3.5) +
    geom_hline(yintercept = 1)


  #Otra manera de diferenciar los brazos puede ser poniendo linea de puntos en cada chr:
  p <-p+ geom_vline(xintercept=abs((cent$left+cent$right)/2), linetype=3, colour="orange")


  #Pinta de colores azul o rojo las variaciones segun el eje y (> o < 1)
  if(nrow(Todos_cnvs)>0){
    if(any(Todos_cnvs$reads.ratio>1)){
      print("aca")
      p <- p + geom_segment(data=subset(Todos_cnvs,reads.ratio>1), aes(x = start, y = reads.ratio-0.01, xend = end, yend = reads.ratio+0.01), col="green", size=2)
    }
    if(any(Todos_cnvs$reads.ratio<1)){
      print("aca2")
      p <- p + geom_segment(data=subset(Todos_cnvs,reads.ratio<1), aes(x = start, y = reads.ratio-0.01, xend = end, yend = reads.ratio+0.01), col="green", size=2)
    }
  }

  #Indicadores de cromosoma
  indicator<-data.frame(x=1:22,y=0, chr="chr ")
  indicator$x[1]<-div_chr[1]/2
  indicator$chr[1]<-"chr1"
  for(i in 2:22){
    indicator$x[i]<-(div_chr[i]+div_chr[i-1])/2
    indicator$chr[i]<-paste("chr", i, sep="")
  }

  p<- p + annotate("text", x=indicator$x, y=-0.2, label= indicator$chr, size= 4, angle=90)

  return(p)

}

#Lo que sigue no se si hacerlo o no

#marca los genes


#' @title Generation of Links
#' @description Generates links to databases for each variation
#' @param CNV_calls
#' @return CNV_calls with 2 new columns "Franklin" and "Varsome"
#' @export GenerateLinks

GenerateLinks<-function(CNV_calls){

  for (i in 1:nrow(CNV_calls)){
    mutF<-CNV_calls$id[i]
    mutF<- gsub(':', '-', mutF)
    mutV<-gsub("-", "%3A", mutF)
    addressF<-paste("https://franklin.genoox.com/clinical-db/variant/sv/",mutF, sep="")
    addressV<-paste("https://varsome.com/cnv/hg19/", mutV, sep="")

    if(CNV_calls$type[1] == "deletion"){
      addressF<-paste(addressF, "-DEL", sep="")
      addressV<-paste(addressV, "%3ADEL?", sep="")
    }else{
      addressF<-paste(addressF, "-DUP", sep="")
      addressV<-paste(addressV, "%3ADUP?", sep="")
    }

    CNV_calls$Franklin[i] <-addressF
    CNV_calls$Varsome[i] <-addressV

  }

  CNV_calls<-CNV_calls[,c(1,2,3,4,5,6,7,12,13,8,9,10,11)]
  return(CNV_calls)
}



# Mostramos el dataframe
#df<-datatable(df, escape = FALSE)


#ruta_completa <- normalizePath("mi_archivo.xlsx")


#' @title Generation of report
#' @description Generates the excel with the results of the patient's analisis
#' @param CNV_calls dataframe with the information of the mutations
#' @param graph ggplot
#' @param id
#' @return workbook
#' @export Report

Report<- function(CNV_calls, graph, id){

  wb <- openxlsx::createWorkbook()

  COLNAMES_STYLE <- createStyle(
    fontSize = 12,
    textDecoration = "bold",
    halign = "center", valign = "center", border = "TopBottom",
    borderColour = "black")


  hyperlink_style <- openxlsx::createStyle(fontColour = "#0000FF")

  addWorksheet(wb, paste(id, "CNVAnalisis"))
  addWorksheet(wb, paste(id, "CNVplots"))
  writeData(wb, sheet = 1, x = CNV_calls,
            headerStyle = COLNAMES_STYLE,
            borderStyle = "dashed",
            borders = "columns", borderColour = "black")


  setColWidths(wb, sheet = 1, cols = 1:4, widths = "auto")
  setColWidths(wb, sheet = 1, cols = 6:ncol(CNV_calls), widths = "auto")
  setColWidths(wb, sheet = 1, cols = 5, widths = 16)

  for(i in 1:nrow(CNV_calls)){
    x<-CNV_calls$Franklin[i]
    y<-CNV_calls$Varsome[i]
    names(x) <- c("View on Franklin DB")
    names(y) <- c("View on Varsome DB")
    class(x) <- "hyperlink"
    class(y) <- "hyperlink"

    writeData(wb, sheet = 1, x = x, startRow= i+1, startCol = 8)
    writeData(wb, sheet = 1, x = y, startRow= i+1, startCol = 9)

  }

  #Colour cells if the size of the variation is >500

  for (i in 1:nrow(CNV_calls)){
    if (CNV_calls$size[i] > 500){
      addStyle(wb, 1, style = openxlsx::createStyle(fgFill = "thistle"), rows = i+1, cols = 1:ncol(CNV_calls))
    }
  }

  #Add Plot
  print(graph_)
  insertPlot(wb, 2, xy = c(1, 1), width = 25, height = 18, fileType = "png", units = "cm")

  #saveWorkbook(wb, paste("~/CNVAnalisis-Package/CNVAnalisis/Results/",id,"_report.xlsx", sep=""))
  openxlsx::saveWorkbook(wb, paste("~/Results/",id,"_report.xlsx", sep=""),  overwrite = TRUE)


}



#' @title Add another patient's information to Exons DB
#' @description Adds the information about the CNVs of exons of a new patient
#' and keeps track of the frequency and type of variation of each exon.
#' @param CNV calls dataframe which contains the information about the CNVs of the new patient
#' @param id of the patient
#' @return Exons DB updated
#' @export UpdateExonsDB
#' @examples
#' ExonsDB<-UpdateExonsDB(CNV_calls)

UpdateExonsDB<-function(CNV_calls, id, e){
  x<-c()
  exones_var<-c()
  for (i in CNV_calls$name.gene){
    print(i)
    x<-append(x,strsplit(i, split=","))
  }
  for (l in 1:length(x)){
    for (g in 1:length(x[[l]])){
      exones_var<-append(exones_var, x[[l]][g])
    }
  }

  Completo<- data.frame("Exon"=exones_var)
  Completo$Gen<- sub("_.*", "", Completo$Exon)

  #Poner la cantidad de duplicaciones o deleciones por cada gen
  todos_dup<-c()
  todos_del<-c()

  #Lista de los exones duplicados y otra con los deleteados:
  exones_duplicados<- unlist(CNV_calls[CNV_calls$type == "duplication", "name.gene" ])

  if(!(length(exones_duplicados)==0)){
    todos_dup<-c()
    q<-c()
    for (i in exones_duplicados){
      q<-append(q, strsplit(i, split=","))
    }
    for (l in 1:length(q)){
      for (g in 1:length(q[[l]])){
        todos_dup<-append(todos_dup, q[[l]][g])
      }
    }
  }

  exones_delecionados<- unlist(CNV_calls[CNV_calls$type == "deletion", "name.gene" ])
  if(!(length(exones_delecionados)==0)){
    todos_del<-c()
    w<-c()
    for (i in exones_delecionados){
      w<-append(w,strsplit(i, split=","))
    }
    for (l in 1:length(w)){
      for (g in 1:length(w[[l]])){
        todos_del<-append(todos_del, w[[l]][g])
      }
    }
  }

  #Para cada paciente(columna): marco si el exon (fila) fue duplicado o deleteado:
  Completo$Paciente[Completo$Exon %in% todos_dup ]<- "duplication"
  Completo$Paciente[Completo$Exon %in% todos_del ]<- "deletion"



  if(max(dim(ExonsDB))==0){#If its the first patient in Exons DB
    ExonsDB<-data.frame(Completo)
  }else{
    p<-as.data.frame(Completo[,c(1,3)])
    ExonsDB<-merge(ExonsDB, p, by="Exon", all=T)
  }
  names(ExonsDB)[names(ExonsDB) == "Paciente"] <- id
  ExonsDB$Gen<- sub("_.*", "", ExonsDB$Exon)

  if("Frec_dup" %in% colnames(ExonsDB)){
    ExonsDB<-select(ExonsDB, -c("Frec_dup", "Frec_del", "Frec_mutac"))
  }

  ExonsDB$Frec_dup<-0
  ExonsDB$Frec_del<-0

  for(i in 1:nrow(ExonsDB)){
    for (j in 3:(ncol(ExonsDB)-2)){
      if(!(is.na(ExonsDB[i,j])) & ExonsDB[i,j]=="deletion" ){
        ExonsDB$Frec_del[i]<-ExonsDB$Frec_del[i]+1
      }else if(!(is.na(ExonsDB[i,j])) & ExonsDB[i,j]=="duplication"){
        ExonsDB$Frec_dup[i]<-ExonsDB$Frec_dup[i]+1
      }
    }
  }


  #Total de mutaciones por exon:
  ExonsDB$Frec_mutac<- ExonsDB$Frec_dup + ExonsDB$Frec_del


  #Save the updated Exons DB
  print("Updating the db")

  saveRDS(ExonsDB, paste("~/DataBases/", e, ".RDS", sep=""))

  return(ExonsDB)

}

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#                 FASTQ
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------



#DOWNLOADS----
# BWA
downloadBWA <- function(){
  createFile("~/MitoRSoftware", "BWA")

  setwd("~/MitoRSoftware/BWA")
  system2("wget", "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/x86_64/bwa-0.7.17-lp154.6.1.x86_64.rpm", wait=TRUE)
  system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.x86_64.rpm | cpio -id", wait=TRUE)


  createFile("~/MitoRSoftware/BWA", "bwa-0.7.17-lp154.6.1.src")
  setwd("~/MitoRSoftware/BWA/bwa-0.7.17-lp154.6.1.src")
  system2("wget", "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/src/bwa-0.7.17-lp154.6.1.src.rpm", wait=TRUE)
  system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.src.rpm | cpio -i", wait=TRUE)

  setwd("~/MitoRSoftware/BWA")
  system2("wget", "https://download.opensuse.org/repositories/home:/vojtaeus/15.4/i586/bwa-0.7.17-lp154.6.1.i586.rpm", wait=TRUE)
  system2("rpm2cpio", "bwa-0.7.17-lp154.6.1.i586.rpm | cpio -id", wait=TRUE)
  # No entiendo porque no entiende el valor de ~ dentro del directorio
  setwd("~")
  parche <- getwd()
  return(sprintf('%s/MitoRSoftware/BWA/usr/bin/bwa', parche))
}

# GATK
downloadGATK <- function(){
  setwd("~/MitoRSoftware")
  URL <- 'https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip'
  system2("wget", URL, wait = TRUE)
  file <- basename(URL)

  # Decompression
  filedc <- substr(file,start=0,stop=(nchar(file)-4))
  system2("unzip", sprintf(file), wait = TRUE)

  setwd("~")
  parche <- getwd()
  setwd(sprintf("~/MitoRSoftware/%s", filedc))
  return(sprintf('%s/MitoRSoftware/%s/gatk-package-4.3.0.0-local.jar', parche, filedc))
}

# PICARD
downloadPICARD <- function(){
  setwd("~/MitoRSoftware")
  system2("wget", "https://github.com/broadinstitute/picard/archive/refs/tags/2.27.5.tar.gz")
  system2("gzip" , "-d ~/MitoRSoftware/2.27.5.tar.gz")
  system2("tar" , "-xvf ~/MitoRSoftware/2.27.5.tar")
  setwd("~/MitoRSoftware/picard-2.27.5")
  system2("wget", "https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar")

  setwd("~")
  parche <- getwd()
  return(sprintf('%s/MitoRSoftware/picard-2.27.5/picard.jar', parche))
}

# REFERENCIAS
downloadHG38 <- function(){
  createFile("~/MitoRSoftware", "RefHG38")
  setwd("~/MitoRSoftware")
  setwd("~/MitoRSoftware/RefHG38")
  URL <- "https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"
  referencia <- substr(basename(URL),1,nchar(basename(URL))-3)
  system2("wget", URL, wait = TRUE)
  system2("gzip" , "-d Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz")
  referenciaN <- str_replace(referencia, "fa", "fasta")
  system2("mv", sprintf("Homo_sapiens.GRCh38.dna.chromosome.MT.fa %s", referenciaN))

  system2("bwa", sprintf("index -a is %s", referenciaN))
  system2("java", sprintf("-jar %s CreateSequenceDictionary -R %s", GATK, referenciaN))
  system2("samtools", sprintf("faidx %s", referenciaN))

  setwd("~")
  parche <- getwd()
  return(sprintf("%s/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta", parche))
}

downloadHG19 <- function(){
  createFile("~/MitoRSoftware", "RefHG19")
  setwd("~/MitoRSoftware/RefHG19")
  URL <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz"
  referencia <- substr(basename(URL),1,nchar(basename(URL))-3)
  system2("wget", URL, wait = TRUE)
  system2("gzip" , sprintf("-d %s", basename(URL)))
  referenciaN <- str_replace(referencia, "fa", "fasta")
  system2("mv", sprintf("%s %s", referencia, referenciaN))
  system2(BWA, sprintf("index -a is %s", referenciaN))
  system2("java", sprintf("-jar %s CreateSequenceDictionary -R %s", GATK, referenciaN))
  system2("samtools", sprintf("faidx %s", referenciaN))

  setwd("~")
  parche <- getwd()
  return(sprintf("%s/MitoRSoftware/RefHG19/%s", parche, referenciaN))
}

downloadSamtools <- function(){
  dir.create("~/MitoRSoftware/Samtools")
  setwd("~")
  setwd("~/MitoRSoftware/Samtools")
  URL <- "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2"
  system2("wget", URL, wait = TRUE)
  system2("bzip2", c("-d", basename(URL)), wait = TRUE)
  system2("tar" , c("-xvf", sprintf("~/MitoRSoftware/Samtools/%s",  list.files("~/MitoRSoftware/Samtools"))))
  file.remove(list.files("~/MitoRSoftware/Samtools")[stringr::str_detect(list.files("~/MitoRSoftware/Samtools"), ".tar")])
  setwd(sprintf("%s", list.files("~/MitoRSoftware/Samtools")))
  system2("./configure")
  system2("make")

  setwd("~")
  parche <- getwd()
  return(sprintf("%s/MitoRSoftware/Samtools/samtools-1.16.1/samtools", parche))
}




# FUNCIONES ----
existingFiles <- function(path){
  setwd(path)
  files <- as.list(system2("ls", stdout=TRUE))
  return(files)
}

createFile <- function(path, name){
  files <- existingFiles(path)
  if (name %in% files){
    setwd(sprintf("%s/%s", path, name))
  } else{
    system2("mkdir", sprintf('%s/%s', path, name))
    setwd(sprintf('%s/%s', path, name))
  }
}


downloadMitoRSoftwares <-function(){
  tryCatch(
    expr = {
      if(file.exists("~/MitoRSoftware/gatk-4.3.0.0") & "gatk" %in% existingFiles("~/MitoRSoftware/gatk-4.3.0.0")){
        message("GATK is already installed")
      }else{
        GATK <<- downloadGATK()
        system2("java", c(sprintf("-jar %s", GATK), "-h"))
        system2("rm", "~/MitoRSoftware/gatk-4.3.0.0.zip")
      }

    },
    error = function(e){
      message("An error occured while perfomring the GATK download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install zip gzip tar
------------------------------------------------------")
      print(e)
    },
    warning = function(w){
      message("An error occured while perfomring the GATK download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
          sudo apt install zip gzip tar rpm2cpio
------------------------------------------------------")
      print(w)
    },
    finally = {
      message("-.Message from GATK")
    }
  )

  tryCatch(
    expr = {
      if("usr" %in% existingFiles("~/MitoRSoftware/BWA")){
        message("BWA is already installed")
      }else{
      BWA <<- downloadBWA()
      system2(BWA)
      system2("rm", "~/MitoRSoftware/bwa-0.7.17.tar")
      }
    },
    error = function(e){
      message("An error occured while perfomring the BWA download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install zip gzip tar
------------------------------------------------------")
      print(e)
    },
    warning = function(w){
      message("An error occured while perfomring the BWA download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install zip gzip tar
------------------------------------------------------")
      print(w)
    },
    finally = {
      message("-.Message from BWA")
    }
  )

  tryCatch(
    expr = {
      if("picard.jar" %in% existingFiles("~/MitoRSoftware/picard-2.27.5")){
        message("Picard is already installed")
      }else{
      PICARD <<- downloadPICARD()
      system2("java", sprintf("-jar %s -h", PICARD))
      }
    },
    error = function(e){
      message("An error occured while perfomring the PICARD download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install git
------------------------------------------------------
You also must register your git account to this local computer. It will allow the cloning
              of PICARD's repository")
      print(e)
    },
    warning = function(w){
      message("An error occured while perfomring the PICARD download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install git
------------------------------------------------------
You also must register your git account to this local computer. It will allow the cloning
              of PICARD's repository")
      print(w)
    },
    finally = {
      message("-.Message from PICARD")
    }
  )

  tryCatch(
    expr = {
      if("Homo_sapiens.GRCh38.dna.chromosome.MT.fasta" %in% existingFiles("~/MitoRSoftware/RefHG38")){
        message("Reference hg38 is already downloaded")
      }else{
      #referencia <<- downloadHG38()
      system2((sprintf("%s", BWA)), sprintf("index -a is %s", referencia))
      }
    },
    error = function(e){
      message("An error occured while perfomring the HG38 download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install zip gzip tar
------------------------------------------------------")
      print(e)
    },
    warning = function(w){
      message("An error occured while perfomring the HG38 download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install zip gzip tar
------------------------------------------------------")
      print(w)
    },
    finally = {
      message("-.Message from BWA")
    }
  )

  tryCatch(
    expr = {
      referencia19 <<- downloadHG19()
      #system2((sprintf("%s", BWA)), sprintf("index -a is %s", referencia19))
    },
    error = function(e){
      message("An error occured while perfomring the HG19 download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install zip gzip tar
------------------------------------------------------")
      print(e)
    },
    warning = function(w){
      message("An error occured while perfomring the HG19 download and installation.
              Please remember that some packages are required and you must download them by yourself.
              On the Linux command-line print:
                ------------------------------------------------------
                sudo apt install zip gzip tar
              ------------------------------------------------------ ")
      print(w)
    },
    finally = {
      message("-.Message from NCBI")
    }
  )

  tryCatch(
    expr = {
      SAMTOOLS <<- downloadSamtools()
      system2(SAMTOOLS)
      system2("rm", "~/MitoRSoftware/Samtools/samtools-1.16.1/samtools")
    },
    error = function(e){
      message("An error occured while perfomring the GATK download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
            sudo apt install zip gzip tar
------------------------------------------------------")
      print(e)
    },
    warning = function(w){
      message("An error occured while perfomring the GATK download and installation.
      Please remember that some packages are required and you must download them by yourself.
On the Linux command-line print:
------------------------------------------------------
          sudo apt install zip gzip tar rpm2cpio
------------------------------------------------------")
      print(w)
    },
    finally = {
      message("-.Message from GATK")
    }
  )
}


searchPatient <- function(path){
  files <- unlist(existingFiles(path))
  pacienteR1 <- files[str_detect(files, "R1")]
  pacienteR2 <- files[str_detect(files, "R2")]
  paciente <- unlist(str_split(substr(str_extract(pacienteR1, pattern = "^(.+?)R1"), 1, nchar(str_extract(pacienteR1, pattern = "^(.+?)R1"))-2), ""))

  if (!(grepl('^[A-Za-z0-9]+$', paciente[length(paciente)], ignore.case = TRUE))){
    paciente <- paste(paciente[-(length(paciente))], sep = "", collapse = "")
  }
  if (is.null(paciente) | paciente == ""){
    paciente <- basename(path) # Llegado el caso que el archivo se llame solamente R1 o R2
  }
  return(list(pacienteR1, pacienteR2, paciente))
}


InstallationsforFASTQ<- function(){
  if (!("MitoRSoftware" %in% existingFiles("~"))){
  createFile("~", "MitoRSoftware")
  setwd("~/MitoRSoftware")
  downloadMitoRSoftwares()

  } else{
    tryCatch(expr = {
      setwd("~")
      parche <- getwd()
      BWA <- sprintf("%s/MitoRSoftware/BWA/usr/bin/bwa", parche)
      PICARD <- sprintf("%s/MitoRSoftware/picard-2.27.5/picard.jar", parche)
      GATK <-sprintf("%s/MitoRSoftware/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar", parche)
      referenciaHG38 <- "~/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta"
      referenciaHG19 <- "~/MitoRSoftware/RefHG19/chrM.fasta"
      system2(sprintf("%s", BWA))
      system2("java", sprintf("-jar %s -h", GATK))
      system2("java", sprintf("-jar %s -h", PICARD))
    },
    warning = function(){
      system2("rm", "~/MitoRSoftware")
      createFile("~", "MitoRSoftware")
      setwd("~/MitoRSoftware")
      downloadMitoRSoftwares()
    },
    error = function(){
      system2("rm", "~/MitoRSoftware")
      createFile("~", "MitoRSoftware")
      setwd("~/MitoRSoftware")
      downloadMitoRSoftwares()
    },
    finally = {
      message("Download and Installation of softwares is done")
    }
    )
  }

}

BAMfromFastQ <- function(keepFastQ = TRUE){
  
  # Pide la carpeta donde se encuentran las lecturas R1 y R2
  dirPrinc <<- rstudioapi::selectDirectory()
  pacientes <- searchPatient(dirPrinc)
  pacienteR1 <- sprintf("%s/%s", dirPrinc, pacientes[[1]])
  pacienteR2 <- sprintf("%s/%s", dirPrinc, pacientes[[2]])
  titulo <<- pacientes[[3]]
  
  # Armamos las carpetas dentro del paciente: MitoR
  setwd(dirPrinc)
  system2("mkdir", sprintf('%s/MitoR', dirPrinc), wait = TRUE)
  setwd(sprintf('%s/MitoR', dirPrinc))
  
  referencia <- referenciaHG38
  
  # Buscar la cantidad de threads posibles para la maquina
  cant_thread  <- (as.integer(system2("nproc", "--all", stdout = TRUE))-1)
  
  cant_thread<-3 ##
  # Header que va a tener el .BAM - Puede ser cambiamble
  # Header-Contains information about the entire file, such as sample name,
  # sample length, and alignment method. Alignments in the alignments section
  # are associated with specific information in the header section.
  header <- sprintf('@RG\\tID:%s\\tLB:MitoR\\tSM:bar', titulo)
  
  
  # BWA mem - Mapea la alineacion del R1 y R2 con la referencia
  system2(BWA, sprintf("mem -t %s -o %s/MitoR/%s_R1R2.sam %s %s %s -R '%s'", cant_thread, dirPrinc, titulo, referencia, pacienteR1, pacienteR2, header), stdout = TRUE, wait = TRUE)
  mem -t 3 -o /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/MitoR/37339_R1R2.sam /home/daniela/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/37339_R1.fastq.gz /home/daniela/CNVAnalisisDefinitivo/PacientesMito/37339Mito/37339_R2.fastq.gz -R @RG\\tID:37339\\tLB:MitoR\\tSM:bar 
  
  
  # STEP 3 SAM to BAM
  system2(SAMTOOLS, sprintf("view -bS %s/MitoR/%s_R1R2.sam > %s/MitoR/%s_R1R2.bam", dirPrinc, titulo, dirPrinc, titulo), stdout = TRUE, wait = TRUE)
  
  # STEP 4 SortSam (En nuestro caso seria algo asi como un SortBAM)
  system2("java", sprintf("-jar %s SortSam -I %s/MitoR/%s_R1R2.bam -O %s/MitoR/%s_sortedR1R2.bam --VALIDATION_STRINGENCY SILENT --SORT_ORDER coordinate", PICARD, dirPrinc, titulo, dirPrinc, titulo), stdout = TRUE, wait = TRUE)
  
  #STEP 5 index sorted Bam (created .bam.bai)
  indexBam(sprintf("%s/MitoR/%s_sortedR1R2.bam", dirPrinc, titulo))
  
  if(keepFastQ==FALSE){
    system2("rm", "%s", paciente1)
    system2("rm", "%s", paciente2)
  }
  
  return(sprintf("%s/MitoR/%s_sortedR1R2.bam", dirPrinc, titulo))
}

