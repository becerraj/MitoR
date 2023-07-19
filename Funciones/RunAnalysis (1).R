#Just the first time:
  #check if the packages are already installed
  #create databases empty
if(system.file(package='ExomeDepth') == ""){
  Instalations()
  message("Packages required have been succesfully installed")
  
}

#Load Libraries
Libraries()



#Check if softwares for processing FASTQ files are downloaded

if(!("MitoRSoftware" %in% existingFiles("~"))){
  d<-readline(prompt = "You will need the following softwares:
               BWA, GATK, Picard and Samtools.
               Do you want to start with installation? (Y/N)")

  if(d=="Y"|d=="y"){ InstallationsforFASTQ()
  }else{
      stop()
    }
}

BWA<-"/home/daniela/MitoRSoftware/BWA/usr/bin/bwa"
GATK<-'/home/daniela/MitoRSoftware/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar'
PICARD<- '/home/daniela/MitoRSoftware/picard-2.27.5/picard.jar'
SAMTOOLS<-"/home/daniela/MitoRSoftware/Samtools/samtools-1.16.1/samtools"
referenciaHG38 <- "home/daniela/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta"


#Ask which analysis you want to do:
TYPE<-readline(prompt = "Welcome to CNV Analisis.
             Write N if you want to analyze the nuclear genome.
             Write M if you want to analyze the mitochondrial genome: ")


#Load the databases and bedfiles depending which type was selected
#----------------------------------------------------------------------
#       WHOLE GENOME ANALYSIS
#----------------------------------------------------------------------
if(TYPE=="N"| TYPE== "n"){
  c<-"ControlDB"
  p<-"PatientsDB"
  e<-"ExonsDB"
  include.chr<-TRUE
  
  if(!(file.exists("~/DataBases/ControlDB.RDS"))){
    InitialDB(c,p,e)
    message("Data Bases have been created. You can find them in DataBases folder")
    
  }
  #ControlDB<-read_excel(paste("~/DataBases/",c,".xlsx", sep=""))
  ControlDB<-readRDS(paste("~/DataBases/",c,".RDS", sep=""))
  #ControlDB<-as.data.frame(ControlDB)
  
  #PatientsDB<-read_excel(paste("~/DataBases/",p,".xlsx", sep=""))
  PatientsDB<-readRDS(paste("~/DataBases/",p,".RDS", sep=""))
  #PatientsDB<-as.data.frame(PatientsDB)
  
  #ExonsDB<-read_excel(paste("~/DataBases/",e,".xlsx", sep=""))
  ExonsDB<-readRDS(paste("~/DataBases/",e,".RDS", sep=""))
  #ExonsDB<-as.data.frame(ExonsDB)
  
  bed<-readRDS("~/CNVAnalisisDefinitivo/Data-CNV/bedfileHuman.RDS")
  

#----------------------------------------------------------------------
#       MITOCHONDRIAL ANALYSIS
#----------------------------------------------------------------------
  }else if(TYPE=="M"| TYPE== "m"){ 
    c<-"ControlMitoDB"
    p<-"PatientsMitoDB"
    e<-"GenesMitoDB"
    include.chr<-FALSE
    
    if(!(file.exists("~/DataBases/ControlMitoDB.xlsx"))){
      InitialDB(c, p, e)
      message("Data Bases have been created. You can find them in DataBases folder")
    }
    
    #ControlDB<-read_excel(paste("~/DataBases/",c,".xlsx", sep=""))
    ControlDB<-readRDS(paste("~/DataBases/",c,".RDS", sep=""))
    
    #ControlDB<-as.data.frame(ControlDB)
    #PatientsDB<-read_excel(paste("~/DataBases/",p,".xlsx", sep=""))
    PatientsDB<-readRDS(paste("~/DataBases/",p,".RDS", sep=""))
    #PatientsDB<-as.data.frame(PatientsDB)
    
    #ExonsDB<-read_excel(paste("~/DataBases/",e,".xlsx", sep=""))
    ExonsDB<-readRDS(paste("~/DataBases/",e,".RDS", sep=""))
    #ExonsDB<-as.data.frame(ExonsDB)
    
    bed<-readRDS("~/CNVAnalisisDefinitivo/Data-CNV/bedfileMito.RDS")
    
}



# OFFER 3 FUNCTIONS:
  #- Add a normal patient to control db
  #- Add a patient to db
  #- Analyze a patient: ACORDARSE que las DB tienen que ser dataframes!

a<-readline(prompt = "Let's start!
            Press AC if you want to add a normal patient for control DB.
            Press RC if you want to remove a normal patient from control DB.
            Press AP if you want to add a patient to the data base.
            Press RP if you want to remove a patient from the data base.
            Press A if you want to analyze a patient: ")

id<-readline(prompt = "Patient ID? ")


if(a=="AC"|a=="ac"){
  ControlDB<-AddControl(ControlDB, PatientsDB, id, bed, c, include.chr)

}else if (a=="AP"|a=="ap"){
  PatientsDB<-AddPatient(ControlDB, PatientsDB, id, bed, p, include.chr)

}else if (a=="RC"|a=="rc"){
  id<-"334"
  ControlDB<-RemovePatient(ControlDB, id)
  saveRDS(ControlDB, paste("~/DataBases/", c, ".RDS", sep=""))
  
}else if (a=="RP"|a=="rp"){
  PatientsDB<-RemovePatient(PatientsDB, id)
  saveRDS(PatientsDB, paste("~/DataBases/", p, ".RDS", sep=""))
  
}else if(a=="A"|a=="a"){
  if(max(dim(ControlDB))==0){
    stop("You need at least 1 control patient to do the analisis")
  }
  outAnalisis<-AnalyzePatient(ControlDB, PatientsDB, id, bed, p, include.chr)
  PatientsDB<-outAnalisis[[1]]
  CNV_calls<-outAnalisis[[2]]
  all.exons<-outAnalisis[[3]]

  #Update completo
  ExonsDB<-UpdateExonsDB(CNV_calls, id, e)
  message("Exons DB has been succesfully updated")

  #OrderCNV
  if(TYPE=="N"|TYPE=="n"){
    CNV_calls<-OrderCNV(CNV_calls)  
  }else if(TYPE=="M"|TYPE=="m"){
    CNV_calls <- subset(CNV_calls, select=c("chromosome", "type", "size", "id", "genes", "nGenes", "nexons","BF", "reads.expected", "reads.observed", "reads.ratio"))
    colnames(CNV_calls)[1]<-"chr"
  }
  

  #Visualization
  geneTable <- openxlsx::read.xlsx( "~/CNVAnalisis-Package/CNVAnalisis/Data-CNVs/GenesTargetCNVs.xlsx")
  load("~/CNVAnalisis-Package/CNVAnalisis/Data-CNVs/Centromere.RData") #indicates the location of each arm of the chromosome
  geneTable$chromosome <- stringr::str_remove_all(geneTable$chromosome,"[a-z]")

  out_data<-DataforVisualization(all.exons)
  Todos_puntos<-out_data[[1]]
  Todos_lp<-out_data[[2]]
  Todos_lq<-out_data[[3]]
  Todos_cnvs<-out_data[[4]]
  div_chr<-out_data[[5]]
  cent<-out_data[[6]]

  graph_<-GenerateGraph(Todos_puntos, Todos_lp, Todos_lq, div_chr, cent)
  graph_

  CNV_calls<-GenerateLinks(CNV_calls)


  Report(CNV_calls, graph=graph_, id)

  message("The analisis has finished! Your will find the report on Results folder
        and the 3 databases(Patients, Control and Exons) on DataBases folder")
}else{
  warning("Insert a valid letter")
}






#_________________________TESTS________________________________________________

#TEST 1: CNVfromBAM
#Comparacion de conteos entre conteos original y el obtenido de CNVfromBAM:

#comparacion<-cts[,5] == DB[,2]
#table(comparacion)["FALSE"] #0


#TEST 2: CNVfrom Counts

#Comparaciòn entre cnvs obtenidos desde bam y desde counts
#CNV_calls==CNV_calls_

#Comparaciòn de tamaños de regiones con variacion segun diferentes minoverlaps:
#out_counts_0.0001<-CNVfromCounts(id="30452", DB, minoverlap = 0.0001)
#CNV_calls_counts_0.0001<-out_counts_0.0001[[1]]
#all.exons_counts<-out_counts_0.0001[[2]]
#Annot_counts<- out_counts_0.0001[[3]]

#different minoverlap
#out_counts_0.001<-AnalisisCNV(id="30452", DB, minoverlap= 0.001)
#CNV_calls_counts_0.001<-out_counts_0.001[[1]]
#all.exons<-out_counts_0.001[[2]]

#CNV_calls_counts_0.0001$size==CNV_calls_counts_0.001$size
  #Mismos maximos y minimos de los tamaños de las regiones con variaciòn



#save(CNV_calls, all.exons, file="out_30452.RData")
#load("out_30452.RData")
#_____________________________________________________________________________



