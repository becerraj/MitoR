install.packages("miniUI", dependencies = TRUE)
install.packages("pkgdown", dependencies = TRUE)
install.packages("systemfonts", type = "source", repos = "http://cran.us.r-project.org", INSTALL_opts = "--no-lock", dependencies = TRUE)

#sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-devfrd
install.packages("devtools", dependencies = TRUE)

library(devtools)
devtools::create(path = "/home/daniela/TestLibrary") 

setwd("/home/daniela/TestLibrary")
devtools::build()
devtools::load_all() #para asegurarte de que se está ejecutando correctamente
devtools::document()
devtools::install()

#Testeo de funciones en un mismo .R
library(TestLibrary)
SayHello()
SayBye()

remove.packages("TestLibrary")

install.packages("TestLibrary")
#Testeo de import entre funciones del mismo paquete
TestImport("/home/daniela/carpeta1/archivo1.txt")

#Testeo de import entre funciones de diferentes paquetes
detach("package:ggplot2", unload = TRUE)

graph <-TestImportExterno()
print(graph)


#Testeo cargar datos al paquete
  #1. Creo carpeta data
dir.create("data")

  #2. Guardo el archivo en la carpeta data dentro de mi paquete

  #3. Creo una funcion que lea el archivo

data(bedfileMito)



#-----------------------------------------
#TESTEO FUNCIONES
mitor_files <- sprintf("%s/CNVMito", Sys.getenv('R_LIBS_USER'))
mitor_files<-"~"
BWA <<- "/home/daniela/MitoRSoftware/BWA/usr/bin/bwa"
Picard <<- "/home/daniela/MitoRSoftware/picard-2.27.5/picard.jar"
GATK <<- "/home/daniela/MitoRSoftware/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar"
Samtools <<- "/home/daniela/MitoRSoftware/Samtools/samtools-1.16.1/samtools"
referencia <<- "~/MitoRSoftware/RefHG38/Homo_sapiens.GRCh38.dna.chromosome.MT.fasta"


help(package = "CNVMito")
help(InitializeDB)
help("AnalyzePatient")

data("bedfileMito")

EjPacientesMito <- read_excel("~/CNVMito/Funciones/EjPacientesMito.xlsx")
PatientsDB <- subset(EjPacientesMito, select=2:4)
ControlDB<-subset(PatientsDB, select=2)
GenesDB<- data.frame()
saveRDS(list(ControlDB, PatientsDB, GenesDB), file = "~/DataBases/DBs.RDS")

initial <- initialize_db()
ControlDB <- initial[[1]]
PatientsDB <- initial[[2]]
GenesDB <- initial[[3]]


ControlDB <- AddControl(id = "37339", path_file = "~/Paciente/37339/37339_R1.fastq.gz")

PatientsDB <- AddPatient(id = "37339", path_file = "~/Paciente/47286/BM22-47286_R1.fastq.gz") 

PatientsDB <- RemoveAnyPatient(id="2.0", n= "Patients")

out_counts <- CNVfromCounts(id = "2.0", bed, ControlDB, PatientsDB, transit=0.9)
CNV_calls <- out_counts[[1]]

out <- AnalyzePatient(path_file = "~/Paciente/47286/BM22-46608_R1.fastq.gz")
PatientsDB <- out[[1]]
CNV_calls <- out[[2]]
GenesDB <- out[[3]]
CNVs<-out[[4]]
graph <- out[[5]]
print(graph)

IndicateFunction(a="AC", path_file = "~/Paciente/37019/37019_R1.fastq" )
IndicateFunction(a="CNV", path_file = "~/Paciente/47286/BM22-47286_R1.fastq" )


#No funciona:
GenerateMitoGraph(CNV_calls, bed)




unzip("~/BWA.zip", exdir="~/BWA")

#---------------------------------------------------------------------------------
library(devtools)
devtools::create(path = "/home/daniela/MitoR") 

setwd("/home/daniela/MitoR")
devtools::build()

devtools::load_all() #para asegurarte de que se está ejecutando correctamente
devtools::document()
#devtools::install_local() #Para que no se haga el lock
devtools::install()

#Para el final para guardarlo en github:
devtools::install_github("DaniOrschanski/MitoR")

library(MitoR)

initial <- initialize_db()
ControlDB <- initial[[1]]
PatientsDB <- initial[[2]]
GenesDB <- initial[[3]]

ControlDB <- add_control_cnv("~/Paciente/47286")
PatientsDB <- add_patient_cnv("~/Paciente/37339")

out <- analyze_cnv(path_dir = "~/Paciente/37339")
PatientsDB <- out[[1]]
CNV_calls <- out[[2]]
GenesDB <- out[[3]]
CNVs<-out[[4]]
graph <- out[[5]]
print(graph)

id<-"47286"
bam.control<- "~/Paciente/47286/MitoR/38069_sortedR1R2.bam"
cts <- getBamCounts(bed.frame = bed, 
                    bam.files = bam.control ,
                    include.chr = F,
                    referenceFasta = NULL)

Pcounts <- data.frame(Pcounts=cts[, ncol(cts)])
PatientsDB <- add_column(PatientsDB, Pcounts)
names(PatientsDB)[names(PatientsDB) == "Pcounts"] <- id

ControlDB<-PatientsDB
#Save the updated Patients DB


v_DB<-readRDS("~/DataBases/DBs.RDS")
v_ControlDB <- v_DB[[1]]
v_PatientsDB <- v_DB[[2]]
v_GenesDB <- v_DB[[3]]

ControlDB<- v_ControlDB
PatientsDB<-v_PatientsDB
GenesDB<- v_GenesDB

saveRDS(list(ControlDB, PatientsDB, GenesDB), file = paste(Sys.getenv('R_LIBS_USER'), "/mitorDB/cnvDB.RDS", sep=""))

#Subir a github
install.packages("usethis")
library(usethis)
create_github_token()

install.packages("gitcreds")
library(gitcreds)
gitcreds_set()
#Token: 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

use_git()
use_github()


#---------------------------------------------------------------------------------
library(devtools)
devtools::create(path = "/home/daniela/mitor") 

setwd("/home/daniela/mitor")
devtools::build()

devtools::load_all() #para asegurarte de que se está ejecutando correctamente
devtools::document()
devtools::install_local() #Para que no se haga el lock

library(mitor)
