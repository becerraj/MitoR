#install.packages("devtools")
#install.packages("rlang")

#Restart R sesion

# Guardar la nueva ruta en una variable
nueva_ruta <- "/home/efernandez/R/x86_64-pc-linux-gnu-library/4.2"

# Establecer la nueva ruta como la primera en la lista de rutas de la biblioteca
.libPaths(c(nueva_ruta, .libPaths()))

# Verificar que la nueva ruta haya sido establecida correctamente
.libPaths()

# Configurar la nueva ruta en la variable de entorno R_LIBS_USER
Sys.setenv(R_LIBS_USER = nueva_ruta)
Sys.getenv('R_LIBS_USER')

library(devtools)

setwd("~/RStudio/MitoR")
devtools::build()
devtools::document()
devtools::install()
library(MitoR)

#out_cnv<- analyze_cnv()

libPath <- "/home/efernandez/R/x86_64-pc-linux-gnu-library/4.2"

#mitor_sof <- sprintf("%s/mitorDB/Softwares", Sys.getenv('R_LIBS_USER'))
#mitor_sof <-"/home/efernandez/R/x86_64-pc-linux-gnu-library/4.2/mitorDB/Softwares"

path_dir <- file.path("/mnt/data/Mitocondrial/Pacientes", basename(select))
add_control_cnv(path_dir)

library(remotes)
install.packages("ExomeDepth")
install_github("DaniOrschanski/MitoR")


#ERRORRES

#37019: no esta en patientsdb
#algo del index
#NO tiene los fasta

#37739: esta en patientsdb
#SNPs&INDELs' analysis has been done!
#Error in file(description = xlsxFile) : invalid 'description' argument.
#Eso se arreglo porque cuando leia el wb era un dataframe. El error estaba al final de mitor_analysis
#Funciona OK

#40702; No tiene fasta

#46608: no estaba en patients db. Funciona OK. Lo agregue a control DB

#44562: funciona OK.
