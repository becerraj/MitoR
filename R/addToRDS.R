#' @title Add to RDS
#' @description MitoR-.
#' @param patient_XLSX_file es un ..
#' @import readxl
#' @return frequency_of_interest frecuencias de 

add_to_RDS <- function(patient_XLSX_file){
  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))

  newPatient <- read_excel(patient_XLSX_file, sheet = 1) # ["Freq_MitoR"]
  
  numPaciente <- basename(strsplit(patient_XLSX_file, "_")[[1]][1]) # [Nombre de la lista]
  
  filter_params <- as.character(read_excel(patient_XLSX_file, sheet = 2)["Filter_Params"]) # ["Filter_Params"]
  
  date <- as.character(read_excel(patient_XLSX_file, sheet = 2)["Date"]) # ["Date"]
  mutations <- paste(newPatient[, "REF"], newPatient[, "POS"], newPatient[, "ALT"], sep = "/") # ["Mutations"]

  # Una vez que tenemos todos los datos ya guardados en variables, vamos a guardarlo.
  if (!"RDS_DB.rds" %in% list.files(sprintf("%s/MitoRSoftware", mitor_files))){
    # Agregamos la columna de frecuencia en MitoR con valor 1 por ser el primero y ordenamos las columnas
    newPatient <- cbind(newPatient, Freq_MitoR = "100 %", Count = 1)

    newPatient <- newPatient[c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Freq_MitoR", "Count", "clinVAR",	"dbSNP",	"MitoMap",	"Omim",	"Disease", "Franklin", "VarSome")]

    # Agregamos todos los valores de la lista a la base de datos RDS
    RDS_DB <- list(Freq_MitoR = newPatient, list(XLSX_Path = sprintf("%s", patient_XLSX_file), Filter_Params = filter_params, Date = format(Sys.Date(), "%d/%m/%Y"),
                                                 Mutations = mutations))
    names(RDS_DB)[length(RDS_DB)] <- numPaciente
    frequency_of_interest <- RDS_DB$Freq_MitoR$Freq_MitoR

    saveRDS(RDS_DB, "~/MitoRSoftware/RDS_DB.rds")
  } else {
    RDS_DB <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))
    # Agregamos el nuevo paciente en un nuevo lugar de la lista. Si ya estaba de antes, tenemos que sacar el anterior y agregar el nuevo
    if (!RDS_DB[sprintf("%s", numPaciente)] %in% RDS_DB){
      RDS_DB <- readRDS("~/MitoRSoftware/RDS_DB.rds")

      # Agregamos el nuevo paciente junto con su nombre de elemento dentro de la lista
      RDS_DB[[length(RDS_DB)+1]] <- list(XLSX_Path = sprintf("%s", patient_XLSX_file), Filter_Params = filter_params, Date = format(Sys.Date(), "%d/%m/%Y"), Mutations = mutations)
      names(RDS_DB)[length(RDS_DB)] <- numPaciente

      results <- addFreq_MitoR(newPatient)
      Freq_MitoR <- results[1]
      frequency_of_interest <- results[2]

      RDS_DB$Freq_MitoR <- Freq_MitoR
      saveRDS(RDS_DB, "~/MitoRSoftware/RDS_DB.rds")

    } else { # Si ya estaba previamente, vamos a eliminar el anterior primero y agregar el nuevo analisis
      Freq_MitoR <- deleteFreq_MitoR(numPaciente)
      RDS_DB[[sprintf("%s", numPaciente)]] <- NULL

      results <- addFreq_MitoR(newPatient)
      Freq_MitoR <- results[1]
      frequency_of_interest <- results[2]

      RDS_DB$Freq_MitoR <- Freq_MitoR
      saveRDS(RDS_DB, "~/MitoRSoftware/RDS_DB.rds")
    }
  }
  return(frequency_of_interest) #Frecuencias de MitoR
}

#' @title Add Freq
#' @description MitoR-.
#' @param newPatient es un ..
#' @return Freq_MitoR y frequency_of_interest

addFreq_MitoR <- function(newPatient){
  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%/MitoRSoftware/RDS_DB.rds", mitor_files))
  Freq_MitoR <- as.data.frame(RDS_DB[1])
  colnames(Freq_MitoR) <- c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Freq_MitoR", "Count", "clinVAR", "dbSNP", "MitoMap", "Omim", "Disease", "Franklin", "VarSome")
  newPatient <- cbind(newPatient, Count = 1, Freq_MitoR = "100 %")
  newPatient <- newPatient[c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Freq_MitoR", "Count", "clinVAR", "dbSNP", "MitoMap", "Omim", "Disease",  "Franklin", "VarSome")]
  newAmount <- length(RDS_DB) # No le restamos uno porque recien al volver a la funcion anterior se guarda el nuevo paciente en el RDS

  # newMut sera un dataframe con las mutaciones nuevas. Que no estaban previamente en la DB
  newMut <- data.frame()
  # frequency_of_interest sera el vector que usaremos para cargarlo luego en el XLSX del paciente como Freq_MitoR
  frequency_of_interest <- c()
  # Inicializamos algunos parametros para poder recorrer la DB anterior y el paciente que ingresa
  i <- 1 # Relacionado a las filas de la DB existente (Freq_MitoR)
  j <- 1 # Relacionado a las filas del nuevo paciente (newPatient)
  Freq_MitoR$POS <- as.integer(Freq_MitoR$POS)
  newPatient$POS <- as.integer(newPatient$POS)

  # Vamos a recorrer la DB y la tabla del nuevo paciente simultaneamente. Para simplificar,
  # llamaremos a las dos tablas como tabla A y tabla B. Chequeamos segun las POS (posicion) de
  # las mutaciones: Se comparan las dos primeras filas, donde la tabla A tiene una mutacion
  # en una POS mas alta que la tabla B, entonces la tabla B pasara a la siguiente fila,
  # comparara de nuevo, y seguira haciendo esto hasta que la POS de la tabla B sea mayor
  # a la POS de la tabla A. En este caso, se invierten los roles, y la tabla A sera la
  # que comparara con su siguiente mutacion. Esto se repite hasta que ambas tablas llegan a
  # su valor de POS maximo.
  while ((i <= nrow(Freq_MitoR)) && (j <= nrow(newPatient))){

    # Este es el caso en que la tabla de la DB ha aumentado su POS hasta haber pasado a la
    # tabla del nuevo paciente sin haber encontrado la mutacion del paciente nuevo en su DB.
    # Esto quiere decir que la mutacion del paciente no se encuentra previamente en la DB, y
    # tiene que ser agregada como una mutacion nueva, entonces se agrega a la tabla newMut.
    if (Freq_MitoR$POS[i] > newPatient$POS[j]){
      newPatient$Freq_MitoR[j] <- sprintf("%s%s", 100/(newAmount), "%")
      frequency_of_interest <- c(frequency_of_interest, newPatient$Freq_MitoR[j])

      newMut <- rbind(newPatient[j, ], newMut)
      j <- j + 1

      # En este caso la POS del DB es menor a la POS de la mutacion del paciente nuevo, como
      # simplemente no ha coincidido todavia con la mutacion del paciente nuevo, avanza en POS
    } else if (Freq_MitoR$POS[i] < newPatient$POS[j]){
      Freq_MitoR$Freq_MitoR[i] <- sprintf("%s%s", (Freq_MitoR$Count[i]*100)/(newAmount), "%")
      i <- i + 1

      # Cuando ambas POS coinciden, se cambian los valores de la DB y pasamos a la siguiente
      # fila de las dos tablas.
    } else if (Freq_MitoR$POS[i] == newPatient$POS[j]){ # Chequea solamente la posicion para no generar tanto dato
      if ((Freq_MitoR$ALT[i] == newPatient$ALT[j]) && (Freq_MitoR$REF[i] == newPatient$REF[j])){ # Chequea el resto
        # Cambiamos la frecuencia alelica por las dudas haya cambiado el HMTVAR
        Freq_MitoR$All_freq_h[i] <- newPatient$All_freq_h[j]

        # Aumenta en uno el Count y cambia la frecuencia de nuestra DB
        Freq_MitoR$Count[i] <- (Freq_MitoR$Count[i] + 1)
        Freq_MitoR$Freq_MitoR[i] <- sprintf("%s%s", (Freq_MitoR$Count[i]*100)/(newAmount), "%")
        frequency_of_interest <- c(frequency_of_interest, Freq_MitoR$Freq_MitoR[i])
        j <- j + 1
        i <- i + 1
      }
    }
  }
  # Llegado el caso que llega al final de las POS del nuevo paciente, quiere decir que no hay
  # mas para comparar. En ese caso acomodamos los valores de frecuencia de las filas que no
  # llegaron a entrar dentro de la comparacion.
  if (i < nrow(Freq_MitoR)){
    Freq_MitoR$Freq_MitoR[i:nrow(Freq_MitoR)] <- sprintf("%s%s", (Freq_MitoR$Count[i:nrow(Freq_MitoR)]*100)/newAmount, "%")
  }

  # Si todas las mutaciones se repetian, newMut va a ser nulo y es un gasto de procesamiento innecesario
  if (length(newMut) != 0){
    colnames(newMut) <- c("CHROM", "POS", "REF", "ALT", "Freq_HMTVAR", "Freq_MitoR", "Count", "clinVAR", "dbSNP", "MitoMap", "Omim", "Disease", "Franklin", "VarSome")
    Freq_MitoR <- rbind(Freq_MitoR, newMut)
  }

  Freq_MitoR <- Freq_MitoR[order(Freq_MitoR$POS), ]
  Freq_MitoR$Freq_HMTVAR[Freq_MitoR$Freq_HMTVAR == "-"] <- "0"
  return(list(Freq_MitoR, frequency_of_interest))
}

#' @title Delete freq
#' @description MitoR-.
#' @param numPaciente es un ..
#' @return Freq_MitoR 
deleteFreq_MitoR <- function(numPaciente){
  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  RDS_DB <- readRDS(sprintf("%/MitoRSoftware/RDS_DB.rds", mitor_files))
  Freq_MitoR <- as.data.frame(RDS_DB[1])
  # Dividimos la columna de las mutaciones en ALT POS REF
  toDelete <- unlist(RDS_DB[[numPaciente]]["Mutations"], use.names = FALSE) %>%
    strsplit("/")
  newAmount <- length(RDS_DB) - 2

  for (i in (1:length(toDelete))){
    Freq_MitoR[(Freq_MitoR$POS == toDelete[[i]][2]) & (Freq_MitoR$ALT == toDelete[[i]][3]), ]$Count <- (Freq_MitoR[(Freq_MitoR$POS == toDelete[[i]][2]) & (Freq_MitoR$ALT == toDelete[[i]][3]), ]$Count - 1)

    # Si estamos eliminando la unica mutacion que habia, eliminamos esa fila
    if (Freq_MitoR[(Freq_MitoR$POS == toDelete[[i]][2]) & (Freq_MitoR$ALT == toDelete[[i]][3]), ]$Count == 0){
      Freq_MitoR <- Freq_MitoR[!((Freq_MitoR$POS == toDelete[[i]][2]) & (Freq_MitoR$ALT == toDelete[[i]][3])), ]
    }
  }
  # Vuelvo a configurar las frecuencias de las mutaciones dentro de la DB
  Freq_MitoR$Freq_MitoR <- sprintf("%s%s", (Freq_MitoR$Count * 100)/newAmount, "%")

  return(Freq_MitoR)
}
