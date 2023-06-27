#' @title Pile up Plot
#' @description MitoR-.
#' @param numPaciente, plot_mutation = TRUE, range = 5
#' @import stringr

pileupPlot <- function(numPaciente, plot_mutation = TRUE, range = 5) {
  results <- generateTXT_gene(numPaciente)
  positions <- results[1]
  mutations_at_genes <- results[2]

  mitor_files <- sprintf("%s/mitor", Sys.getenv('R_LIBS_USER'))
  mutations <- readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))[[sprintf("%s", numPaciente)]][4] %>%
    unlist() %>%
    stringr::str_split("/")

  path_to_MitoR <- substr(positions, start = 0, stop = (nchar(positions)-nchar(basename(positions))))

  plots_to_xlsx <- list()
  for (i in list.files(positions)) {
    # De aca hay que hacer el awk{$} solamente con el valor de la posicion y el nchar de las lecturas que nos corresponden
    system2(SAMTOOLS, sprintf("mpileup -f %s --positions %s/%s %s%s_dedup_reads.bam | awk '{print $1,$2,$3,$4,$5}' > %s%s_region.bam", referencia, positions, i, path_to_MitoR, numPaciente, path_to_MitoR, numPaciente), stdout = TRUE)
    data_pileup <- read.table(sprintf("%s/%s_region.bam", path_to_MitoR, numPaciente), header=FALSE, sep=" ")

    colnames(data_pileup) <- c("CHROM", "POS", "REF", "Total Counts", "REF Counts")

    # Contamos la cantidad de lecturas de REF (la cantidad de "." marca la cantidad de REF)
    for (j in (1:nrow(data_pileup))) {
      data_pileup[j,5] <- (sum(charToRaw(data_pileup[j,5]) == charToRaw(".")) + sum(charToRaw(data_pileup[j,5]) == charToRaw(",")))
    }

    # Voy a cambiarlo aca asi ya tenemos una tabla mas simplificada, sin los puntos y comas
    data_plot <- data_pileup[, c("POS", "REF", "Total Counts", "REF Counts")]
    data_plot$`REF Counts` <- as.integer(data_plot$`REF Counts`)
    data_plot$`ALT Counts` <- as.integer(data_plot$`Total Counts`) - data_plot$`REF Counts`

    ### Averiguar porque el ultimo dato del DF da un valor super grande ###
    data_plot <- data_plot[1:(nrow(data_plot)-1), ]

    # Guardamos las mutaciones encontradas en el gen que analizamos
    mutations_in_gen <- mutations_at_genes[stringr::str_remove(i, ".txt") == mutations_at_genes$name, "mut"] %>%
      strsplit(" ")

    # Aca hay que guardarlo en una lista
    gene_plot <- plot_gene_design(data_plot, mutations_in_gen)

    # Guardo en la lista con el nombre del gen
    plots_to_xlsx[[stringr::str_remove(i, ".txt")]] <- gene_plot

    # PLOTMUT---------
    if (plot_mutation == TRUE) {
      mutations_to_plot <- mutations_at_genes[mutations_at_genes$name == stringr::str_remove(i, ".txt"), "mut"] %>%
        strsplit(" ") %>% unlist()

      for (mutation in mutations_to_plot) {
        params <- strsplit(mutation, ">") %>% unlist()
        position <- params[2]
        alternative_base <- params[3]

        inferior_row <- which(data_plot$POS == (as.integer(position) - range))
        superior_row <- which(data_plot$POS == (as.integer(position) + range))
        data_mutation_plot <- data_plot[inferior_row : superior_row,]

        mutation_plot <- plot_mutation_design(data_mutation_plot, alternative_base)

        # Guardo este plot dentro de la Nested del gen
        plots_to_xlsx[[stringr::str_remove(i, ".txt")]][position] <- mutation_plot
      }
    }
  }
  file.remove(sprintf("%s%s_dedup_reads.bam", path_to_MitoR, numPaciente))
  file.remove(sprintf("%s%s_region.bam", path_to_MitoR, numPaciente))
}

#' @title Generation of TXT
#' @description MitoR-.
#' @param numPaciente blablaba
#' @return path_for_positions, mutations_at_genes

generateTXT_gene <- function(numPaciente) {
  # Load mutations
  path <- as.character(readRDS(sprintf("%s/MitoRSoftware/RDS_DB.rds", mitor_files))[[sprintf("%s", numPaciente)]][1])
  mutations <- readRDS("~/MitoRSoftware/RDS_DB.rds")[[sprintf("%s", numPaciente)]][4]
  mutations <-  stringr::str_split(mutations[[1]],"/", simplify= TRUE) %>% as.data.frame(row.names = NULL)
  colnames(mutations) <- c("ALT", "POS", "REF")

  # Load BED file with: Gene name, start position, end position
  gene_start_end <- readRDS("/home/juan/MitoRSoftware/bedfileMito.RDS")

  # Dealing with insertions and deletions on the length of the gene
  gene_start_end <- indels_in_gene(mutations, gene_start_end)

  # Funcion que busca en que genes se encuentran las mutaciones
  mutations_at_genes <- genes_of_mutations(mutations, gene_start_end)

  # Create a directory to store the genes sequences in txt files
  path_for_positions <- substr(path,start=0,stop=(nchar(path)-nchar(basename(path)))) %>%
    paste("GenesSequence", sep = "")

  if (file.exists(path_for_positions)){
    unlink(path_for_positions, recursive = TRUE)
  }

  dir.create(path_for_positions)

  # Store each DF with the genes positions on a txt
  for (j in (1:nrow(mutations_at_genes))) {
    gene_positions <- data.frame(mutations_at_genes[1, "chromosome"], (mutations_at_genes[j, "start"] : mutations_at_genes[j, "end"]))
    pathTXT <- sprintf("%s/%s.txt", path_for_positions, mutations_at_genes[j, "name"])
    write.table(gene_positions, pathTXT, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }

  return(path_for_positions, mutations_at_genes)
}

# Dandole la fila del datos_plot["REF"], nos va a decir cual de las cuatro bases es la que mas se repite. En nuestro caso va a
# ser siempre el ALT, porque seleccionamos las filas donde hay una mutacion.

searchALT <- function(base_read_String) {
  base_read_String <- toupper(base_read_String)
  counting <- table(strsplit(base_read_String, ""))
  return(names(counting)[which.max(counting)])
}

# Depending on the position of the instertion or deletion of the mutation, we change the length of the gene.
indels_in_gene <- function(mutations, start_end_gene) {
  j <- 1
  for (i in (1:nrow(mutations))) {
    indel <- nchar(mutations[i, "ALT"]) - nchar(mutations[i, "REF"])
    if (indel == 0) {
      # Continues
    } else {
      # Insertion or Deletion
      while(!((mutations[i, "POS"] > start_end_gene[j, "start"]) && (mutations[i, "POS"] < start_end_gene[j, "end"]))){
        j <- j + 1
      }
      start_end_gene[j, "end"] <- start_end_gene[j, "end"] + indel
      if (j < nrow(start_end_gene)){
        start_end_gene[j+1, "start"] <- start_end_gene[j+1, "start"] + indel
      }
    }
  }
  return(start_end_gene)
}

genes_of_mutations <- function(mutations, gene_start_end) {
  position_of_mutations <- as.integer(mutations$POS)
  i <- 1 # Filas de las posiciones
  j <- 1 # Filas de genes
  genes_of_mutations <- data.frame()

  # Como comentario esta la opcion para agregar esas mutaciones al primer gen. El tema es
  # que ni siquiera se van a ver reflejados en el plot porque no son parte de ese gen
  while (position_of_mutations[i] < gene_start_end$start[1]) {
    #if (!(gene_start_end$name[j] %in% genes_of_mutations$name)) {
    #  genes_of_mutations <- rbind(genes_of_mutations, gene_start_end[j, ])
      # Lo agregamos. No pasamos a la siguiente fila porque capaz tenemos mas dentro de la misma fila
    #}
    i <- i + 1
  }

  while (j < nrow(gene_start_end)) {
    while ((position_of_mutations[i] > gene_start_end$start[j]) && (position_of_mutations[i] < gene_start_end$end[j])) {
      if (!(gene_start_end$name[j] %in% genes_of_mutations$name)) {
        genes_of_mutations <- rbind(genes_of_mutations, cbind(gene_start_end[j, ], mut = paste0(mutations[i,], sep = "", collapse = ">")))
        # Lo agregamos. No pasamos a la siguiente fila porque capaz tenemos mas dentro de la misma fila
      } else {
        genes_of_mutations[nrow(genes_of_mutations), "mut"] <- paste(genes_of_mutations[nrow(genes_of_mutations), "mut"], paste0(mutations[i,], sep = "", collapse = ">"))
      }
      i <- i + 1
    }
    j <- j + 1
  }

  # Si siguen quedando mutaciones sin haber sido encasilladas en ningun gen es porque no son parte de los genes
  # del cromosoma mitocondrial.

  return(genes_of_mutations)
}


plot_gene_design <- function(data_plot, mutations_in_gen) {

  Position <- data_plot$POS
  Counts <- data_plot$`Total Counts`
  y2 <- data_plot$`ALT Counts`
  par(xpd = TRUE)
  # Primera línea
  plot(Position, Counts, type = "l", col = "black", main = "ND1", ylim = c(0, max(Counts)))

  # Segunda línea
  lines(Position, y2, type = "o", col = "red") # Mismos valores de X

  # Leyenda
  legend("bottomleft", legend = c("Total Counts", "ALT Counts"), col = c("black", "red"), lty = 1, cex = 0.5)

  # Mutaciones
  for (i in (1:length(mutations_in_gen[[1]]))) {
    mtext(mutations_in_gen[[1]][i], side = 3, adj = 1,
          line = (length(mutations_in_gen[[1]]) - i)) # Ponemos la anotacion que queramos segun las coordenadas
  }


  plot_of_gene <- recordPlot()

  return(plot_of_gene)
}

plot_mutation_design <- function(data_mutation_plot, alternative_base) {
  # X axis: positions
  Position <- data_mutation_plot$POS

  # Y axis: Counts
  Counts <- data_mutation_plot$`REF Counts`
  y2 <- data_mutation_plot$`ALT Counts`

  # ALT sequence. Only the mutation has its letter
  ALT_letters <- c(rep(".", range), alternative_base, rep(".", range))

  # Title: mutation
  mutation <- sprintf("%s>%s>%s", data_mutation_plot$REF[(nrow(data_mutation_plot)/2)+0.5], data_mutation_plot$POS[(nrow(data_mutation_plot)/2)+0.5], alternative_base)

  # First line of plot: REF sequence
  plot(Position, Counts, type = "b", pch = data_mutation_plot$REF, col = "black", main = mutation, ylim = c(-500, max(Counts)*1.1))

  # Second line of plot: ALT sequence. Only the mutation has its letter
  lines(Position, y2, type = "b", pch = ALT_letters, col = "red")

  # Legend always on the right because the mutation will always be at the centre
  legend("right", legend = c("REF Counts", "ALT Counts"), col = c("black", "red"), lty = 1, cex = 0.5)

  plot_of_mutation <- recordPlot()

  return(plot_of_mutation)
}


