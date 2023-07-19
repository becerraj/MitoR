#' @title Generation of the plot for CNVs.
#' @description Generates a graph that represents the CNVs detected. The user will be able to
#' visualize the distribution, sizes, and positions of the variations within the whole mtDNA
#' and the genes that are affected by them.
#' @param CNV_calls dataframe that contains the information about the CNVs detected.
#' @return graph
#' @export
#' @import ggplot2
#' @import tidyverse
#' @import ggtext
#' @import dplyr
#' @import showtext

plot_cnv <- function(CNV_calls) {

  data("bedfileMito")
  cnv <- subset(CNV_calls, select = c("start", "end", "name.gene", "type"))
  colnames(cnv)[3] <- "name"
  cnv$col <- ifelse(cnv$type == "duplication", "green", "red")
  bed$col <- "black"

  todo <- rbind(bed[, 2:5], cnv[, c(1, 2, 3, 5)])
  todo <- todo[order(todo$start), ]
  todo$yname <- ifelse(todo$col == "black", ifelse(seq_along(todo$col) %% 2 == 0, 0.25, -0.25), ifelse(todo$col == "green", 0.75, -0.75))
  todo$y <- ifelse(todo$col == "black", 0, ifelse(todo$col == "green", 1, -1))


  #Ordenado segun start:
  gene_order <- unique(todo[order(todo$start), ]$name)
  todo$name <- factor(todo$name, levels = gene_order)
  #todo$col <- factor(todo$col)

  #Delete the rows which names are duplicated but keeping the one that is coloured:
  todo <- ungroup(filter(group_by(todo, name), !(n() > 1 & col == "black")))

  showtext_auto()
  graph <- ggplot(todo, aes(x = start, y = name)) +
    geom_segment(aes(xend = end, yend = name), color = todo$col,  size = 1.2) +
    scale_color_manual(values = c("red", "green", "black")) +
    scale_fill_manual(values = c("red", "green", "black")) +
    labs(x = "Posición", y = "Genes", color = "Color", fill = "Color") +
    theme_classic() +
    theme(axis.text.y = element_markdown(color = todo$col),
          axis.text.y.right = element_blank(),
          axis.text.y.left = element_text(vjust = 1, size= 14),
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 25),
          axis.text.x = element_text(size = 25) )
  graph
  #plot_cnv<- recordPlot()
  return(graph)
}


