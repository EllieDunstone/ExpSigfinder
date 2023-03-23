#' Plot mutational catalog
#' @param muts_basis 96-channel mutational catalog
#' @param colnum number of columns in the plot
#' @param h height of plot
#' @param w width of plot
#' @param outputname name of output
#' @export
plotCountbasis <- function(muts_basis, colnum, h, w, outputname) {
  
  muts_basis_melt <- reshape2::melt(muts_basis, "MutationType")
  names(muts_basis_melt) <- c("MutationType", "sample", "count")
  muts_basis_melt$mutation <- substr(muts_basis_melt$MutationType, 3, 5)
  
  mutation_order <- muts_basis[,c("MutationType", "MutationType")]
  names(mutation_order) <- c("MutationType", "mutation")
  mutation_order$mutation <- substr(mutation_order$MutationType, 3, 5)
  mutation_order <- mutation_order[order(mutation_order$mutation),]
  mypalette <- c("sky blue", "black", "red", "grey", "seagreen2", "lightpink2")
  grDevices::pdf(file = outputname, onefile = TRUE, width = w, height = h)
  
  p <- ggplot2::ggplot(data = muts_basis_melt, ggplot2::aes(x = MutationType, y = count, fill = mutation)) + 
    ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.8) + ggplot2::xlab("Mutation Types") + ggplot2::ylab("Count")
  
  p <- p + ggplot2::scale_x_discrete(limits = as.character(mutation_order$MutationType)) + ggplot2::ggtitle(outputname)
  p <- p + ggplot2::scale_fill_manual(values = mypalette)
  p <- p + ggplot2::theme(#axis.text.x=element_blank(),
    axis.text.x = ggplot2::element_text(size = 5, angle = 90, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
    axis.title.x = ggplot2::element_text(size = 15),
    axis.title.y = ggplot2::element_text(size = 15),
    plot.title = ggplot2::element_text(size = 10),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_rect(fill = "white"),
    panel.border = ggplot2::element_rect(colour = "black", fill = NA))
  p <- p + ggplot2::facet_wrap(~sample, ncol = colnum, scales = "free")
  
  print(p)
  grDevices::dev.off()
}

#' Plot mutational signature (percentage)
#' @param muts_basis 96-channel mutational catalog
#' @param colnum number of columns in the plot
#' @param h height of plot
#' @param w width of plot
#' @param outputname name of output
#' @export
plotPercentagebasis <- function(muts_basis, colnum, h, w, outputname) {
  
  muts_basis[,-which(colnames(muts_basis) == "MutationType")] <- muts_basis[,-which(colnames(muts_basis) == "MutationType")] / colSums(muts_basis[,-which(colnames(muts_basis) == "MutationType")])[col(muts_basis[,-which(colnames(muts_basis) == "MutationType")])]
  muts_basis_melt <- reshape2::melt(muts_basis, "MutationType")
  names(muts_basis_melt) <- c("MutationType", "sample", "count")
  muts_basis_melt$mutation <- substr(muts_basis_melt$MutationType, 3, 5)
  
  mutation_order <- muts_basis[,c("MutationType", "MutationType")]
  names(mutation_order) <- c("MutationType", "mutation")
  mutation_order$mutation <- substr(mutation_order$MutationType, 3, 5)
  mutation_order <- mutation_order[order(mutation_order$mutation),]
  
  mypalette <- c("sky blue", "black", "red", "grey", "seagreen2", "lightpink2")
  
  grDevices::pdf(file = outputname, onefile = TRUE, width = w, height = h)
  
  p <- ggplot2::ggplot(data = muts_basis_melt, ggplot2::aes(x = MutationType, y = count, fill = mutation)) + ggplot2::geom_bar(stat = "identity", position = "dodge", width = 0.8) + ggplot2::xlab("Mutation Types") + ggplot2::ylab("Percentage")
  p <- p + ggplot2::scale_x_discrete(limits = as.character(mutation_order$MutationType)) + ggplot2::ggtitle(outputname) + ggplot2::scale_y_continuous(limits = c(0, 0.5),breaks = seq(0, 0.5, 0.1), labels = scales::percent)
  p <- p + ggplot2::scale_fill_manual(values = mypalette)
  p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 5, colour = "black"),
               axis.text.y=ggplot2::element_text(size=15, colour = "black"),
               axis.title.x = ggplot2::element_text(size=15),
               axis.title.y = ggplot2::element_text(size=15),
               plot.title = ggplot2::element_text(size=10),
               panel.grid.minor.x=ggplot2::element_blank(),
               panel.grid.major.x=ggplot2::element_blank(),
               panel.grid.major.y = ggplot2::element_blank(),
               panel.grid.minor.y = ggplot2::element_blank(),
               panel.background = ggplot2::element_rect(fill = "white"),
               panel.border = ggplot2::element_rect(colour = "black", fill=NA))
  p <- p+ggplot2::facet_wrap(s~sample,ncol=colnum,scales = "free")
  
  print(p)
  grDevices::dev.off()
}
