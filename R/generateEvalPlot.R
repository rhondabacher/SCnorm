#' @title Internal plotting function.
#' @inheritParams plotCountDepth
#' @param MedExpr non-zero median expression for all genes.
#' @param SeqDepth sequencing depth for each cell/sample.
#' @param Slopes per gene estimates of the count-depth relationship.
#' @param Name name for plot title.
#' @param BeforeNorm whether dat have already been normalized.

#' @description Genes are divided into NumExpressionGroups = 10 equally sized
#'    groups based on their non-zero median expression. Slope densities are plot
#'    for each group. 
#' @return a plot of the un-normalized slope densities. 
#' @author Rhonda Bacher
#' 
#' @import grDevices

generateEvalPlot <- function(MedExpr, SeqDepth, Slopes, Name, 
                             NumExpressionGroups = 10, BeforeNorm = TRUE) {
    
    colors <- grDevices::colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), 
                bias=2)(n = NumExpressionGroups)
        

    ExprGroups <- splitGroups(MedExpr[intersect(names(MedExpr), names(Slopes))], 10)
    DensH <- getDens(ExprGroups, Slopes, "Height")

    YMax <- pmin(round(max(DensH), 2) + .2, 10) #just for plotting

    groupGenes <- unlist(lapply(seq_len(NumExpressionGroups), function(x) {
                   Y <- rep(x, length(ExprGroups[[x]]))
                   names(Y) <- names(ExprGroups[[x]])
                   return(Y)}))
    meltedGroups <- data.table::melt(groupGenes, value.name="Group")
    meltedGroups <- data.table::data.table(Group=meltedGroups$Group, Gene = rownames(meltedGroups))
    Slopes <- data.table::data.table(Slope = Slopes, Gene = names(Slopes))
    longData <- merge(meltedGroups, Slopes, by="Gene")
    
    MEDS <- round(vapply(ExprGroups, median, numeric(1)), 2)
    MEDS[1] <- paste0(MEDS[1], " (lowest)")
    MEDS[NumExpressionGroups] <- paste0(MEDS[NumExpressionGroups], " (highest)")
    
    
    QQ <- ggplot2::ggplot(longData, ggplot2::aes(x=longData$Slope, colour = factor(longData$Group))) +
          ggplot2::geom_density(size=1.2, ggplot2::aes(fill=factor(longData$Group)), alpha=0) + 
          ggplot2::scale_colour_manual( values = colors, guide = FALSE)+ 
          ggplot2::scale_fill_manual( values = colors, guide=FALSE, labels=MEDS)+ 
          ggplot2::theme_bw()+ ggtitle(Name)+
          ggplot2::labs(x="Slope", y="Density") + 
          ggplot2::guides(fill=ggplot2::guide_legend(override.aes = list(colour=colors, alpha=1),
                          title="Expression Group Medians"))+
          ggplot2::theme(
              axis.text=ggplot2::element_text(size=15,colour = "black"), 
      	       legend.key.size=ggplot2::unit(1.5,"lines"),
      	       text=ggplot2::element_text(size=15), 
              panel.border = ggplot2::element_rect(size=1, colour = "black"),
              legend.text=ggplot2::element_text(size=13,colour = "black"), 
              legend.title=ggplot2::element_text(size=15,colour = "black"))
        
    if(BeforeNorm == TRUE) {QQ <- QQ + ggplot2::geom_vline(xintercept=1, size=1, color="black")}
    if(BeforeNorm == FALSE) {QQ <- QQ + ggplot2::geom_vline(xintercept=0, size=1, color="black")}
    
  plot(QQ) 
  return(PLOTDATA = as.data.frame(longData))
}
    