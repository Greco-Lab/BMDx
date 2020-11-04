#' TODO input validation
#' plotting function
#' function calls
#' @import dplyr
#' @param path_mat
#' @param title, deafulat ""
#' @param experiment_ann, default c()
#' @param  gene_group, default NULL
#' @param discrete, default F
#' @param square_colors, default c()
#' @param color_leg, default c()
#' @param treat_text_size, default 8
#' @param asRatio, default TRUE
#' @return ggplot object, I think tuple/ list of plots/ combined plots
#' @export
#'


#Plot grid map for genes of a selected pathway
plot_grid_genes <- function(path_mat, title="", experiment_ann=c(), gene_group=NULL, discrete=F, square_colors=c(), color_leg=c(), treat_text_size=8, asRatio=TRUE){
  
  kegg_melt <- melt(path_mat)
  colnames(kegg_melt)[1:2] = c("Var1","Var2")
  
  #if discrete create value-color mapping using input values or this default scale
  if (discrete){
    if(length(square_colors)==0){
      colors <- c("-1"="darkgreen","0"="white","1"="red")
      color_leg <- c("-1"= "negative","0"="neutral","1"="positive")
    }else{
      colors <- square_colors
    }
    #if discrete we need to convert values in the matrix in factors to esure ggplot will plot on a discrete scale
    kegg_melt$value <- as.factor(kegg_melt$value)
  }
  
  # if a grouping for the sample is provided assign it to the melted dataframe rows (ncol by ncol)
  #this will be used in faceting
  
  # print("experiment_ann")
  #print(experiment_ann)
  if(length(experiment_ann[,1])>0){
    
    
    x = experiment_ann[which(experiment_ann[,2] %in% rownames(path_mat)),1]
    
    
    kegg_melt$experiment <- rep(x,ncol(path_mat))
  }else{
    #print("Creating fake 'experiment'!")
    #otherwise assign a default group to all samples
    kegg_melt$experiment <- as.factor("treatment")
  }
  
  #using tile plot to render thje matrix as a grid of squares
  
  
  ukm = unique(kegg_melt$experiment)
  
  
  for(i in 1:nrow(kegg_melt)){
    kegg_melt[i,"experiment"] = which(ukm %in% kegg_melt[i,"experiment"])
  }
  
  
  #Ordinamento per mantenere ordinamento i gruppi ordinati secondo ordine numerico
  
  labelx = sort(unique(as.numeric(as.character(kegg_melt$experiment))))
  
  
  kegg_melt$experiment = factor(kegg_melt$experiment, levels=labelx)
  
  #Ordiniamo i campioni con lo stesso ordine delle foglie dell'albero gerarchico
  
  kegg_melt$Var1 = factor(kegg_melt$Var1, levels=unique(kegg_melt$Var1))
  
  if(discrete == FALSE){
    #the values in kegg_melt must be numeric, not factors
    kegg_melt$value = as.numeric(as.vector(kegg_melt$value))
  }
  
  if(!is.null(gene_group)){
    
    kegg_melt$gene_group = gene_group
  }
  
  #Make pathway melt
  kegg_melt_path <- dplyr::filter(kegg_melt, gene_group=="Pathway")
  
  #Pathway map
  pPath <- ggplot(kegg_melt_path, aes(Var1, Var2, fill=value)) +
    facet_grid(.~experiment, scales="free", space="free") +
    geom_tile(colour="black", size=0.1) +
    ggtitle(title) +
    theme(plot.title=element_text(hjust=0.5)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0, size=treat_text_size),
          axis.text.y = element_text(size=treat_text_size, color="black"),
          strip.text.y=element_text(angle=0, size=9, face="bold")) +
    scale_x_discrete(position="top") +
    labs(x="", y="")
  
  if(asRatio) pPath = pPath+theme(aspect.ratio = 1)
  
  if (discrete){
    pPath <- pPath + scale_fill_manual(values=colors,labels=color_leg,na.value = 'gray50')
  }else{
    pPath <- pPath + scale_fill_gradient(low = "darkgreen", high = "red",na.value = 'gray50')
  }
  
  #afrer faceting we need to open the ggplot object to put different row/column label colors to different facets
  
  #png(file.path(tempdir(),"juen.png"))
  gplotPath <- ggplotGrob( pPath )
  #dev.off()
  nms <- lapply( gplotPath$grobs , function(x) names( x[]$children ) )
  #we search for axis.line.x and axis.line.y containers
  grbs_x_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.x" , x ) ) , any ) == 1 )
  grbs_y_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.y" , x ) ) , any ) == 1 )
  #grob variable iterates over blocks of colums of the ggplot
  for (grob in (1:length(grbs_x_id))){
    to_select <- grbs_x_id[grob]
    if (length(to_select) >= 1) {
      item_at_grob <- gplotPath$grobs[[to_select]]
      if (! is.null(item_at_grob)) {
        item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
      }
    }
  }
  
  #Make genes melt
  kegg_melt_gene <- dplyr::filter(kegg_melt, gene_group=="Genes")
  
  #print("str(kegg_melt_gene)")
  #print(str(kegg_melt_gene))
  
  ##Genes map
  p <- ggplot(kegg_melt_gene, aes(Var1, Var2, fill=value)) +
    facet_grid(.~experiment, scales="free", space="free") +
    geom_tile(colour="black", size=0.1) +
    ggtitle(title) +
    theme(plot.title=element_text(hjust=0.5)) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=0, size=treat_text_size),
          axis.text.y = element_text(size=treat_text_size, color="black"),
          strip.text.y=element_text(angle=0, size=9, face="bold")) +
    scale_x_discrete(position="top") +
    labs(x="", y="")
  
  
  if(asRatio) p = p+theme(aspect.ratio = 1)
  
  if (discrete){
    p <- p + scale_fill_manual(values=colors,labels=color_leg,na.value = 'gray50')
  }else{
    p <- p + scale_fill_gradient(low = "darkgreen", high = "red",na.value = 'gray50')
  }
  
  #afrer faceting we need to open the ggplot object to put different row/column label colors to different facets
  png(file.path(tempdir(),"juen.png"))
  gplot <- ggplotGrob( p )
  dev.off()
  nms <- lapply( gplot$grobs , function(x) names( x[]$children ) )
  #we search for axis.line.x and axis.line.y containers
  grbs_x_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.x" , x ) ) , any ) == 1 )
  grbs_y_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.y" , x ) ) , any ) == 1 )
  #grob variable iterates over blocks of colums of the ggplot
  for (grob in (1:length(grbs_x_id))){
    to_select <- grbs_x_id[grob]
    if (length(to_select) >= 1) {
      item_at_grob <- gplot$grobs[[to_select]]
      if (! is.null(item_at_grob)) {
        item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
      }
    }
  }
  
  
  #Combine ggplot for pathway and genes
  #gplotComb <- arrangeGrob(gplotPath, gplot, nrow=2, ncol=1, widths=1)
  gplotComb <- rbind(gplotPath, gplot, size="first")
  #print("str(gplotComb)")
  #print(str(gplotComb))
  return(gplotComb)
}
