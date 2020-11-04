#plots a matrix having colums for wich a hierachy is provided
# path_mat is matrix to be plotted
# path_hier is the hierarchy defined over columns
# experiment_ann is a vector of the same length as # of rows of path_mat defining grouping for samples
# discrete tells the function if the value are to be plotted using a countinuos scale or a discrete scale
# square_colors if a discrete scale is chosen, than the colours for each possible value have to be provided
# color_leg     if a discrete scale is chosen, than the colour legend for each possible value has to be provided
# level_col   level (column number) from the hierarchy used to group columns (pathways)

plot_grid <- function(path_mat,path_hier, title="", experiment_ann=c(),discrete=FALSE,square_colors=c(),color_leg=c(),level_col=1,treat_text_size=8,path_text_size=6, asRatio = TRUE) {
  #save(path_mat, path_hier,experiment_ann,title,discrete,square_colors,color_leg,
  #     level_col,treat_text_size,path_text_size,file="demo/demo_plot.RData")
  
  #path_mat = path_mat[rownames(experiment_ann),]
  #define the groups from the hierarchy and the chosen level
  path_col <- factor(path_hier[,level_col], levels = unique(path_hier[,level_col]))
  #prepare a set of colors to assign to each group
  #darkcols <- c(brewer.pal(8, "Dark2"),brewer.pal(8, "Accent")[-4],brewer.pal(12, "Paired")[-11])
  
  #using ggplot2 we need to melt the input matrix
  kegg_melt <- reshape::melt(path_mat)
  colnames(kegg_melt)[1:2] = c("Var1","Var2")
  #if discrete create value-color mapping using input values or this default scale
  if (discrete){
    if(length(square_colors)==0){
      colors <- c("-1"="darkgreen","0"="white","1"="red")
      color_leg <- c("-1"= "negative","0"="neutral","1"="positive")
    }else {colors <- square_colors}
    #if discrete we need to convert values in the matrix in factors to esure ggplot will plot on a discrete scale
    kegg_melt$value <- as.factor(kegg_melt$value)
  }
  # if a grouping for the sample is provided assign it to the melted dataframe rows (ncol by ncol)
  #this will be used in faceting
  if(length(experiment_ann[,1])>0){
    #print(path_mat)
    #print(kegg_melt)
    #print(experiment_ann)
    #print(ncol(path_mat))
    x = experiment_ann[which(experiment_ann[,2] %in% rownames(path_mat)),1]
    #print(x)
    kegg_melt$experiment <- rep(x,ncol(path_mat))
  }else{
    #otherwise assign a default group to all samples
    kegg_melt$experiment <- as.factor("treatment")
  }
  
  #assign group fromf rom the hierarchy and the chosen level (path_col) to the melted dataframe rows (nrow by nrow)
  #this will be used in faceting
  kegg_melt$path_group <- factor(as.character(rep(path_col,each=nrow(path_mat))),levels=unique(as.character(rep(path_col,each=nrow(path_mat)))))
  
  #prepare a set of colors to assign to each group
  darkcols <- randomcoloR::randomColor(length(unique(kegg_melt$path_group)),luminosity = "dark")
  
  
  # prepare the ggplot instance -- to render groups we will use faceting
  #using tile plot to render thje matrix as a grid of squares
  
  
  
  #print(kegg_melt$experiment)
  #print(head(kegg_melt))
  ukm = unique(kegg_melt$experiment)
  #print(ukm)
  
  for(i in 1:nrow(kegg_melt)){
    kegg_melt[i,"experiment"] = which(ukm %in% kegg_melt[i,"experiment"] )
  }
  #print(head(kegg_melt))
  
  #Ordinamento per mantenere ordinamento i gruppi ordinati secondo ordine numerico
  labelx = sort(unique(as.numeric(as.character(kegg_melt$experiment))))
  #print(labelx)
  kegg_melt$experiment = factor(kegg_melt$experiment,levels = labelx)
  
  #Ordiniamo i campioni con lo stesso ordine delle foglie dell'albero gerarchico
  kegg_melt$Var1 = factor(kegg_melt$Var1,levels = unique(kegg_melt$Var1))
  
  if(discrete == FALSE){
    #the values in kegg_melt must be numeric, not factors
    kegg_melt$value = as.numeric(as.vector(kegg_melt$value))
  }
  
  shorten = function(string){return(string)}
  
  plot_obj <- ggplot(kegg_melt, aes(Var1, Var2, fill = value)) +
    facet_grid(path_group~experiment,scales="free",space="free",labeller = labeller(path_group=shorten)) +
    # theme(axis.text.x = element_text(colour = as.numeric(as.factor(experiment_ann))),
    #       axis.text.y = element_text(colour = darkcols[as.numeric(path_col)]))   +
    geom_tile(colour = "black",size=0.1) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    #coord_equal() +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0, size=treat_text_size),
          axis.text.y = element_text(size=path_text_size),
          #plot.margin = margin(1, 1, 1, 1, "cm"),
          strip.text.y = element_text(angle = 0,size = 9,face = "bold")) +
    scale_x_discrete(position = "top") +
    labs(x = "", y="")
  
  if(asRatio) {
    plot_obj = plot_obj+theme(aspect.ratio = 1)
  }
  
  if (discrete){
    plot_obj <- plot_obj + scale_fill_manual(values=colors,labels=color_leg,na.value = 'gray50')
  }else{
    #limits=c(min(kegg_melt$value[is.na(kegg_melt$value)==FALSE]), max(kegg_melt$value[is.na(kegg_melt$value)==FALSE]))
    plot_obj <- plot_obj + scale_fill_gradient(low = "darkgreen", high = "red",na.value = 'gray50')
  }
  
  #after faceting we need to open the ggplot object to put different row/column label colors to different facets
  gplot <- ggplotGrob( plot_obj )
  nms <- lapply( gplot$grobs , function(x) names( x[]$children ) )
  #we search for axis.line.x and axis.line.y containers
  grbs_x_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.x" , x ) ) , any ) == 1 )
  grbs_y_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.y" , x ) ) , any ) == 1 )
  
  #grob variable iterates over blocks of colums of the ggplot
  if(length(length(grbs_x_id))>0){
    for (grob in (1:length(grbs_x_id))){
      to_select <- grbs_x_id[grob]
      if (length(to_select) >= 1) { # Check added because there was a 0 length vector. Luca
        item_at_grob <- gplot$grobs[[to_select]]
        if (! is.null(item_at_grob)) { # Check added because there where nulls. Luca
          item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
        }
      }
    }
  }
  
  #grob variable iterates over blocks of rows of the ggplot
  if(length(length(grbs_y_id))>0){
    for (grob in (1:length(grbs_y_id))){
      to_select <- grbs_y_id[grob]
      if (length(to_select) >= 1) {
        item_at_grob <- gplot$grobs[[to_select]]
        if (! is.null(item_at_grob)) {
          item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=darkcols[grob]
        }
      }
    }
  }
  
  return(gplot)
}



# #plots a matrix having colums for wich a hierachy is provided
# # path_mat is matrix to be plotted
# # path_hier is the hierarchy defined over columns
# # experiment_ann is a vector of the same length as # of rows of path_mat defining grouping for samples
# # discrete tells the function if the value are to be plotted using a countinuos scale or a discrete scale
# # square_colors if a discrete scale is chosen, than the colours for each possible value have to be provided
# # color_leg     if a discrete scale is chosen, than the colour legend for each possible value has to be provided
# # level_col   level (column number) from the hierarchy used to group columns (pathways)
# 
# plot_grid <- function(path_mat,path_hier, title="", experiment_ann=c(),discrete=FALSE,square_colors=c(),color_leg=c(),level_col=1,treat_text_size=8,path_text_size=6, asRatio = TRUE) {
#   #save(path_mat, path_hier,experiment_ann,title,discrete,square_colors,color_leg,
#   #     level_col,treat_text_size,path_text_size,file="demo/demo_plot.RData")
#   
#   #path_mat = path_mat[rownames(experiment_ann),]
#   #define the groups from the hierarchy and the chosen level
#   path_col <- factor(path_hier[,level_col], levels = unique(path_hier[,level_col]))
#   #prepare a set of colors to assign to each group
#   #darkcols <- c(brewer.pal(8, "Dark2"),brewer.pal(8, "Accent")[-4],brewer.pal(12, "Paired")[-11])
#   
#   #using ggplot2 we need to melt the input matrix
#   kegg_melt <- reshape::melt(path_mat)
#   colnames(kegg_melt)[1:2] = c("Var1","Var2")
#   #if discrete create value-color mapping using input values or this default scale
#   if (discrete){
#     if(length(square_colors)==0){
#       colors <- c("-1"="darkgreen","0"="white","1"="red")
#       color_leg <- c("-1"= "negative","0"="neutral","1"="positive")
#     }else {colors <- square_colors}
#     #if discrete we need to convert values in the matrix in factors to esure ggplot will plot on a discrete scale
#     kegg_melt$value <- as.factor(kegg_melt$value)
#   }
#   # if a grouping for the sample is provided assign it to the melted dataframe rows (ncol by ncol)
#   #this will be used in faceting
#   if(length(experiment_ann[,1])>0){
#     #print(path_mat)
#     #print(kegg_melt)
#     #print(experiment_ann)
#     #print(ncol(path_mat))
#     x = experiment_ann[which(experiment_ann[,2] %in% rownames(path_mat)),1]
#     #print(x)
#     kegg_melt$experiment <- rep(x,ncol(path_mat))
#   }else{
#     #otherwise assign a default group to all samples
#     kegg_melt$experiment <- as.factor("treatment")
#   }
#   
#   #assign group fromf rom the hierarchy and the chosen level (path_col) to the melted dataframe rows (nrow by nrow)
#   #this will be used in faceting
#   kegg_melt$path_group <- factor(as.character(rep(path_col,each=nrow(path_mat))),levels=unique(as.character(rep(path_col,each=nrow(path_mat)))))
#   
#   #prepare a set of colors to assign to each group
#   darkcols <- randomcoloR::randomColor(length(unique(kegg_melt$path_group)),luminosity = "dark")
#   
#   
#   # prepare the ggplot instance -- to render groups we will use faceting
#   #using tile plot to render thje matrix as a grid of squares
#   
#   
#   
#   #print(kegg_melt$experiment)
#   #print(head(kegg_melt))
#   ukm = unique(kegg_melt$experiment)
#   #print(ukm)
#   
#   for(i in 1:nrow(kegg_melt)){
#     kegg_melt[i,"experiment"] = which(ukm %in% kegg_melt[i,"experiment"] )
#   }
#   #print(head(kegg_melt))
#   
#   #Ordinamento per mantenere ordinamento i gruppi ordinati secondo ordine numerico
#   labelx = sort(unique(as.numeric(as.character(kegg_melt$experiment))))
#   #print(labelx)
#   kegg_melt$experiment = factor(kegg_melt$experiment,levels = labelx)
#   
#   #Ordiniamo i campioni con lo stesso ordine delle foglie dell'albero gerarchico
#   kegg_melt$Var1 = factor(kegg_melt$Var1,levels = unique(kegg_melt$Var1))
#   
#   if(discrete == FALSE){
#     #the values in kegg_melt must be numeric, not factors
#     kegg_melt$value = as.numeric(as.vector(kegg_melt$value))
#   }
#   
#   shorten = function(string){return(string)}
#   
#   plot_obj <- ggplot(kegg_melt, aes(Var1, Var2, fill = value)) +
#     facet_grid(path_group~experiment,scales="free",space="free",labeller = labeller(path_group=shorten)) +
#     # theme(axis.text.x = element_text(colour = as.numeric(as.factor(experiment_ann))),
#     #       axis.text.y = element_text(colour = darkcols[as.numeric(path_col)]))   +
#     geom_tile(colour = "black",size=0.1) +
#     ggtitle(title) +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     #coord_equal() +
#     theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0, size=treat_text_size),
#           axis.text.y = element_text(size=path_text_size),
#           #plot.margin = margin(1, 1, 1, 1, "cm"),
#           strip.text.y = element_text(angle = 0,size = 9,face = "bold")) +
#     scale_x_discrete(position = "top") +
#     labs(x = "", y="")
#   
#   if(asRatio) {
#     plot_obj = plot_obj+theme(aspect.ratio = 1)
#   }
#   
#   if (discrete){
#     plot_obj <- plot_obj + scale_fill_manual(values=colors,labels=color_leg,na.value = 'gray50')
#   }else{
#     #limits=c(min(kegg_melt$value[is.na(kegg_melt$value)==FALSE]), max(kegg_melt$value[is.na(kegg_melt$value)==FALSE]))
#     plot_obj <- plot_obj + scale_fill_gradient(low = "darkgreen", high = "red",na.value = 'gray50')
#   }
#   
#   #after faceting we need to open the ggplot object to put different row/column label colors to different facets
#   gplot <- ggplotGrob( plot_obj )
#   nms <- lapply( gplot$grobs , function(x) names( x[]$children ) )
#   #we search for axis.line.x and axis.line.y containers
#   grbs_x_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.x" , x ) ) , any ) == 1 )
#   grbs_y_id <- which( sapply( lapply( nms , function(x) grepl( "axis.line.y" , x ) ) , any ) == 1 )
#   
#   #grob variable iterates over blocks of colums of the ggplot
#   for (grob in (1:length(grbs_x_id))){
#     to_select <- grbs_x_id[grob]
#     if (length(to_select) >= 1) { # Check added because there was a 0 length vector. Luca
#       item_at_grob <- gplot$grobs[[to_select]]
#       if (! is.null(item_at_grob)) { # Check added because there where nulls. Luca
#         item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=grob
#       }
#     }
#   }
#   #grob variable iterates over blocks of rows of the ggplot
#   for (grob in (1:length(grbs_y_id))){
#     to_select <- grbs_y_id[grob]
#     if (length(to_select) >= 1) {
#       item_at_grob <- gplot$grobs[[to_select]]
#       if (! is.null(item_at_grob)) {
#         item_at_grob$children$axis$grobs[[1]]$children[[1]]$gp$col=darkcols[grob]
#       }
#     }
#   }
#   
#   return(gplot)
# }
