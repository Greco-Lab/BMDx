#this function takes in input the filtered list of bmd results and compute the intersection of the genes at different timepoints
venn_diagram_bmd_genes_across_time_point = function(bmd_list){
  GL = list()
  for(i in 1:length(bmd_list)){
    BMD_tab <- bmd_list[[i]][[1]]
    if(!is.null(BMD_tab) & nrow(BMD_tab)>0){
      GL[[names(bmd_list)[i]]] =  BMD_tab[,"Gene"]
    }
  }
  
  ItemsList = venn(GL)
  ItemsList = attr(ItemsList,"intersections")
  
  return(list(GL = GL, ItemsList=ItemsList))
}

# given the list of genes with bmd effect in the venn diagram, it build a dataframe with a list of genes in every column
# corresponding to the genes in the different intersections part
build_dataframe_from_genes_in_venn_diagrams = function(ItemsList){
  df = matrix("", nrow = max(unlist(lapply(ItemsList, length))), ncol = length(ItemsList))
  colnames(df) = names(ItemsList)
  for(i in 1:length(ItemsList)){
    if(length(ItemsList[[i]])>0){
      df[1:length(ItemsList[[i]]),i] = ItemsList[[i]]
    }
  }
  
  return(df)
}

# this function cluster the genes that show bmd at all time point by computing their correlations across their bmd values at the different time points
create_gene_bmd_dataframe_and_cluster_genes_by_bmd = function(bmd_list,ItemsList,hmethod, nclust, intersectionName){
  #innerset = unlist(strsplit(x = names(ItemsList)[[length(ItemsList)]],split = ":"))
  innerset = unlist(strsplit(x = intersectionName,split = ":"))
  
  BMD_gene = c() # it will contain a matrix with a number of rows equal to the number of genes and number of column equal to the timepoints.
                 # the matrix will contain the bmd value of every gene at different timepoint and it will be used to cluster the genes accoring 
                 # to their bmd patterns.
  XX = c() # it will contain a dataframe with the following column: Gene, BMD, TimePoint, Clustering
  for(i in innerset){
    x = bmd_list[[i]][[1]]
    rownames(x) = as.character(x[,1])
    XX = rbind(XX,cbind(x[ItemsList[[intersectionName]],c("Gene","BMD")],i))
    BMD_gene = cbind(BMD_gene,x[ItemsList[[intersectionName]],"BMD"])
  }
  
  colnames(XX) = c("Gene","BMD","TimePoint")
  XX = as.data.frame(XX)  
  
  colnames(BMD_gene) = innerset
  rownames(BMD_gene) = as.character(x[ItemsList[[intersectionName]],"Gene"])
  
  if(nrow(BMD_gene)>1 && ncol(BMD_gene)>1){ #if there are at least two genes and two timepoints I can cluster genes based on their BMD across time
    # hierarchical clustering of the genes
    hls = hclust(as.dist(1-cor(t(BMD_gene))), method = hmethod)
    #plot(hls)
    
    cls = cutree(hls, nclust)
    XX = cbind(XX, cls[as.character(as.vector(XX[,1]))])
    colnames(XX)[4] = "Cluster"
  }else{ #if there is one gene and multiple time point, there is only 1 cluster with the single gene
    if(ncol(BMD_gene)>1){
      XX = cbind(XX, 1)
      colnames(XX)[4] = "Cluster"
      hls = NULL
    }else{ #if there is only one time point, cannot cluster genes, I will make a scatterplot instead
      XX = cbind(XX, 0)
      colnames(XX)[4] = "Cluster"
      hls = NULL
    }

  }
  
  
  
  return(list(XX=XX, hls = hls, BMD_gene = BMD_gene))
}

