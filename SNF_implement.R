

source('./snf_aux_shiny.R')


check=function(features_list){
fl=lapply(features_list,read_shiny_files)
}

pprocess_files=function(features_list,rm,rmv,im,imv,nr){


 #remove rows and columns with more than C missing values
 if(rm==TRUE) features_list=lapply(features_list,remove_missing,rmv)


  #impute. Method based on Troyanskaya et al
 if(im==TRUE) features_list=lapply(features_list,function(x,imv) t(impute.knn(t(x),k=imv)$data),imv)#####Catch errors!!!!!!!!!!!!!!!!!!!!!!!!!

 if(nr==TRUE) features_list=lapply(features_list,standardNormalization)

 return(features_list)
}

makeWf=function(K,alpha,tau,features_list){
#return(as.matrix(features_list))


  ###patients common among the different data
  cp=as.matrix(common_patients(features_list))

  features_list=lapply(features_list,function(x,p) x[p,],cp)



  ##generate affinity matrices
  Ws = mapply(aff_matrices,features_list,rep("cont",length(features_list)),K=K,alpha=alpha,SIMPLIFY=F)

    #generate the overall fused matrix
  Wf <-SNF(Ws, K, T)

  #assign to rows and colums the patient names (3 first components of TCGA barcodes)
  rownames(Wf)=cp
  colnames(Wf)=cp
  return(Wf)
}

number_of_clusters=function(Wf){
#estimate number of clusters by using eigen/gaps and rotation cost
  clusters=as.matrix(estimateNumberOfClustersGivenGraph(Wf))
  clusters=cbind(c("eigen-gaps best estimate","eigen-gaps 2nd best estimate","rotation cost best estimate","rotation cost 2nd best estimate"),clusters)
  #rownames(clusters)=c("eigen-gaps best estimate","eigen-gaps 2nd best estimate","rotation cost best estimate","rotation cost 2nd best estimate")
  colnames(clusters)=c("Method","Estimate number of clusters given graph")
  return(clusters)
}

#Cluster groups by using spectral clustering
cluster_groups=function(Wf,c){
  group = spectralClustering(Wf,c);
  names(group)=rownames(Wf)
  #recalculate clusters to match Wang paper
  group=order_group(group)
}



#Display heatmap with clustering results
display_heatmap=function(Wf,group){

  cols=c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
#   ##Visualize clusters
  cg=getColorsForGroups(group,colors=cols)#get group labels
  cl=colorRampPalette(c("black","cyan"))(1000)
#   #displayClustersWithHeatmap(Wf,group,cg,col=cl,labCol=NA,labRow="",cex.lab=3)#display heatmap
#
#   ###post-process files following BOs instructions
  Wc=post_process(Wf)
  displayClustersWithHeatmap(Wc,group,cg,col=cl,labCol=NA,labRow="",cex.lab=3)#display heatmap

}

#plo
nw_to_plot=function(nw,groups,top_percentage){
  colnames(nw)=c( "source" ,"target" ,"value" )
  top_cols=round(nrow(nw)*top_percentage)
  nw$source=groups$ID[match(nw$source,groups$name)]-1
  nw$target=groups$ID[match(nw$target,groups$name)]-1
  nw=nw[rev(order(nw$value)),]
  nw=nw[1:top_cols,]

}

 