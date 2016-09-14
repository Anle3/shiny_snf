

#### order groups based on their size
order_group=function(group,descending=T,letter=F){
  if (descending==T) ord=rev(order( table(group))) else ord=order( table(group))

  tmp=group

  if(letter==F){
    k=1
    for(i in ord){
      tmp[group==i]=k
      k=k+1
    }
  }else{
    k=letters
    for(i in ord){

      tmp[group==i]=k[i]}
  }
  group=tmp
  rm(tmp)
  return(group)
}


remove_missing=function(table,cutoff=0.2){
###this function removes rows and columns from a data.frame with more than x% missing values where x% is defined by cutoff
###columns###
  col_cutoff=round(dim(table)[1]*cutoff)
  include_columns=which(!colSums(is.na(table))>col_cutoff)
  table=table[,include_columns]

###rows###
  row_cutoff=round(dim(table)[2]*cutoff)
  include_rows=which(!rowSums(is.na(table))>row_cutoff)
  table=table[include_rows,]
  return(table)
}

#fread seems problematic, have to give ful path to work






aff_matrices = function(feature_matrix,type=c("cont","disc","bin"),K,alpha,norm=T){

  ###Arguments##
  ##The feature_matrix:the feature table
  ##type:feature data type e.g continuous discrete
  ##K,alpha: parameters for the generation of affinity matrix
  ##norm: normalize data before generating affinity matrix
  ####

  ##The function returns a weight matrix (affinity matrix)



  type = match.arg(type)

  if(type == "cont"){

    Dist = dist2(as.matrix(feature_matrix),as.matrix(feature_matrix));
  }else if(type == "disc"){
    Dist = chiDist2(as.matrix(feature_matrix),as.matrix(feature_matrix));
  }

  W = affinityMatrix(Dist, K, alpha)

  return(W)

}

###get patient names from TCGA barcodes

patients=function(x){
  x=(do.call(rbind, strsplit(x,"\\-|\\."))[,1:3])
  x=apply(x,1,function(x) paste(x[1:3],collapse="-"))
  return(x)
}
#Get sample names from barcodes
samples=function(x){
  x=(do.call(rbind, strsplit(x,"\\-|\\."))[,1:4])
  x=apply(x,1,function(x) paste(c(x[1:3],substr(x[4],1,2)),collapse="-"))
  return(x)
}

#sample types
sample_types=function(x){
  x=(do.call(rbind, strsplit(x,"\\-|\\."))[,4])
  x=substr(x,1,2)
  return(unique(x))
}

#generate network out of similarity matrix
make_nw=function(W){


  normalize <- function(X) X / rowSums(X)
  diag(W) = 0
  W = normalize(W);
  W = W + t(W);
  W[upper.tri(W)] = NA
  x=melt(as.matrix(W))
  Nw=x[!is.na(x$value),]
  Nw=Nw[Nw[,1]!=Nw[,2],]

}

###get common TCGA patients from the SNF input list
common_patients=function(features_list){
  t=lapply(features_list,rownames)
  return(Reduce("intersect",t))
}

plot_surv=function(surv,title="",xlab="",xscale=1){
  dist_group=sort(unique(surv$groups))#distinct groups
  surv=surv
  cols=c("red","green","blue","yellow","orange")
  legends=c("group1","group2","group3","group4","group5")
  fit=survfit(Surv(Survival)~groups,data=surv)
  plot(fit,col=cols[dist_group],main=title,cex.lab=2,cex.axis=1.8,xlab=xlab,lwd=4,xscale=xscale)
  sdf=survdiff(Surv(Survival,Death==1)~groups,data=surv)
  p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  legend("topright",legend=legends[dist_group],col=cols[dist_group],pch=19,cex=1.8)
  text(40,0.8,paste("p-value=",round(p.val,4),sep=""),cex=2)
}

compare=function(file1,file2,name,sampling=0){
  l1=deparse(substitute(file1))
  l2=deparse(substitute(file2))
  rownames(file2)=patients(rownames(file2))
  rownames(file1)=patients(rownames(file1))

  introws=rownames(file1)[rownames(file1) %in% rownames(file2)]
  intcols=colnames(file1)[colnames(file1) %in% colnames(file2)]

  x=c(file1[introws,intcols])#all common datapoints in file2 set
  y=c(file2[introws,intcols])#all common datapoints in file1 set

  if(sampling!=0){
    rsample=sample(1:length(x),sampling)#plot only a subset for performance reasons

  plot(x[rsample],y[rsample],main=name,xlab=l1,ylab=l2)
  legend("topleft",legend=paste("Rsq=",round(summary(lm(x[rsample]~y[rsample ]))[[8]],3)))}
  else{
    plot(x ,y,main=name,xlab=l1,ylab=l2)
    legend("topleft",legend=paste("Rsq=",round(summary(lm(x~y))[[8]],3)))
  }


  legend("bottomright",c(paste(length(introws),"common"),paste(length(rownames(file1)),l1),paste(length(rownames(file2)),l2)))
  abline(0,1)
}

format_ranks=function(features,ranks){

##This function formats the feature ranks in a tabular form and assigns the correct feature names to the scores
  features_ranked=list()
  for(i in seq(length(features))){
  features_ranked[[i]] = as.data.frame(cbind(ranks[[1]][[i]],as.integer(ranks[[2]][[i]])))
  colnames(features_ranked[[i]]) = c("score","ranking")
  rownames(features_ranked[[i]]) = colnames(features[[i]])
  features_ranked[[i]] = features_ranked[[i]][order(features_ranked[[i]]$ranking),]
  }
  return(features_ranked)

}

##This function is post processing the Similarity network for better visualization
post_process=function(W){
  x=W

  for(i in seq(dim(x)[1])){x[i,i]=NA}
  x=t(apply(x,1, function(x) {x[x<mean(x,na.rm=T)]=NA
                              return(x)}))
  x[is.na(x)]=0
  xNorm=(x+t(x))/2
  return(xNorm)
}


## Arguments:
## data: a list, where each item in the list is a matrix of values for each data type
## W: the target network for which the NMI is calculated against for each feature
##ncl:number of clusters
## Details:
## NMI is calculated based on the clustering assignments using spectral clustering
## The number of clusters is set based on the estimateNumberOfClustersGivenGraph on the target matrix
## using default parameters.
##
## Outputs:
## A list that contains the NMI score for each feature and their ranks from highest to lowest
## output[[1]] is the NMI score
## output[[1]][[1]] is the NMI score of first data type
## output[[1]][[1]][1] is the NMI score of the first feature of the first data type
## similarly for output[[2]]... except it is the rank instead of the score

rankFeaturesByNMI.2 <- function(data_list, W,ncl=4,cl=cl,num_of_clusters)
{
  stopifnot(class(data_list) == "list")

  feature_ranks=function(data,clustering_fused,num_of_clusters_fused){
    scores=vector(mode="numeric", length=dim(data)[2])
    for (feature_ind in seq(dim(data)[2]))
    {
      affinity_matrix <- affinityMatrix(
        dist2(as.matrix(data[, feature_ind]), as.matrix(data[, feature_ind])))
      clustering_single_feature <- spectralClustering(affinity_matrix, num_of_clusters_fused)
      scores[feature_ind] <- calNMI(clustering_fused, clustering_single_feature)
    }

    return(scores)
  }

  NUM_OF_DATA_TYES <- length(data_list)
  NMI_scores <- vector(mode="list", length=NUM_OF_DATA_TYES)
  NMI_ranks <- vector(mode="list", length=NUM_OF_DATA_TYES)


  clustering_fused <- spectralClustering(W, num_of_clusters)

  #clusterCall(cl, function() library(SNFtool))
  for (data_type_ind in 1:NUM_OF_DATA_TYES)
  {
    NUM_OF_FEATURES <- dim(data_list[[data_type_ind]])[2]
    NMI_scores[[data_type_ind]] <- vector(mode="numeric", length=NUM_OF_FEATURES)
    seq_colnames=colnames(data_list[[data_type_ind]])
    chunk_size=ceiling(length(seq_colnames)/ncl)#size of feature chunks that will be processed in diferent clusters
    col_split=split(seq_colnames, ceiling(seq_along(seq_colnames)/chunk_size))#split colnames in different chunks keeping the order
    data_split=lapply(col_split,function(x,data) data[,x],data_list[[data_type_ind]])# split feature in different chunks
    NMI_scores[[data_type_ind]]=do.call(c, parLapply(cl,  data_split,feature_ranks,clustering_fused, num_of_clusters))
    NMI_ranks[[data_type_ind]] <- rank(-NMI_scores[[data_type_ind]], ties.method="first")
  }
  stopCluster(cl)
  return(list(NMI_scores, NMI_ranks))
}