

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

 
#generate network out of similarity matrix
make_nw=function(W){


  normalize <- function(X) X / rowSums(X)
  #assign diagonal value to 0
  diag(W) = 0
  W = normalize(W);
  W = W + t(W);
  W[upper.tri(W)] = NA
  x=melt(as.matrix(W))
  Nw=x[!is.na(x$value),]
  Nw=Nw[Nw[,1]!=Nw[,2],]

}

###get common patient names from the tables in input list
common_patients=function(features_list){
  t=lapply(features_list,rownames)
  return(Reduce("intersect",t))
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

 