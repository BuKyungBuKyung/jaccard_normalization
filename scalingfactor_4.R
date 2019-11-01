sourceCpp('/home/node01/seurat normalization test/neighbors.cpp')
sourceCpp('/home/node01/seurat normalization test/scaling2.cpp')
jac_neighbors<-function(obs, lib.size=NULL, method='generalized_jaccard'){ #지금은jaccard threshold가manual한값
  # obs<-countunnorm
  
  
  x <- as.matrix(obs)
  if(any(is.na(x))) stop("NA counts not permitted")
  nsamples <- ncol(x)
  
  if(is.null(lib.size)) {
    lib.size <- colSums(x)
  } else {
    if(anyNA(lib.size)) stop("NA lib.sizes not permitted")
    if(length(lib.size) != nsamples) {
      if(length(lib.size) > 1L) warning("calcNormFactors: length(lib.size) doesn't match number of samples",call.=FALSE)
      lib.size <- rep(lib.size,length=nsamples)
    }
  }
  
  #	Remove all zero rows
  allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
  if(any(allzero)) x <- x[!allzero,,drop=FALSE]
  
  
  non_zero<-list()
  for(i in 1:ncol(x)) non_zero[[i]] <- which(x[,i] != 0L)
  # sourceCpp('/home/node01/seurat normalization test/neighbors.cpp')
  jaccard_neighbors<-list()
  # jaccard_neighborsrj<-r_jaccard_neighbors(non_zero, ncol(x))
  # jaccard_neighborsrgj<-r_generalized_jaccard_neighbors(x, non_zero)
  # jaccard_neighborscj<-rcpp_jaccard_neighbors(non_zero, ncol(x))
  # jaccard_neighborscgj<-rcpp_generalized_jaccard_neighbors(x, non_zero, ncol(x))
  # jaccard_neighbors<-rcpp_jaccard_neighbors(non_zero, ncol(x))
  if(method=='generalized_jaccard'){
    jaccard_neighbors<-rcpp_generalized_jaccard_neighbors(x, non_zero, ncol(x))
  }else if(method=='jaccard'){
    jaccard_neighbors<-rcpp_jaccard_neighbors(non_zero, ncol(x))
  }
  # jaccard_neighbors<-rcpp_generalized_jaccard_neighbors(x, non_zero, ncol(x))
  
  return(jaccard_neighbors)
}


r_jaccard_neighbors<-function(non_zero, colx){
  jaccard_neighbors=list()
  for(i in 1:colx){
    # print(paste(i, 'jaccard', sep=' '))
    
    temp_jaccard_neighbor<-vector()
    for(j in 1:colx){
      temp_jaccard_neighbor<-append(temp_jaccard_neighbor, jaccard_similarity(non_zero[[i]],non_zero[[j]]))
    }
    jaccard_neighbors[[i]]<-temp_jaccard_neighbor
  }
  return(jaccard_neighbors)
}

r_generalized_jaccard_neighbors<-function(x,non_zero){
  jaccard_neighbors=list()
  for(i in 1:ncol(x)){
    # print(paste(i, 'jaccard', sep=' '))
    temp_jaccard_neighbor<-vector()
    for(j in 1:ncol(x)){
      inter=intersect(non_zero[[i]],non_zero[[j]])
      print(length(inter))
      temp_jaccard_neighbor<-append(temp_jaccard_neighbor, (sum(x[inter,i])+sum(x[inter,j]))/(sum(x[,i])+sum(x[,j])))
    }
    jaccard_neighbors[[i]]<-temp_jaccard_neighbor
  }
  return(jaccard_neighbors)
}



scalingfactor<-function(obs, lib.size=NULL, neighbors){
  # obs<-countunnorm
  # lib.size=NULL
  # neighbors=neighbors
  x <- as.matrix(obs)
  if(any(is.na(x))) stop("NA counts not permitted")
  nsamples <- ncol(x)

  if(is.null(lib.size)) {
    lib.size <- colSums(x)
  } else {
    if(anyNA(lib.size)) stop("NA lib.sizes not permitted")
    if(length(lib.size) != nsamples) {
      if(length(lib.size) > 1L) warning("calcNormFactors: length(lib.size) doesn't match number of samples",call.=FALSE)
      lib.size <- rep(lib.size,length=nsamples)
    }
  }

  #	Remove all zero rows
  allzero <- .rowSums(x>0, nrow(x), nsamples) == 0L
  if(any(allzero)) x <- x[!allzero,,drop=FALSE]



  #n neighbor로바뀐다면
  scf<-list()
  #scf<-rscale(x, neighbors)
  scf<-rcpp_scale(x,neighbors, ncol(x))

  
  #아래부분은 차후 rcpp scale이나 rscale함수에 통합시킬것.
  scf_n<-list()
  for(i in 1:length(scf)){
    scf_n[[i]]<-scf[[i]]
    scf_n[[i]]<-2^mean(log2(scf_n[[i]]))
  }
  scf_n<-unlist(scf_n)
  scf_n
}


rscale<-function(x,neighbors){
  scf<-list()
  for(i in 1:ncol(x)){
    # print(paste(i, 'neighbors', sep=' '))
    nodescf<-vector()
    for(j in neighbors[[i]]){
      obs<-as.numeric(x[,j])
      ref<-as.numeric(x[,i])
      nO <- sum(obs)
      nR <- sum(ref)                           # 여기까지만 실행시켰을때 rscale 시간 178ms, c에서 numeric vector 등으로 여기까지실행시켰을때 걸린시간 38ms.
      logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size

      fin<-is.finite(logR)
      f <- mean(logR[fin], na.rm=TRUE)         #r에서는 narm으로 0으로 나눠주거나 0을 나눠준 값을 걸러낼수있어요. 그래서 0 index를 따로구하거나 할필요가 없는데 c에서는 직접구해줘야했어요..
      nodescf<-append(nodescf, sqrt(2^f))
      # nodescf<-append(nodescf, nO)
    }
    scf[[i]]<-nodescf
  }
  return(scf)
}



# 
# rscale2<-function(x,neighbors){
#   scf<-list()
#   for(i in 1:ncol(x)){
#     # print(paste(i, 'neighbors', sep=' '))
#     nodescf<-vector()
#     for(j in neighbors[[i]]){
#       obs<-as.numeric(x[,j])
#       ref<-as.numeric(x[,i])
#       nO <- sum(obs)
#       nR <- sum(ref)
#       logR <- log2(obs/ref)
#       # logR <- log2((obs/nO)/(ref/nR))          # log ratio of expression, accounting for library size
#       # 
#       fin<-is.finite(logR)
#       # f <- mean(logR[fin], na.rm=TRUE)
#       # nodescf<-append(nodescf, sqrt(2^f))
#       nodescf<-append(nodescf, nO)
#     }
#     scf[[i]]<-nodescf
#   }
#   return(scf)
# }





jaccard_similarity<-function(vec1, vec2){
  return(length(intersect(vec1, vec2))/length(union(vec1,vec2)))
}



