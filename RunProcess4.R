library(dplyr)
library(Seurat)
library(ggplot2)
library(Rcpp)
library(sctransform)
library(SCnorm)
library(scran)
library(cluster)
library(biomaRt)
sourceCpp('/home/node01/seurat normalization test/scaling2.cpp')
sourceCpp('/home/node01/seurat normalization test/neighbors.cpp')
source('/home/node01/seurat normalization test/scalingfactor_4.R')





#generate data
a<-read.table('/home/node01/seurat normalization test/data/transformed/raw/xin_transformed.txt',stringsAsFactors = FALSE)
countdata<-a
label<-read.table('/home/node01/seurat normalization test/data/transformed/raw/xin_labels.txt',stringsAsFactors = FALSE)
label=as.character(label[,1])

# If data is raw_count
# compare the number of genes with gene length in grch37 archive and grch38 and automatically choose the db with higher number of common genes.
# You can manually assign version(37 or 38)
ensembl<-findmart(countdata, species='human', vers=NULL)

                                


# 
# #basic tutorial
# # scaled 는 rpkm인지 umi count인지 none(raw count)인지 적어주는거에요.
# project<-processtoseurat(countdata = countdata, projectname = 'yan_embryo',norm.method='scnorm', neighbor.method = 'No jaccard', jaccard_qthreshold = 'none', scaled='none')
# project<-runFinal(project,npcs=30)
# 
# # if you have cell pre-annotation vector
# cellannotation=c(rep('Oocyte',3),rep('Zygote',3),rep('2cell',6), rep('4cell',12), rep('8cell',20), rep('Morulae',16),rep('Late.blastocyst',30),rep('hESC.passage.0',8),rep('hESC.passage.10',26))
# # if you have preprocessed yan embryo and cells are cut out.
# if(project@misc$cutoutcells!='none'){
#   cellannotation<-cellannotation[-project@misc$cutoutcells]
# }
# 
# # jaccard_qthreshold default=0.3
# # if you don't have cell pre-annotation vector, cellannot default is null.
# # if you use cellscaling with jaccard, norm.method should be previous norm.method+'_'+neighbor.method
# # if you don't use cellscaling with jaccard, jaccard_qthreshold is not used at the result session. It's just for printing plots.
# Result(project, cellannot = cellannotation, norm.method='sctransform' ,jaccard_qthreshold = 0.3, method='tsne')
# Result(project, cellannot = cellannotation, norm.method='sctransform' ,jaccard_qthreshold = 0.3, method='umap')













tsne_res<-vector()
umap_res<-vector()
pca_res<-vector()
rname<-vector()

#for other previous methods
#cell scaling 할때만 biomart를 사용해줬습니다.
for(norm.methods in c( 'scran','sctransform','scnorm','lognormalize')){
  project<-processtoseurat(countdata = countdata, projectname = 'xin_pancreas', norm.method=norm.methods, neighbor.method = 'No jaccard', jaccard_qthreshold = 'none', scaled='none')
  project<-runFinal(project,npcs=30)
  
  # if you have cell pre-annotation vector
  cellannotation=label
  # if you have preprocessed yan embryo and cells are cut out.
  if(project@misc$cutoutcells!='none'){
    cellannotation<-cellannotation[-project@misc$cutoutcells]
  }
  
  # jaccard_qthreshold default=0.3
  # if you don't have cell pre-annotation vector, cellannot default is null.
  tsne_res<-append(tsne_res, Result(project, cellannot = cellannotation, norm.method=norm.methods ,jaccard_qthreshold = 0.3, method='tsne'))
  umap_res<-append(umap_res, Result(project, cellannot = cellannotation, norm.method=norm.methods ,jaccard_qthreshold = 0.3, method='umap'))
  pca_res<-append(pca_res, Result(project, cellannot = cellannotation, norm.method=norm.methods ,jaccard_qthreshold = 0.3, method='pca', dims=30))
  rname<-append(rname, norm.methods)
}
                                                                                                                                        

#for cell scaling so as to test various jaccard qthreshold
for(iter in c(1)){
  for(neighbor.method in c('jaccard')){
    for(jaccard_qthreshold in c(0.25)){
      project<-processtoseurat(countdata = countdata, projectname = 'xin_pancreas', norm.method='cellscaling', neighbor.method = neighbor.method, jaccard_qthreshold = jaccard_qthreshold, scaled='none', iterator = iter, mart=ensembl)
      project<-runFinal(project,npcs=30)
      
      # if you have cell pre-annotation vector
      cellannotation=label
      # if you have preprocessed yan embryo and cells are cut out.
      if(project@misc$cutoutcells!='none'){
        cellannotation<-cellannotation[-project@misc$cutoutcells]
      }
      
      # jaccard_qthreshold default=0.3
      # if you don't have cell pre-annotation vector, cellannot default is null.
      tsne_res<-append(tsne_res, Result(project, cellannot = cellannotation, norm.method=paste('cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='tsne' , iterator=iter))
      umap_res<-append(umap_res, Result(project, cellannot = cellannotation, norm.method=paste('cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='umap' , iterator=iter))
      pca_res<-append(pca_res, Result(project, cellannot = cellannotation, norm.method=paste('cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='pca',dims=30, iterator=iter))
      rname<-append(rname,paste(iter,'_iterated ', neighbor.method,' ', jaccard_qthreshold, 'quantile'))
    }
  }
}




cname<-c('tsne','umap','pca')                                                                                                                                                                                                                    
yanE<-matrix(c(tsne_res,umap_res,pca_res),ncol=3)
colnames(yanE)<-cname
rownames(yanE)<-rname
yanE<-round(yanE, digits=4)
write.table(yanE, file ="/home/node01/seurat normalization test/result/xinP.csv",col.names = T, row.names = T )




# 
# #create seuratobject
# project<-CreateSeuratObject(counts=countdata, project='insert name', min.cells = 3, min.features = 200)
# project<-preprocess(project, processed=F) #지금 preprocess를 진행하면 preannotation vector를 수정해줘야하는데 preprocess cutoff를 모르겠어서 일단 rna level에 따른 cutoff는 제외해뒀습니다. 스텝진행하지않는것으로해뒀습니다.
# #혹시 preprocess 진행하면 잘려나간 cell들의 index를 확인하고 cell pre annotation vector에서 제거해줘야합니다. processed=F면 preprocess진행됩니다.
# project<-normalization(project, method='cellscaling', neighbor.method='jaccard' , jaccard_qthreshold = 0.3, logcount=F) 
# project<-runFinal(project,npcs=30)
# 
# # if you have cell pre-annotation vector
# cellannotation=c(rep('Oocyte',3),rep('Zygote',3),rep('2cell',6), rep('4cell',12), rep('8cell',20), rep('Morulae',16),rep('Late.blastocyst',30),rep('hESC.passage.0',8),rep('hESC.passage.10',26))
# # if you have preprocessed yan embryo and cells are cut out.
# if(project@misc$cutoutcells!='none'){
#   cellannotation<-cellannotation[-project@misc$cutoutcells]
# }
# 
# # jaccard_qthreshold default=0.3
# # if you don't have cell pre-annotation vector, cellannot default is null.
# Result(project, cellannot = cellannotation, jaccard_qthreshold = 0.3, method='tsne')
# Result(project, cellannot = cellannotation, jaccard_qthreshold = 0.3, method='umap')













processtoseurat<-function(countdata, projectname='xin_p', norm.method='cellscaling', neighbor.method='jaccard', jaccard_qthreshold=0.3, scaled='none', iterator=1, mart=NULL){
  if(norm.method == 'scran'){
    
    tempproject<-CreateSeuratObject(counts=countdata, project=projectname, min.cells = 3, min.features = 200)
    tempproject<-preprocess(tempproject, processed=F) #지금 preprocess를 진행하면 preannotation vector를 수정해줘야하는데 preprocess cutoff를 모르겠어서 일단 rna level에 따른 cutoff는 제외해뒀습니다. 스텝진행하지않는것으로해뒀습니다.
    sce <- SingleCellExperiment(list(counts=as.matrix(GetAssayData(object = tempproject, slot='counts'))))
    sce <- normalize(sce)
    logcounts<-sce@assays$data$logcounts
    logcounts<-log1p(2^logcounts)
    sce@assays$data$logcounts<-logcounts
    project<-as.Seurat(sce, counts = "counts", data = "logcounts")
    project@project.name=projectname
    project@misc$cutoutcells<-tempproject@misc$cutoutcells
    
  }else if(norm.method == 'sctransform'){
    project<-CreateSeuratObject(counts=countdata, project=projectname, min.cells = 3, min.features = 200)
    project<-preprocess(project, processed=F) #지금 preprocess를 진행하면 preannotation vector를 수정해줘야하는데 preprocess cutoff를 모르겠어서 일단 rna level에 따른 cutoff는 제외해뒀습니다. 스텝진행하지않는것으로해뒀습니다.
    project<-SCTransform(project)
  }else if(norm.method == 'scnorm'){
    tempproject<-CreateSeuratObject(counts=countdata, project=projectname, min.cells = 3, min.features = 200)
    tempproject<-preprocess(tempproject, processed=F) #지금 preprocess를 진행하면 preannotation vector를 수정해줘야하는데 preprocess cutoff를 모르겠어서 일단 rna level에 따른 cutoff는 제외해뒀습니다. 스텝진행하지않는것으로해뒀습니다.
    sce <- SingleCellExperiment(list(counts=as.matrix(GetAssayData(object = tempproject, slot='counts'))))
    Conditions = rep(c(1), each= ncol(sce@assays$data$counts))
    DataNorm<-SCnorm(Data = sce, Conditions = Conditions, PrintProgressPlots = TRUE, FilterCellNum = 10, K =1, NCores=1, reportSF = TRUE)
    project<-as.Seurat(DataNorm, counts = "normcounts", data = "normcounts")
    project@project.name=projectname
    #neighbor.method='No jaccard'
    #norm.method='scnorm'
    project<-normalization(project, method=norm.method, neighbor.method = neighbor.method, jaccard_qthreshold = jaccard_qthreshold, logcount=F, scaled=scaled)
    project@misc$cutoutcells<-tempproject@misc$cutoutcells
  }else{
    project<-CreateSeuratObject(counts=countdata, project=projectname, min.cells = 3, min.features = 200)
    project<-preprocess(project, processed=F)
    #혹시 preprocess 진행하면 잘려나간 cell들의 index를 확인하고 cell pre annotation vector에서 제거해주도록했습니다. 이 index는 project@misc에 넣었습니다. processed=F면 preprocess진행됩니다.

    project<-normalization(project, method=norm.method, scaled=scaled, neighbor.method=neighbor.method , jaccard_qthreshold = jaccard_qthreshold, logcount=F, iterator=iterator, mart=mart) 
  }
  return(project)
}





preprocess<-function(project, processed=F){
  mito.genes <- grep(pattern = "^MT-", x = rownames(project@assays[["RNA"]]), value = TRUE)
  percent.mito <- Matrix::colSums(project@assays[["RNA"]][mito.genes, ])/Matrix::colSums(project@assays[["RNA"]])
  project <- AddMetaData(object = project, metadata = percent.mito, col.name = "percent.mito") 
  project$percent.mito <- percent.mito
  if(!processed){
    #cutoutcells<-(which(!(project$nFeature_RNA > quantile(project$nFeature_RNA,0.001) & project$nFeature_RNA < quantile(project$nFeature_RNA,0.999) & project$percent.mito >  -Inf & project$percent.mito < 0.05)))
    cutoutcells<-which(!(project$percent.mito >  -Inf & project$percent.mito < 0.05))
    if(length(cutoutcells)>0){
      project@misc$cutoutcells=cutoutcells
    }else{
      project@misc$cutoutcells='none'
    }
    # project <- subset(x = project, subset = nFeature_RNA > quantile(project$nFeature_RNA,0.001) & nFeature_RNA < quantile(project$nFeature_RNA,0.999) & percent.mito >  -Inf & percent.mito < 0.05 )
    project <- subset(x = project, subset = percent.mito >  -Inf & percent.mito < 0.05 )
  }
  return(project)
}


# method=c('lognormalize','cellscaling','none)
# scaled=c('none', 'UMI', 'RPKM', 'TPM', 'FPKM')
# iterator is only for jaccard cell scaling
normalization<-function(project, scaled='none', method='lognormalize', neighbor.method='jaccard' , jaccard_qthreshold=0.3, logcount=FALSE, iterator=1, mart=NULL){
  if(method=='scran'){
    logcount=T
  }
  if(method=='lognormalize'){
    if(scaled %in% c('RPKM','TPM','FPKM', 'RPM')){
      countunnorm<-as.matrix(GetAssayData(object = project, slot='counts'))
      countnorm<-log1p(countunnorm)
      project<-SetAssayData(object=project, slot='data', new.data = countnorm)
    }else{
      NormalizeData(object = project, normalization.method = "LogNormalize", scale.factor = 10000)
    }
  }else if(method=='cellscaling'){
    countunnorm<-as.matrix(GetAssayData(object = project, slot='counts'))
    count_normalizing<-countunnorm
    libsize<-colSums(count_normalizing)
    if(scaled=='none'){
      
      g_inter<-intersect(rownames(count_normalizing), mart$external_gene_name)
      inter_m_index<-match(g_inter,mart$external_gene_name)
      mart<-mart[inter_m_index,]
      count_normalizing<-count_normalizing[match(g_inter, rownames(count_normalizing)),]
      cutoutcells<-project@misc$cutoutcells
      project <- CreateSeuratObject(count_normalizing)
      project@misc$cutoutcells<-cutoutcells
      count_normalizing<-(count_normalizing*(10^3))/mart$gene_length
      
    }
    for(i in 1:iterator){
      if(scaled%in%c('UMI','RPKM','TPM','FPKM','RPM')){
        jaccard_neighbors<-jac_neighbors(count_normalizing, method=neighbor.method)
        jaccard_threshold=quantile(unlist(jaccard_neighbors),jaccard_qthreshold)
        neighbors<-lapply(jaccard_neighbors, function(x)which(x>jaccard_threshold))
        norm.factor<-scalingfactor(count_normalizing, neighbors = neighbors)
      }else if(scaled=='none'){
        jaccard_neighbors<-jac_neighbors(count_normalizing, method=neighbor.method)
        jaccard_threshold=quantile(unlist(jaccard_neighbors),jaccard_qthreshold)
        neighbors<-lapply(jaccard_neighbors, function(x)which(x>jaccard_threshold))
        norm.factor<-scalingfactor(count_normalizing, neighbors = neighbors)
      }
      count_normalizing<-t(t(count_normalizing)*(norm.factor))
    }
    if(scaled%in%c('RPKM','TPM','FPKM', 'RPM')){
    }else{
      count_normalizing<-t(t(count_normalizing)/(libsize))
      count_normalizing<-count_normalizing*10^4
    }
    countnorm<-log1p(count_normalizing)
    
    countnorm <- as(object = countnorm, Class = "dgCMatrix")
    project<-SetAssayData(object=project, slot='data', new.data = countnorm)
    # 
  }else if(logcount==T){
    countunnorm<-as.matrix(GetAssayData(object = project, slot='counts'))
    countnorm <- as(object = countunnorm, Class = "dgCMatrix")
    project<-SetAssayData(object=project, slot='data', new.data = countnorm)
  }else{
    countunnorm<-as.matrix(GetAssayData(object = project, slot='counts'))
    countnorm<-(countunnorm/colSums(countunnorm))*10^4
    countnorm<-log1p(countnorm)
    
    countnorm <- as(object = countnorm, Class = "dgCMatrix")
    project<-SetAssayData(object=project, slot='data', new.data = countnorm)
  }
  return(project)
}


double_normalization<-function(project, scaled='none', first.method='lognormalize', second.neighbor.method='jaccard' , jaccard_qthreshold=0.3, logcount=FALSE, iterator=1, mart=NULL){
  if(method=='scran'){
    logcount=T
  }
  if(method=='lognormalize'){
    if(scaled %in% c('RPKM','TPM','FPKM', 'RPM')){
      countunnorm<-as.matrix(GetAssayData(object = project, slot='counts'))
      countnorm<-log1p(countunnorm)
      project<-SetAssayData(object=project, slot='data', new.data = countnorm)
    }else{
      NormalizeData(object = project, normalization.method = "LogNormalize", scale.factor = 10000)
    }
}



runFinal<-function(project, npcs=30){
  
  project <- FindVariableFeatures(object = project, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
  project <- ScaleData(object = project, vars.to.regress = c("nCount_RNA", "percent.mito"))
  project <- RunPCA(object = project,  npcs = npcs, verbose = FALSE)
  
  project <- FindNeighbors(project, reduction = "pca", dims = 1:20)
  project <- FindClusters(project, resolution = 0.5, algorithm = 1)
  project <- RunTSNE(object = project, dims.use = 1:20, do.fast = F)
  project <- RunUMAP(project, dims = 1:20)
  
  return(project)
}


Result<-function(project, jaccard_qthreshold=0.3, cellannot=NULL, norm.method='cellscaling', method='tsne', dims=2 ,iterator=1){
  if(!is.null(cellannot)){
    project$seurat_clusters<-cellannot
    project$seurat_clusters<-factor(cellannot)
    project@active.ident<-project$seurat_clusters
  }
  
  dist.matrix <- dist(x = Embeddings(object = project[[method]])[, 1:dims])
  clusters <- project@active.ident
  sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
  silscore<-mean(sil[,3])
  
  
  
  png(paste('/home/node01/seurat normalization test/result/',project@project.name,norm.method,method,jaccard_qthreshold*100,'percent quantile.png',sep=' '))
  print(DimPlot(project, reduction = method, pt.size = 2) + labs(title=paste(project@project.name, method, sep='_')) + labs(subtitle = if(norm.method=='cellscaling_jaccard'||norm.method=='cellscaling_generalized_jaccard')paste(norm.method, iter, 'iterated' , 'jaccard threshold :', jaccard_qthreshold*100 ,'% quantile \n silhouette score :', round(silscore,digits = 4), sep=' ')else paste(norm.method ,'\n silhouette score :', round(silscore,digits = 4), sep=' ') ))
  dev.off()
  
  return(silscore)
}

findmart<-function(countdata, species, vers=NULL){
  dataset<-switch(species, 'human'="hsapiens_gene_ensembl", 'mouse'='mmusculus_gene_ensembl')
  if(is.null(vers)){
    ensembl37<-useMart('ensembl', dataset=dataset, host='grch37.ensembl.org')
    ensembl38<-useMart('ensembl', dataset=dataset)
    ensembl37 <-getBM(mart=ensembl37,attributes = c('external_gene_name','start_position','end_position'))
    ensembl38 <-getBM(mart=ensembl38,attributes = c('external_gene_name','start_position','end_position'))
    #if annotated gene names have more common in grch37, return ensembl dataset from grch37 annotation  
    
    if(length(intersect(rownames(countdata),ensembl37$external_gene_name))>length(intersect(rownames(countdata),ensembl38$external_gene_name))){
      ensembl<-ensembl37
    }else{
      ensembl<-ensembl38
    }
  }else if(vers==37){
    ensembl37<-useMart('ensembl', dataset=dataset, host='grch37.ensembl.org')
    ensembl37 <-getBM(mart=ensembl37,attributes = c('external_gene_name','start_position','end_position'))
    ensembl<-ensembl37
  }else if(vers==38){
    ensembl38<-useMart('ensembl', dataset=dataset)
    ensembl38 <-getBM(mart=ensembl38,attributes = c('external_gene_name','start_position','end_position'))
    ensembl<-ensembl38
  }else{
    print('unknown vers')
    return(NULL)
  }
  ensembl<-ensembl[which(ensembl$start_position!="" & ensembl$end_position!=""),]
  gene_length<-(ensembl$end_position-ensembl$start_position)+1
  ensembl <-cbind(ensembl,gene_length)
  ensembl<-ensembl[,c(1,4)]
  ensembl_max<-vector()
  for(i in unique(ensembl$external_gene_name)){
    ind<-which(ensembl$external_gene_name==i)
    ensembl_max<-append(ensembl_max, max(ensembl$gene_length[ind]))
  }
  ensembl_ind<-match(unique(ensembl$external_gene_name),ensembl$external_gene_name)
  ensembl<-ensembl[ensembl_ind,]
  ensembl$gene_length<-ensembl_max
  return(ensembl)
}


# Result<-function(project, jaccard_qthreshold=0.3, cellannot=NULL){
#   if(!is.null(cellannot)){
#     project$seurat_clusters<-cellannot
#     project$seurat_clusters<-factor(cellannot)
#     project@active.ident<-project$seurat_clusters
#   }
#   
#   dist.matrix_tsne <- dist(x = Embeddings(object = project[['tsne']])[, 1:2])
#   clusters <- project@active.ident
#   siltsne <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
#   silscore_tsne<-mean(siltsne[,3])
#   
#   dist.matrix_umap <- dist(x = Embeddings(object = project[['umap']])[, 1:2])
#   silumap <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
#   silscore_umap<-mean(silumap[,3])
#   
#   # png(paste('/home/node01/seurat normalization test/result/',project@project.name,'tsne',jaccard_qthreshold*100,'percent quantile','.png',sep=''))
#   DimPlot(project, reduction = "tsne", pt.size = 2) + labs(title=paste(project@project.name, 'tsne', sep='_')) + labs(subtitle = paste('jaccard threshold :', jaccard_qthreshold*100 ,'% quantile \n silhouette score :', round(silscore_tsne,digits = 4), sep=' ') )
#   # dev.off()
#   # png(paste('/home/node01/seurat normalization test/result/',project@project.name,'umap',jaccard_qthreshold*100,'percent quantile','.png',sep=''))
#   DimPlot(project, reduction = "umap", pt.size = 2) + labs(title=paste(project@project.name, 'umap', sep='_')) + labs(subtitle = paste('jaccard threshold :', jaccard_qthreshold*100 ,'% quantile \n silhouette score :', round(silscore_umap,digits=4), sep=' ') )
#   # dev.off()
# }
# 
