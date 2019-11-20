library(dplyr)
library(Seurat)
library(ggplot2)
library(Rcpp)
library(sctransform)
library(SCnorm)
library(scran)
library(cluster)
library(biomaRt)
sourceCpp('scaling2.cpp')
sourceCpp('neighbors.cpp')
source('scalingfactor_4.R')
source('basefunctions.R')




res_dir='/home/node01/seurat normalization test/result/'
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
projectname='xin_pancreas'

#for other previous methods
#cell scaling 할때만 biomart를 사용해줬습니다.
for(norm.methods in c( 'scran','sctransform','scnorm','lognormalize')){
  project<-processtoseurat(countdata = countdata, projectname = projectname, norm.method=norm.methods, neighbor.method = 'No jaccard', jaccard_qthreshold = 'none', scaled='none')
  project<-runFinal(project,npcs=30)
  
  # if you have cell pre-annotation vector
  cellannotation=label
  # if you have preprocessed yan embryo and cells are cut out.
  if(project@misc$cutoutcells!='none'){
    cellannotation<-cellannotation[-project@misc$cutoutcells]
  }
  
  # jaccard_qthreshold default=0.3
  # if you don't have cell pre-annotation vector, cellannot default is null.
  tsne_res<-append(tsne_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=norm.methods ,jaccard_qthreshold = 0.3, method='tsne'))
  umap_res<-append(umap_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=norm.methods ,jaccard_qthreshold = 0.3, method='umap'))
  pca_res<-append(pca_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=norm.methods ,jaccard_qthreshold = 0.3, method='pca', dims=30))
  rname<-append(rname, norm.methods)
}
                                                                                                                                        

#for cell scaling so as to test various jaccard qthreshold
for(iter in c(1,3,5)){
  for(neighbor.method in c('jaccard', "generalized_jaccard")){
    for(jaccard_qthreshold in c(0.25,0.5,0.75,0.85)){
      project<-processtoseurat(countdata = countdata, projectname = projectname, norm.method='cellscaling', neighbor.method = neighbor.method, jaccard_qthreshold = jaccard_qthreshold, scaled='none', iterator = iter, mart=ensembl)
      project<-runFinal(project,npcs=30)
      
      # if you have cell pre-annotation vector
      cellannotation=label
      # if you have preprocessed yan embryo and cells are cut out.
      if(project@misc$cutoutcells!='none'){
        cellannotation<-cellannotation[-project@misc$cutoutcells]
      }
      
      # jaccard_qthreshold default=0.3
      # if you don't have cell pre-annotation vector, cellannot default is null.
      tsne_res<-append(tsne_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=paste('cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='tsne' , iterator=iter))
      umap_res<-append(umap_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=paste('cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='umap' , iterator=iter))
      pca_res<-append(pca_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=paste('cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='pca',dims=30, iterator=iter))
      rname<-append(rname,paste(iter,'_iterated ', neighbor.method,' ', jaccard_qthreshold, 'quantile'))
    }
  }
}

#for double scaling
for(iter in c(1,3,5)){
  for(neighbor.method in c('jaccard', "generalized_jaccard")){
    for(jaccard_qthreshold in c(0.25,0.5,0.75,0.85)){
      for(norm.methods in c( 'scran','sctransform','scnorm','lognormalize')){
        print(paste(iter, neighbor.method,jaccard_qthreshold,norm.methods,sep=' '))
        project<-processtoseurat(countdata = countdata, projectname = projectname, norm.method=norm.methods, neighbor.method = neighbor.method, jaccard_qthreshold = jaccard_qthreshold, scaled='none', iterator = iter, mart=ensembl, doublenorm = T)
        project<-runFinal(project,npcs=30)
        
        # if you have cell pre-annotation vector
        cellannotation=label
        # if you have preprocessed yan embryo and cells are cut out.
        if(project@misc$cutoutcells!='none'){
          cellannotation<-cellannotation[-project@misc$cutoutcells]
        }
        
        # jaccard_qthreshold default=0.3
        # if you don't have cell pre-annotation vector, cellannot default is null.
        tsne_res<-append(tsne_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=paste(norm.methods,'cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='tsne' , iterator=iter))
        umap_res<-append(umap_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=paste(norm.methods,'cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='umap' , iterator=iter))
        pca_res<-append(pca_res, Result(project, result.dir=res_dir, cellannot = cellannotation, norm.method=paste(norm.methods,'cellscaling',neighbor.method,sep='_') ,jaccard_qthreshold = jaccard_qthreshold, method='pca',dims=30, iterator=iter))
        rname<-append(rname,paste(norm.methods, iter,'_iterated ', neighbor.method,' ', jaccard_qthreshold, 'quantile'))
      }
    }
  }
}



cname<-c('tsne','umap','pca')                                                                                                                                                                                                                    
yanE<-matrix(c(tsne_res,umap_res,pca_res),ncol=3)
colnames(yanE)<-cname
rownames(yanE)<-rname
yanE<-round(yanE, digits=4)
write.table(yanE, file =paste(res_dir,projectname,'.csv',sep=''),col.names = T, row.names = T )




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











