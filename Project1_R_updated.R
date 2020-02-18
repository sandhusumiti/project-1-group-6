if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)

data=ReadAffy(celfile.path = '/projectnb/bf528/users/group6/project_1/samples',compress = TRUE)
norm=rma(data)

fit=fitPLM(data, normalize=TRUE, background=TRUE)
RLE(fit, main='Microbiome dataset RLE')

RLE_median=RLE(fit,type='stats')
RLE_median=data.frame(RLE_median)
median_RLE=RLE_median[-c(2),]
median_RLE=as.numeric(median_RLE[1,])
hist(median_RLE, main='Median RLE', xlab='Median RLE')



NUSE_median=NUSE(fit,type='stats')
NUSE_median=data.frame(NUSE_median)
median_NUSE=NUSE_median[-c(2),]
median_NUSE=as.numeric(median_NUSE[1,])
hist(median_NUSE, main='Median NUSE', xlab='Median NUSE')

proj_metadata=read.csv('/project/bf528/project_1/doc/proj_metadata.csv',header=TRUE, stringsAsFactors = FALSE)
norm_data=exprs(norm)

batch_effect=proj_metadata$normalizationcombatbatch
#mod=proj_metadata$normalizationcombatmod
mod1=model.matrix(~proj_metadata$normalizationcombatmod)


ComBat(as.matrix(norm_data),batch_effect,mod=mod1)

Combat_data= ComBat(dat=norm_data,batch=batch_effect,mod=mod1)
write.csv(Combat_data, "expression_data.csv")

PCA_comp=t(Combat_data)
PCA_comp=scale(PCA_comp)
PCA_comp=t(PCA_comp)
PCA_comp=prcomp(PCA_comp, center=FALSE,scale=FALSE)

PCA_comp
prcomp_val=as.data.frame(PCA_comp$rotation)
prcomp_val
#plot(prcomp_val[,1], prcomp_val[,2], main="PCA plot", xlab="PC 1", ylab="PC 2",  col = )
#components <- as.data.frame(PCA_comp$x)
PCA_comp$importance

library(ggplot2)
summary(PCA_comp)

ggplot(data = prcomp_val, mapping = aes(x = PC1, y = PC2)) +
  geom_point() +  
  theme_bw() +
  labs(title="PCA plot")

ggplot(data = prcomp_val, mapping = aes(y = PC1)) +
  geom_boxplot() +  
  theme_bw() +
  labs(title="PC1 histogram")

ggplot(data = prcomp_val, mapping = aes(y = PC2)) +
  geom_boxplot() +  
  theme_bw() +
  labs(title="PC2 histogram")


id = which(!(prcomp_val$PC1 > mean(prcomp_val$PC1) + 3*sd(prcomp_val$PC1) | prcomp_val$PC1 < mean(prcomp_val$PC1) - 3*sd(prcomp_val$PC1) | 
               prcomp_val$PC2 > mean(prcomp_val$PC2) + 3*sd(prcomp_val$PC2) | prcomp_val$PC2 < mean(prcomp_val$PC2) - 3 * sd(prcomp_val$PC2))) 


#length(id) is 133

Combat_data.filtered <- Combat_data[id,]
write.csv(Combat_data.filtered, "expression_data.csv_updated")


PCA_comp_new=t(Combat_data.filtered)
PCA_comp_new=scale(PCA_comp_new)
PCA_comp_new=t(PCA_comp_new)
PCA_comp_new=prcomp(PCA_comp_new, center=FALSE,scale=FALSE)
summary(PCA_comp_new)

new_prcomp_val=as.data.frame(PCA_comp_new$rotation)
new_prcomp_val

ggplot(data = new_prcomp_val, mapping = aes(x = PC1, y = PC2)) +
  geom_point() +  
  theme_bw() +
  labs(title="New PCA plot")
