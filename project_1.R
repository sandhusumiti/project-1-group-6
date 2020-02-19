#installing required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("hgu133plus2.db"))
library(hgu133plus2.db)

data <- read.table("welch_results_allgenes.csv", head= TRUE, sep=",") 
head(data)
#test <- select(hgu133plus2.db, keys = data[,1], columns = ("SYMBOL"))
class(data[,1])
#dim(data)
data <-na.omit(data)
#dim(data)

#matching the probes to the gene symbols 
test <- select(hgu133plus2.db, keys = as.character(data[,1]), columns = ("SYMBOL"))
dim(test)
head(test)

test2 = test[ !duplicated(test[1]), ]
duplicated(test[1])


symbol_col <- test2$SYMBOL
data <- cbind(data, symbol_col)
data = data[ !duplicated(data$symbol_col), ]
dim(data)

write.csv(data,"6_2.csv")

###### getting the genesets from genesetcollections
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSEABase")
library(GSEABase)
gene_set_kegg<- getGmt("c2.cp.kegg.v7.0.symbols.gmt",
                       collectionType=BroadCollection(category="c3"),
                       geneIdType=SymbolIdentifier())
gene_set_kegg

gene_set_GO<- getGmt("c5.all.v7.0.symbols.gmt",
                     collectionType=BroadCollection(category="c3"),
                     geneIdType=SymbolIdentifier())
gene_set_GO

gene_set_Hallmark<- getGmt("h.all.v7.0.symbols.gmt",
                           collectionType=BroadCollection(category="c3"),
                           geneIdType=SymbolIdentifier())
gene_set_Hallmark

# genes from kegg
genes_kegg <- geneIds(gene_set_kegg)

# genes from GO
genes_GO <- geneIds(gene_set_kegg)

#genes from hallmark
genes_hallmark <- geneIds(gene_set_Hallmark)

#function for fisher's exact test


fisher <- function(geneset){ 
  new_data <- read.table(file="6_2.csv", TRUE, sep=",")
  #n_pathway <- dim(summary(geneset)) 
  n_pathway <- length(geneset)
  c=0
  p_val <- vector("numeric", length(geneset))
  est <-vector("numeric", length(geneset))
  for(pathway in c(1:n_pathway)){ 
    genes <- geneset[[pathway]]
    no_genes <-length(genes)
    #print(genes)
    #print(no_genes)
    gene_df <- data.frame(genes)
    #print(gene_df)
    #for (genes in c(1:no_genes)) {
    total_in_geneset<- sum(new_data$symbol_col %in% gene_df$genes)
    #print(total_in_geneset)
    not_in_geneset <- 1077-total_in_geneset
    #print(not_in_geneset)
    first <- sum(new_data$adj_p<0.05 & new_data$symbol_col %in% gene_df$genes)
    #print(first)
    second<- sum(new_data$adj_p<0.05 & !(new_data$symbol_col %in% gene_df$genes))
    #print(second)
    third <- total_in_geneset - first
    #print(third)
    fourth <- not_in_geneset - second
    #print(fourth)
    
    result <- fisher.test(matrix(c(first,second,third,fourth),nrow=2))
    
    #print(result$p.value)
    #print(result$estimate)
    p <-result$p.value
    e <- result$estimate
    p_val[pathway] <- p
    est[pathway]<- e
    #print(p_val)
    #print(est)
    c=c+1
    #print(c)
    
    data_geneset<-data.frame(names(geneset),p_val,est)
    #print(data_kegg)
    write.csv(data_geneset,"fishers_table.csv")
    #}
  }
  
}

#calling function for each geneset collection
#fisher(genes_kegg)
#fisher(genes_GO)
#fisher(genes_hallmark)

#reading data gathered from all three comparisons
data_kegg <- read.table(file="kegg.csv", TRUE, sep=",")
data_go <- read.table(file="GO.csv", TRUE, sep=",")
data_hallmark <- read.table(file="hallmark.csv", TRUE,sep=",")

#adding adjusted p-value (FDR) to the table containing results from fishers test
adj_p<-p.adjust(data_kegg$p_val,method = "fdr")
data_kegg <- cbind(data_kegg, adj_p)
write.csv(data_kegg,"kegg.csv")
adj_p<-p.adjust(data_go$p_val,method = "fdr")
data_go <- cbind(data_go, adj_p)
write.csv(data_go,"GO.csv")
adj_p<-p.adjust(data_hallmark$p_val,method = "fdr")
data_hallmark <- cbind(data_hallmark, adj_p)
write.csv(data_hallmark,"hallmark.csv")

#ordering by nominal p-value
data_kegg <-data_kegg[order(data_kegg$p_val),]
write.csv(data_kegg,"kegg.csv")

data_go <-data_go[order(data_go$p_val),]
write.csv(data_go,"GO.csv")
#head(data_go)

data_hallmark <-data_hallmark[order(data_hallmark$p_val),]
write.csv(data_go,"GO.csv")
#head(data_go)

#top three enriched genesets
e_kegg<-data_kegg[1:3,]
e_kegg <- e_kegg[ -c(1) ]
write.csv(e_kegg,"enriched_kegg.csv")
e_go<-data_go[1:3,]
e_go<- e_go[ -c(1) ]
write.csv(e_go,"enriched_GO.csv")
e_h<-data_hallmark[1:3,]
e_h<- e_h[ -c(1) ]
write.csv(e_h,"enriched_hallmark.csv")


#up and down-regulated genes
data <- read.table(file="t_test_filter_2.csv", TRUE, na.strings=c("NA", "NA"), sep=",")
data <-data[order(data$t_stat),]


down<-data[1:10,]
down<- down[ -c(1) ]
write.csv(down,"down_regulated_genes.csv")
up<-tail(data,n=10)
write.csv(up,"upregulated_gene.csv")