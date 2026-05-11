if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("AUCell")
library('GEOquery')
setwd("E:/MasterDegree/R/1")
gset = getGEO('GSE12417', destdir=".",getGPL = F)
print(gset)
View(gset)
class(gset)
names(gset)
e2=gset[[2]]
e2.1=gset[["GSE12417-GPL96_series_matrix.txt.gz" ]]
class(e2)
e2.2=gset[["GSE12417-GPL96_series_matrix.txt.gz"]]
gset <- getGEO('GSE12417', destdir=".",getGPL = F)
gset <- getGEO('GSE12417', destdir=".",getGPL = F)

save(gset,file = "gset.rdata")
getGEO('GSE12417', destdir=".",getGPL = F)
save(gset,file = "3.gse12417.rdata")
saveRDS(gset,file = "gse12417.rds")

View(gset)
anno=gset[["GSE12417-GPL96_series_matrix.txt.gz"]]@assayData[["exprs"]]
anno1=anno[c(1,2,3,4),c(1,2,3,4)]
print(anno1)
class(anno1)
library(data.table)
anno=fread("GPL96-57554.txt",sep = "\t",header = T,data.table = F)

View(anno)
colnames(anno)
gene=anno[,c(1,11)]
x1=gene$`Gene Symbol`
class(x1)

gene.1=anno[,c("ID","Gene Symbol")]
###exp是矩阵 merge合并必须是数据框的合并
exp1=as.data.frame(exp)
exp.anno=merge(x=anno1,y=exp1,by.x=1,by.y=0)
x=gene$`Gene Symbol`
a1=strsplit(x,split = " /// ",fixed = T)
print(a1[c(1,2,3,4,5,6,7)])
gene.all = sapply(a1,function(x){x[1]})
View(sapply(a1,function(x){x[2]}))
a3=data.frame(anno$ID,gene.all)
View(a3)
exp2=e2@assayData[["exprs"]]
View(exp2)
View(anno)
exp.merge=merge(x=a3,y=exp2,by.x=1,by.y=0)
View(exp.merge)
exp.distinct=distinct(exp.merge, gene.all, .keep_all = T)
exp3=na.omit(exp.distinct)
?na.omit
rownames(exp3)=exp3$gene.all
View(exp3)
exp3 <- exp3[,-c(1,2)]
View(exp3)
?distinct
