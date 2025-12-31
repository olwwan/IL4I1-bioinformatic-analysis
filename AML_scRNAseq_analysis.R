##数据转换 R4.3.0:Seurat v4.3.0（本次用） R4.4.0: Seurat v5.1.0 ；
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(R.utils)
library(data.table)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(patchwork)
library(devtools)
library(harmony) #install.packages
#library("clusterProfiler")
library(SingleCellExperiment)
#options(Seurat.object.assay.version = "v3")
library(SingleR)
#BiocManager::install("celldex")
library(celldex)
#BiocManager::install("BiocParallel")
library(BiocParallel)



###################################
###################################
###########1. dataset input#############
###################################
###################################
###################################

###1、读取10x代码：数据集4
dir='./rawdata/GSE154109_RAW/'
#每个样本在 一个文件夹下
samples=list.files(dir)
samples
project_id = gsub("_$","", substr(samples,1,10))
project_id
#抽屉AML样本
samples = samples[1:12] #AML:5-12, HC:1:4, ALL bone:13:19
project_id = project_id[1:12]
#合并数据
sce1 = list()
for(i in 1:length(samples)){ 
  pro = samples[i]
  print(pro)
  dir2 = paste0(dir,pro)
  sce1[[i]]=CreateSeuratObject(Read10X(file.path(dir,pro)) ,
                         project = project_id[i],
                         min.cells = 3,
                         min.features = 300)				 
  sce1[[i]][["percent.mt"]] <- PercentageFeatureSet(sce1[[i]], pattern = "^MT-") #线粒体基因
  sce1[[i]][["percent.ribo"]] <- PercentageFeatureSet(sce1[[i]], pattern = "^RP[SL]") #核糖体基因：以RPS或RPL开头的基因
  sce1[[i]][["percent.hb"]] <- PercentageFeatureSet(sce1[[i]], pattern = "^HB[^(P)]") #红细胞基因：除以HBP开头外的HB开头基因
}
head(sce1[[length(samples)]]@meta.data)
save(sce1,file = "1_GSE154109_AMLHCsamples_scRNA.Rdata")



###################################
###################################
###########2. scRNAseq analysis#############
###################################
###################################
###################################
rm(list=ls())

load("1_GSE154109_samples_scRNA.Rdata")

outdir="AML_GSE154109_Rout"

if(!dir.exists(outdir)) dir.create(outdir)
setwd(outdir)
scRNA_merge= sce1
scRNA_merge= merge(scRNA_merge[[1]], y=c(scRNA_merge[-1]))

# 质控结果的可视化
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo","percent.hb")
# 一般用小提琴图对结果进行质控
rp = VlnPlot(scRNA_merge, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3,) + NoLegend() #pt.size = 0.1显示点
pdf("11_raw_samples_features_before_filtration.pdf", height = 12, width = 17)
rp
dev.off()
# 基因数量随细胞数量的变化关系
FeatureScatter(scRNA_merge, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
# 线粒体基因数量随细胞数量的变化关系
FeatureScatter(scRNA_merge, "nCount_RNA", "percent.mt", group.by = "orig.ident", pt.size = 0.5)

#第6,10个样本质量太差，删除
scRNA_merge= sce1
scRNA_merge= merge(scRNA_merge[[1]], y=c(scRNA_merge[c(2:5,7:9,11:12)]))

rp = VlnPlot(scRNA_merge, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3,) + NoLegend() #pt.size = 0.1显示点
pdf("11_raw_samples_features_before_filtration.pdf", height = 12, width = 17)
rp
dev.off()
# 基因数量随细胞数量的变化关系
FeatureScatter(scRNA_merge, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
# 线粒体基因数量随细胞数量的变化关系
FeatureScatter(scRNA_merge, "nCount_RNA", "percent.mt", group.by = "orig.ident", pt.size = 0.5)


###################################
#========过滤1：数据质控===========
###################################

# 对RNA及线粒体数初步过滤
# 对数据结构即subset命令
scRNA_filt <- subset(scRNA_merge,
                         subset = nFeature_RNA > 200& 
                          nFeature_RNA < 5000 &
                          percent.hb <0.1 & 
                          percent.mt < 5 &
                          percent.ribo < 30 &
                          nCount_RNA < 40000)

rp2 = VlnPlot(scRNA_filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3,) + NoLegend() #pt.size = 0.1显示点
pdf("12_raw_samples_features_after_filtration.pdf", height = 12, width = 17)
rp2
dev.off()

table(scRNA_filt@meta.data$orig.ident)
#GSM4664009 GSM4664010 GSM4664011 GSM4664012 GSM4664013 GSM4664015 GSM4664016 GSM4664017 
#       557        191        307        219        780        197        552        165 
#GSM4664019 GSM4664020 
#        87        264 

#添加分组信息
##修改样本名
group= scRNA_filt@meta.data[["orig.ident"]]
group[which(group %in% c("GSM4664009","GSM4664010","GSM4664011","GSM4664012"))]="HC"
group[which(group %in% c("GSM4664013","GSM4664015","GSM4664016","GSM4664017","GSM4664019","GSM4664020"))]="AML"
scRNA_filt$Group = group

scRNA_filt$group = "AML"

###################################
#========过滤2：数据质控===========
###################################

#看高表达异常基因
par(mar = c(4, 8, 2, 1))
C <- scRNA_filt@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
dev.off()  #Rplots.pdf
pdf("1_high_expressedTopGene_after_filtration.pdf", height = 12, width = 12)
print(t)
dev.off()
##基于核糖体与线粒体的过滤指标
#可以从上图看出，在过滤结果中占比较多的基因依旧是与核糖体和线粒体相关的基因，对这些基因进行下一步的过滤的时候，有三种选择可以考虑：
#第一种：如果在第一步基本过滤后，还存在较多的细胞，那么可以把线粒体基因过滤的阈值降到最低，直接过滤到含有线粒体基因的细胞。
#第二种：删除所有与线粒体相关的基因，
#第三种：对每个细胞的线粒体基因含量占比进行计算，删除那些线粒体相关基因含量占比较高的细胞。
#一般采用较低的阈值先过滤异常的细胞，然后在除去线粒体与核糖体的基因，这篇教程里就是这么做的。

#可以看出MALAT1基因的数量占了基因总数的7%，这是一个非常不正常的值，可能是由于实验或测序过程中的技术错误导致的，所以在后续的分析中需要讲这个基因从基因列表中删除，在质控中，这一步很重要，是对质控结果的一次检验，也就是说，质控的结果中，不应该存在一个基因占基因总数较大的百分比。
# Filter MALAT1
scRNA_filt <- scRNA_filt[!grepl("MALAT1", rownames(scRNA_filt)), ]


###################################
###################################
#=========3.数据整合、降维、聚类分群===========
###################################
###################################

# 转换为Seurat v4格式会自动合并层
scRNA_v4 <- as(scRNA_filt, "Seurat")
scRNA_filt = scRNA_v4
scRNA_mfilt <- NormalizeData(scRNA_filt, layer = "counts") %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)

###选择合适的PCs数
plot1 <- DimPlot(scRNA_mfilt, reduction = "pca", group.by="orig.ident") 
###ElbowPlot() 可以快速的检查降维的效果（Q4:ndims需要调整）
plot2 <- ElbowPlot(scRNA_mfilt, ndims=50, reduction="pca") 
###我们一般选择拐点作为降维的度数。
plotc <- plot1+plot2
plotc
pdf("21_pca_elbowplot_pcNu.pdf",height=7,width=14)
plotc 
dev.off()
#后续分析要根据右图选择提取的pc轴数量，一般选择斜率平滑的点之前的所有pc轴，这里选择40
ndims=22

#####去批次效应前
scRNA_mfilt <- FindNeighbors(scRNA_mfilt, dims = 1:ndims) %>% FindClusters(resolution = 0.1) ##FindClusters即编辑ident
table(scRNA_mfilt@meta.data$seurat_clusters)
scRNA_mfilt <- RunUMAP(scRNA_mfilt, dims = 1:ndims)
head(scRNA_mfilt@meta.data)
plota =DimPlot(scRNA_mfilt, reduction = "umap",group.by=c('ident',"orig.ident"),label = T) #聚类分布
plota
#combinate
pdf("22_raw_Nobatch_umap_dimplot.pdf",height=10,width=14)
plota
dev.off()

#####去批次后
# head(scRNA_mfilt@meta.data)
#先用不同来源的做批次效应
#system.time({scRNA_mfilt_harmony <- RunHarmony(scRNA_mfilt, group.by.vars = "Project")})#分组信息只存在于meta.data中 lusc_scRNA_4datasets_filt_harmonyProject.Rdata
system.time({scRNA_mfilt_harmony <- RunHarmony(scRNA_mfilt, group.by.vars = "orig.ident")})

#选择合适的resolution (Q5:关于分辨率的选择)
seq <- seq(0.1, 0.5, by = 0.1)
for(res in seq){
  scRNA_mfilt_harmony <- FindClusters(scRNA_mfilt_harmony, resolution = res)
}
#clustree通过可视化不同分辨率下不同cluster之间的交互关系来帮助我们选择合适的分辨率来进行下游分析。
#remotes::install_version("clustree")
library(clustree)
library(patchwork)
p1 <- clustree(scRNA_mfilt_harmony, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(scRNA_mfilt_harmony, group.by = 'RNA_snn_res.0.1', label = T)
p3 <- DimPlot(scRNA_mfilt_harmony, group.by = 'RNA_snn_res.0.2', label = T)
p4 <- DimPlot(scRNA_mfilt_harmony, group.by = 'RNA_snn_res.0.3', label = T)
p5 <- DimPlot(scRNA_mfilt_harmony, group.by = 'RNA_snn_res.0.4', label = T)
p6 <- DimPlot(scRNA_mfilt_harmony, group.by = 'RNA_snn_res.0.5', label = T)
p1 + p2  +p3 + p4+ p5 + p6+plot_layout(widths = c(3, 2))
ggplot2::ggsave('23_resolution_cuttree.pdf',plot = p1 + p2  +p3 + p4+ p5 + p6+plot_layout(widths = c(2, 1)),he=15,wi=27)


#通过这个图我们会发现当resolution大于0.4后不同cluster在不同的分辨率下会存在越来越多的相互交织，结合文章文群数量，我们选择resolution=0.1来进行后面的分析。
###根据resolution聚类
res=0.1
scRNA_mfilt_harmony <- FindNeighbors(scRNA_mfilt_harmony,reduction = "harmony", dims = 1:ndims) %>% FindClusters(resolution = 0.1)
#细胞亚群数量
table(scRNA_mfilt_harmony@meta.data$seurat_clusters)
#UMAP
scRNA_mfilt_harmony <- RunUMAP(scRNA_mfilt_harmony, reduction = "harmony", dims = 1:ndims)
DimPlot(scRNA_mfilt_harmony, reduction = "umap",group.by=c('ident',"orig.ident"),label = T) 

#看批次效应是否存在 
plot_pca <- DimPlot(scRNA_mfilt_harmony, reduction = "pca", group.by="orig.ident") #pc图
plot1 =DimPlot(scRNA_mfilt_harmony, reduction = "umap",group.by=c('ident',"orig.ident"),label = T) #聚类分布
DimPlot(scRNA_mfilt_harmony, reduction = "umap",group.by="Group",label = T)
plota <- plot_pca+plot1
pdf("24_batchSamples_umap_dimplot.pdf",height=10,width=17)
plota
dev.off()

save(scRNA_mfilt_harmony,file = "1_GSE154109_samples_scRNA_filt_harmony.Rdata")

##################################################################################
#细胞注释
#install.packages("BiocManager")
#BiocManager::install("SingleR")
library(SingleR)
#BiocManager::install("celldex")
library(celldex)
#BiocManager::install("BiocParallel")
library(BiocParallel)
library(Seurat)

#下载
#hpca.se <- HumanPrimaryCellAtlasData()
rm(list=ls())
load("1_GSE154109_samples_scRNA_filt_harmony.Rdata")
load("F:/R_workdir/singleR_databases/HumanPrimaryCellAtlas_hpca.se_human.RData") #之前已经下载好的hpca.se
load("F:/R_workdir/singleR_databases/ref_Hematopoietic.RData") #ref_Hematopoietic

singleRref = ref_Hematopoietic
#获取基因的表达谱的count数据
#singleRdata <- GetAssayData(scRNA_mfilt_harmony, assay = "RNA", data="counts")
singleRdata <- JoinLayers(scRNA_mfilt_harmony, assay = "RNA")  # 先合并layers
singleRdata <- GetAssayData(singleRdata, assay = "RNA")  
#获取聚类的亚群
clusters <- scRNA_mfilt_harmony@meta.data$seurat_clusters
cellpred <- SingleR(test =  singleRdata, 
                    ref = singleRref, 
                    labels = singleRref$label.fine, #注释小类别
					#labels = hpca.se$label.main,注释大类别
                    clusters = clusters, 
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")
					
plotScoreHeatmap(cellpred) #26_singleR_anno_plotScoreHeatmap.pdf
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)

# meta中添加singleR的注释结果
scRNA_mfilt_harmony@meta.data$Celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA_mfilt_harmony@meta.data[which(scRNA_mfilt_harmony@meta.data$seurat_clusters == celltype$ClusterID[i]),'Celltype'] <- celltype$celltype[i]}

#根据singleR的结果找更大的注释
ltype = vector()
t = scRNA_mfilt_harmony@meta.data$Celltype
for (i in 1:length(t)){
	type = strsplit(t[i],split = ":")[[1]][1]
	ltype = c(ltype,type)
}
scRNA_mfilt_harmony@meta.data$Celltype_main = ltype
 
#生成注释图 
plot1 =DimPlot(scRNA_mfilt_harmony, reduction = "umap",group.by='orig.ident',label = F)  
plot2 = DimPlot(scRNA_mfilt_harmony, reduction = "umap", group.by='Celltype',label = T,label.size = 5)
plot3 = DimPlot(scRNA_mfilt_harmony, reduction = "umap", group.by='Celltype_main',label = T,label.size = 5)
 #按照聚类
plota <- plot1+plot2+plot3
pdf("31_cluster_umap_dimplot_singleR_celltype.pdf",height=10,width=17)
plota
dev.off()

save(scRNA_mfilt_harmony,file = "scRNA_merge_celltype_annotated.Rdata")


scRNA=scRNA_mfilt_harmony
scRNA@active.ident=factor(scRNA_mfilt_harmony$Celltype_main,levels = unique(scRNA_mfilt_harmony$Celltype_main))
save(scRNA,file = "3scRNA_merge_annotated.Rdata")


##################################
###################################
########=====4：细胞亚群注释=======
###################################
###################################
rm(list=ls())
load("scRNA1.Rdata")
#setwd("2_updataed_20241101")
scRNA2  <- JoinLayers(scRNA)
markers <- FindAllMarkers(object = scRNA2, test.use = "wilcox",
                          only.pos = TRUE,#只输入表达参数大于0的基因
                          logfc.threshold = 0.25, min.pct=0.25,)  
						  
						  
##markers数据中pct.1的意思是该基因在该cluster中百分之多少的细胞有表达；
##pct.2的意思是该基因在除该cluster之外的百分之多少的细胞有表达
all.markers=markers %>% dplyr::select(gene,everything()) %>% subset(p_val<0.05)
write.csv(all.markers,'4all.markers.csv',row.names = F)
all.markers<-read.csv("4all.markers.csv",header = TRUE)
#rm(list=ls())

top10=all.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
#显示每个类群的前10个markers的热图 
pdf('4_top10.heatmap.pdf', width = 14, height = 20)
DoHeatmap(scRNA2,features=top10$gene) + NoLegend() 
dev.off()

#scRNA1=scRNA2
pdf("4_annotated.clusters.dimplot.pdf",height=7,width=10)
DimPlot(scRNA1, reduction = "umap",label = T) 
dev.off()

scRNA1$celltype=scRNA1@active.ident #也是scRNA1@meta.data$celltype

pdf("4_annotated.clusters_bygroup.dimplot.pdf",height=7,width=14)
DimPlot(scRNA1,split.by = 'Group')
dev.off()

save(scRNA1,file = "scRNA1_annotated.Rdata")		





####################################
###################################
#=========细胞亚群细分 ===========
###################################
###################################
setwd("F:\\helping\\#huangjl\\scRNA\\AML_GSE154109_Rout")
load("scRNA1_annotated.Rdata")

outdir="2_AML_mono_IL4I"
if(!dir.exists(outdir)) dir.create(outdir)
setwd(outdir)

scRNA1@active.ident=factor(scRNA1$celltype,levels = unique(scRNA1$celltype))
scRNA1_mono=subset(scRNA1,ident=c("Monocytes")) #Monocytes = 1485
GeneExpNor <-  GetAssayData(object = scRNA1_mono, slot = 'data') #
IL4I1_exp = GeneExpNor["IL4I1",]
IL4I1_group <- IL4I1_exp
IL4I1_group[IL4I1_exp==0] =  "IL4I1-"
IL4I1_group[IL4I1_exp>0] =  "IL4I1+"

scRNA1_mono$IL4I1_group = IL4I1_group

scRNA1_mono <- NormalizeData(scRNA1_mono) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)
  
###选择合适的PCs数
plot1 <- DimPlot(scRNA1_mono, reduction = "pca", group.by="orig.ident") 
###ElbowPlot() 可以快速的检查降维的效果（Q4:ndims需要调整）
plot2 <- ElbowPlot(scRNA1_mono, ndims=50, reduction="pca") 
###我们一般选择拐点作为降维的度数。
plotc <- plot1+plot2
plotc 
pdf("51_pca_elbowplot_scRNA1_mono_before.pdf",height=7,width=14)
plotc 
dev.off()
#后续分析要根据右图选择提取的pc轴数量，一般选择斜率平滑的点之前的所有pc轴，这里选择40
ndims=14

#?? Provided graph.name not present in Seurat object
#选择合适的resolution (Q5:关于分辨率的选择)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  scRNA1_mono <- FindClusters(scRNA1_mono, resolution = res)
}
#clustree通过可视化不同分辨率下不同cluster之间的交互关系来帮助我们选择合适的分辨率来进行下游分析。
#remotes::install_version("clustree")
library(clustree)
library(patchwork)
p1 <- clustree(scRNA1_mono, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.1', label = T)
p3 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.2', label = T)
p4 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.3', label = T)
p5 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.4', label = T)
p6 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.5', label = T)
p7 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.1', label = T)
p1 + p2  +p3 + p4+ p5 + p6 + p7 + plot_layout(widths = c(6, 4))
ggplot2::ggsave('51_resolution_cuttree_scRNA1_mono.pdf',plot = p1 + p2  +p3 + p4+ p5 + p6 + p7 +plot_layout(widths = c(3, 1)),he=15,wi=15)
p5
resolution = 0.1  # 8   0.7:9

#
scRNA1_mono <- FindNeighbors(scRNA1_mono, dims = 1:ndims) %>% FindClusters(resolution = 0.1) 

table(scRNA1_mono@meta.data$seurat_clusters)
scRNA1_mono <- RunUMAP(scRNA1_mono, dims = 1:ndims)
head(scRNA1_mono@meta.data)
#plot1 = DimPlot(scRNA1_mono, reduction = "umap", group.by='group') 
#plot2 = DimPlot(scRNA1_mono, reduction = "umap", group.by='sample')#样本分布
plot3 = DimPlot(scRNA1_mono, reduction = "umap",group.by='IL4I1_group',label = T) 
plot4 = DimPlot(scRNA1_mono, reduction = "umap",group.by='ident',label = T) #聚类分布

#combinate
plota <- plot3+plot4
plota
pdf("52_Monophage_Nobatch_umap_dimplot.pdf",height=7,width=14)
plota
dev.off()

save(scRNA1_mono, file="scRNA1_mono_raw.RData")

#===========指定分组IL4I1计算差异基因====
load("scRNA1_mono.RData")
Idents(scRNA1_mono) <- "IL4I1_group"
head(scRNA1_mono@active.ident)

#差异和功能分析
CELLDEG <- FindMarkers(scRNA1_mono, ident.1 = "IL4I1+", ident.2 = "IL4I1-", verbose = FALSE)
my_enrichment(CELLDEG,"IL4I1+-_Mono.DEGs","hsa") #见尾部自定义函数

 
#计算平均表达量，输出结果
#提取标准化后的表达矩阵,包括分组信息
GeneExpNor <-  GetAssayData(object = scRNA1_mono, slot = 'data') #等价于GeneExpNor <- as.data.frame(scRNA1_mono@assays$RNA@data)#行为基因，列为细胞
GeneExpNor<-t(GeneExpNor)%>%as.data.frame() #行为细胞，列为基因
group_info <- scRNA1_mono@meta.data[,c("IL4I1_group","celltype")]
GeneExpdat <- cbind(group_info,GeneExpNor)
rownames(GeneExpdat) <- rownames(GeneExpNor) #前三列是分组信息，其余列才是表达数据
#计算平均值		
case_expdat <- GeneExpdat[GeneExpdat$IL4I1_group == "IL4I1+",3:ncol(GeneExpdat)]
control_expdat <- GeneExpdat[GeneExpdat$IL4I1_group == "IL4I1-",3:ncol(GeneExpdat)]
case_avgexp <- apply(case_expdat, 2, mean)
control_avgexp <- apply(control_expdat, 2, mean)
degs_id <- rownames(CELLDEG)
CELLDEG_exp_out <- data.frame(cbind(Gene=rownames(CELLDEG),IL4I1pos_case_xpression=case_avgexp[degs_id],IL4I1neg_control_xpression=control_avgexp[degs_id],CELLDEG))
#输出差异分析结果和表达
write.csv(CELLDEG_exp_out,file="IL4I1+vsIL4I1-_Mono_DEGsa_all.CSV",row.names=F)

#显著差异
CELLDEG_exp_sig <- CELLDEG_exp_out[CELLDEG_exp_out$p_val_adj<=0.05,]
write.csv(CELLDEG_exp_sig,file="IL4I1+vsIL4I1-_Mono_DEG_adjP0.05.CSV",row.names=F)




####################################
###################################
#=========细胞亚群细分 ===========
###################################
###################################
load("scRNA1_annotated.Rdata")

outdir="2_AML_mono_IL4I"
if(!dir.exists(outdir)) dir.create(outdir)
setwd(outdir)

scRNA1@active.ident=factor(scRNA1$celltype,levels = unique(scRNA1$celltype))
scRNA1_mono=subset(scRNA1,ident=c("Monocytes")) #Monocytes = 1485
GeneExpNor <-  GetAssayData(object = scRNA1_mono, slot = 'data') #
IL4I1_exp = GeneExpNor["IL4I1",]
IL4I1_group <- IL4I1_exp
IL4I1_group[IL4I1_exp==0] =  "IL4I1-"
IL4I1_group[IL4I1_exp>0] =  "IL4I1+"

scRNA1_mono$IL4I1_group = IL4I1_group

scRNA1_mono <- NormalizeData(scRNA1_mono) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(verbose=FALSE)
  
###选择合适的PCs数
plot1 <- DimPlot(scRNA1_mono, reduction = "pca", group.by="orig.ident") 
###ElbowPlot() 可以快速的检查降维的效果（Q4:ndims需要调整）
plot2 <- ElbowPlot(scRNA1_mono, ndims=50, reduction="pca") 
###我们一般选择拐点作为降维的度数。
plotc <- plot1+plot2
plotc 
pdf("51_pca_elbowplot_scRNA1_mono_before.pdf",height=7,width=14)
plotc 
dev.off()
#后续分析要根据右图选择提取的pc轴数量，一般选择斜率平滑的点之前的所有pc轴，这里选择40
ndims=14

#?? Provided graph.name not present in Seurat object
#选择合适的resolution (Q5:关于分辨率的选择)
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  scRNA1_mono <- FindClusters(scRNA1_mono, resolution = res)
}
#clustree通过可视化不同分辨率下不同cluster之间的交互关系来帮助我们选择合适的分辨率来进行下游分析。
#remotes::install_version("clustree")
library(clustree)
library(patchwork)
p1 <- clustree(scRNA1_mono, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.1', label = T)
p3 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.2', label = T)
p4 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.3', label = T)
p5 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.4', label = T)
p6 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.0.5', label = T)
p7 <- DimPlot(scRNA1_mono, group.by = 'RNA_snn_res.1', label = T)
p1 + p2  +p3 + p4+ p5 + p6 + p7 + plot_layout(widths = c(6, 4))
ggplot2::ggsave('51_resolution_cuttree_scRNA1_mono.pdf',plot = p1 + p2  +p3 + p4+ p5 + p6 + p7 +plot_layout(widths = c(3, 1)),he=15,wi=15)
p5
resolution = 0.1  # 8   0.7:9

#
scRNA1_mono <- FindNeighbors(scRNA1_mono, dims = 1:ndims) %>% FindClusters(resolution = 0.1) 

table(scRNA1_mono@meta.data$seurat_clusters)
scRNA1_mono <- RunUMAP(scRNA1_mono, dims = 1:ndims)
head(scRNA1_mono@meta.data)
#plot1 = DimPlot(scRNA1_mono, reduction = "umap", group.by='group') 
#plot2 = DimPlot(scRNA1_mono, reduction = "umap", group.by='sample')#样本分布
plot3 = DimPlot(scRNA1_mono, reduction = "umap",group.by='IL4I1_group',label = T) 
plot4 = DimPlot(scRNA1_mono, reduction = "umap",group.by='ident',label = T) #聚类分布

#combinate
plota <- plot3+plot4
plota
pdf("52_Monophage_Nobatch_umap_dimplot.pdf",height=7,width=14)
plota
dev.off()

save(scRNA1_mono, file="scRNA1_mono_raw.RData")

#===========指定分组IL4I1计算差异基因====
load("scRNA1_mono.RData")
Idents(scRNA1_mono) <- "IL4I1_group"
head(scRNA1_mono@active.ident)

#差异和功能分析
CELLDEG <- FindMarkers(scRNA1_mono, ident.1 = "IL4I1+", ident.2 = "IL4I1-", verbose = FALSE)
my_enrichment(CELLDEG,"IL4I1+-_Mono.DEGs","hsa") #见尾部自定义函数

 
#计算平均表达量，输出结果
#提取标准化后的表达矩阵,包括分组信息
GeneExpNor <-  GetAssayData(object = scRNA1_mono, slot = 'data') #等价于GeneExpNor <- as.data.frame(scRNA1_mono@assays$RNA@data)#行为基因，列为细胞
GeneExpNor<-t(GeneExpNor)%>%as.data.frame() #行为细胞，列为基因
group_info <- scRNA1_mono@meta.data[,c("IL4I1_group","celltype")]
GeneExpdat <- cbind(group_info,GeneExpNor)
rownames(GeneExpdat) <- rownames(GeneExpNor) #前三列是分组信息，其余列才是表达数据
#计算平均值		
case_expdat <- GeneExpdat[GeneExpdat$IL4I1_group == "IL4I1+",3:ncol(GeneExpdat)]
control_expdat <- GeneExpdat[GeneExpdat$IL4I1_group == "IL4I1-",3:ncol(GeneExpdat)]
case_avgexp <- apply(case_expdat, 2, mean)
control_avgexp <- apply(control_expdat, 2, mean)
degs_id <- rownames(CELLDEG)
CELLDEG_exp_out <- data.frame(cbind(Gene=rownames(CELLDEG),IL4I1pos_case_xpression=case_avgexp[degs_id],IL4I1neg_control_xpression=control_avgexp[degs_id],CELLDEG))
#输出差异分析结果和表达
write.csv(CELLDEG_exp_out,file="IL4I1+vsIL4I1-_Mono_DEGsa_all.CSV",row.names=F)

#显著差异
CELLDEG_exp_sig <- CELLDEG_exp_out[CELLDEG_exp_out$p_val_adj<=0.05,]
write.csv(CELLDEG_exp_sig,file="IL4I1+vsIL4I1-_Mono_DEG_adjP0.05.CSV",row.names=F)
