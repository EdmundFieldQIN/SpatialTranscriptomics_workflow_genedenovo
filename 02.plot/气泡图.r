#加载相关R包
BiocManager::install("MAST")
library(Seurat)
library(ggplot2)
library(dplyr)
library(MAST)
#读取单细胞转录组数据
getwd()
setwd("/data/h007/workspace/SpatialTranscriptomics/2.Demo/02.plot/")
project.name <- "ST"
output_dir <- "./output/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

load('ST_anno.rda')
dim(st)
##将默认细胞注释由seurat_cluster切换为celltype
st@active.ident<-as.factor(st$celltype)

#基于上调基因分析挑选用于绘图的基因 大约5min
dif<-FindAllMarkers(st,logfc.threshold = 0.6,min.pct = 0.4,only.pos = T,test.use = 'MAST')
#若时间太长可直接载入dif
saveRDS(dif,file = 'dif.rds')
dif <- readRDS('dif.rds')
sig.dif <- dif |> group_by(cluster) |> top_n(n=5,wt=avg_log2FC) |> mutate(diff.pct = pct.1 - pct.2)
genes<-unique(sig.dif$gene)

#根据基因筛选出所有数据
data <- t(FetchData(st,vars = genes))
data <- expm1(data) # 去log化
#计算每个cluster每个基因的平均表达量和表达基因的细胞比例
#按照细胞类型拆分数据
data.list<-lapply(unique(st$celltype),FUN = function(x){
  cellid<-colnames(st)[which(st$celltype==x)] # 提取每个cluster的细胞id
  data.cluster<-data[,cellid] # 提取每个cluster的表达矩阵
  return(data.cluster)
})

names(data.list)<-unique(st$celltype)
# 计算均值和表达细胞比例
data.sum<-list()
for(i in 1:length(data.list)){
  means <- log1p(rowMeans(data.list[[i]])) # 对当前矩阵按行（基因）求平均表达量，得到每个基因的平均表达值。
  
  pro<-apply(data.list[[i]],1,FUN = function(x){
    pro<-sum(x>0.5)/length(x) # 表达量大于0.5的细胞视为表达基因,计算统计当前基因在多少个细胞中表达值大于 0.5（即“表达”的细胞数），除以总细胞数 length(x)，得到表达该基因的细胞比例。
    return(pro)
  })

  gene<-factor(rownames(data.list[[i]]),levels = genes)  # 将基因设置为因子型向量，使后续绘图顺序固定
  celltype<-rep(names(data.list)[i],nrow(data.list[[i]])) # 将细胞类型名重复 nrow 
  stat<-data.frame(gene,celltype,means,pro) # 构建数据框
  data.sum[[i]]<-stat # 将每个细胞类型的统计结果存储在列表中
}
dot.data<-do.call("rbind",data.sum) # 将多个细胞类型的结果合并为一个数据


#数据做归一化处理
data.scale<-lapply(genes,FUN = function(x){
  data1<-dot.data[which(dot.data$gene==x),] # 提取当前基因的数据
  data1$scale<-as.vector(scale(data1$means)) # 对当前基因的平均表达量进行标准化处理，得到一个新的列 scale，表示标准化后的表达水平。scale 函数会将数据中心化（减去均值）并缩放（除以标准差），使得不同基因的表达水平具有可比性。
  return(data1)
})
dot.data<-do.call("rbind",data.scale)


#通过细胞类型和基因确定xy轴关系，绘制初步的散点图
#建立散点大小和表达细胞比例和表达量之间的关系
p <- ggplot(dot.data,aes(x = celltype,y = gene)) +
  geom_point(aes(size = pro,color = scale),shape = 16)  +
  scale_color_gradient(high = "#E9814E",low = "#71B5B4") + # 更改颜色
  labs(x = "Cell Type",y = "Genes",size = "Proportion",color = "Expression",title = "Top 5 of DE") # 设置坐标轴、图例、图表的标题
p

#自定义图表主题
mytheme<-theme_bw() + 
  theme(plot.title = element_text(size = 18,hjust = 0.5,face = "bold"),
                          axis.text.x = element_text(size = 12,color = "black"),
                          axis.text.y = element_text(size = 12,color = "black"),
                          axis.title = element_text(size = 16,color = "black",face = "bold"),
                          legend.text = element_text(size = 12,color = "black"),
                          legend.title = element_text(size = 16,color = "black"))
p5 <- p + mytheme
p5

#可以更改横纵坐标
mytheme1 <- theme_bw() + 
  theme(plot.title = element_text(size = 18,hjust = 0.5,face = "bold"),
                          axis.text.x = element_text(size = 12,color = "black",angle = 45,vjust = 1,hjust = 1),
                          axis.text.y = element_text(size = 12,color = "black"),
                          axis.title = element_text(size = 16,color = "black",face = "bold"),
                          legend.text = element_text(size = 12,color = "black"),
                          legend.title = element_text(size = 16,color = "black"))
p6 <- p + coord_flip() + mytheme1
p6

#保存图片
pdf(paste0(output_dir, project.name, "_dotplot.pdf"), width = 8, height = 6)
p5
dev.off()

pdf(paste0(output_dir, project.name, "_dotplot_flip.pdf"), width = 8, height = 6)
p6
dev.off()