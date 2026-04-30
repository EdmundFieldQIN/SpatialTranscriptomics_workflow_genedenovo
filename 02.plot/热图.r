#加载相关R包
library(Seurat)
library(pheatmap)
library(ggplot2)
library(dplyr)

#读取单细胞转录组数据
load('ST_anno.Rda')
dim(st)
##将默认细胞注释由seurat_cluster切换为celltype
st@active.ident<-as.factor(st$celltype)

#基于上调基因分析挑选用于绘图的基因
dif<-FindAllMarkers(st,logfc.threshold = 0.6,min.pct = 0.4,only.pos = T)
sig.dif<-dif%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
genes<-unique(sig.dif$gene)

#获取绘图使用的表达量信息
st<-ScaleData(st,features = row.names(st))
data<-st@assays$Spatial$scale.data[genes,]

#默认参数绘图
pheatmap(data)

#去除行列聚类、归一化和细胞名称展示，简化图形
pheatmap(data,scale = "none",cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE)

#将细胞按照细胞类型进行排布,并增加细胞注释和上调基因所属细胞类型注释
celltype<-st$celltype #调用细胞注释信息
celltype<-celltype[order(celltype)] #将细胞按照细胞注释进行排序
celltype<-data.frame(celltype)
data<-data[,rownames(celltype)]
sig.dif<-sig.dif[!duplicated(sig.dif$gene),] #保留不重复的基因的信息
gene.anno<-data.frame(gene.anno=sig.dif$cluster,row.names = sig.dif$gene)
pheatmap(data,scale = "none",cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,
         annotation_col = celltype,annotation_row = gene.anno,annotation_names_row = FALSE)

#因为单细胞数据表达量差异较大，我们可以将高于2.5的值设置为同一个颜色，将低于-2.5的值设置为同一个颜色
data[which(data > 2.5)] <- 2.5
data[which(data < -2.5)] <- -2.5
pheatmap(data,scale = "none",cluster_cols = FALSE,cluster_rows = FALSE,show_colnames = FALSE,
         annotation_col = celltype,annotation_row = gene.anno,annotation_names_row = FALSE,
         legend_breaks = seq(-2,2,1),legend_labels = seq(-2,2,1))

#更改渐变颜色和图例颜色
color<-colorRampPalette(c("navy", "white", "firebrick3"))(50)
color1<-c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
         "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
         "#756455")
celltype.col<-color1[1:length(unique(celltype$celltype))]
names(celltype.col)<-unique(celltype$celltype)
col<-list(celltype = celltype.col,gene.anno = celltype.col)
pheatmap(data,scale = "none",
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,
         annotation_col = celltype,annotation_row = gene.anno,annotation_names_row = FALSE,
         legend_breaks = seq(-2,2,1),legend_labels = seq(-2,2,1),
         annotation_colors =  col,color = color)

#为细胞增加样本的注释信息
sample<-st$orig.ident
celltype<-st$celltype
meta<-data.frame(celltype,sample)
sample.col<-c(anterior = "#9ab7d1",posterior = "#efb9a0")
col<-list(celltype = celltype.col,gene.anno = celltype.col,sample = sample.col)
pheatmap(data,scale = "none",
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,
         annotation_col = meta,annotation_row = gene.anno,annotation_names_row = FALSE,
         legend_breaks = seq(-2,2,1),legend_labels = seq(-2,2,1),
         annotation_colors =  col,color = color)

#按照细胞类型切割热图
gap<-c()
stat<-0
celltype<-celltype[order(celltype)]
celltype<-unique(celltype)
for(i in 1:length(celltype)){
  num<-sum(meta$celltype==celltype[i])
  stat<-stat+num
  gap[i]<-stat
}
pheatmap(data,scale = "none",
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,
         annotation_col = meta,annotation_row = gene.anno,annotation_names_row = FALSE,
         legend_breaks = seq(-2,2,1),legend_labels = seq(-2,2,1),
         annotation_colors =  col,color = color,gaps_col = gap)

#更换细胞排布方式，以样本注释为基础注释
#为了保证细胞类型内部细胞按照细胞类型进行排序，先对meta进行拆分
meta.list<-lapply(unique(meta$sample),FUN = function(x){
  a<-meta[which(sample==x),]
  a<-a[order(a$celltype),]
  return(a)
})
meta<-do.call('rbind',meta.list)
meta<-meta[,c(2,1)]
data<-data[,rownames(meta)]
pheatmap(data,scale = "none",
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,
         annotation_col = meta,annotation_row = gene.anno,annotation_names_row = FALSE,
         legend_breaks = seq(-2,2,1),legend_labels = seq(-2,2,1),
         annotation_colors =  col,color = color)

#保存图片
pdf("heatmap.pdf")
pheatmap(data,scale = "none",
         cluster_rows = FALSE,cluster_cols = FALSE,show_colnames = FALSE,
         annotation_col = meta,annotation_row = gene.anno,annotation_names_row = FALSE,
         legend_breaks = seq(-2,2,1),legend_labels = seq(-2,2,1),
         annotation_colors =  col,color = color)
dev.off()
