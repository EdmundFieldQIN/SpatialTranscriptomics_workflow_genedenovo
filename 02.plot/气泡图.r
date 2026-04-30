#加载相关R包
BiocManager::install("MAST")
library(Seurat)
library(ggplot2)
library(dplyr)
library(MAST)
#读取单细胞转录组数据
load('ST_anno.Rda')
dim(st)
##将默认细胞注释由seurat_cluster切换为celltype
st@active.ident<-as.factor(st$celltype)

#基于上调基因分析挑选用于绘图的基因 大约5min
dif<-FindAllMarkers(st,logfc.threshold = 0.6,min.pct = 0.4,only.pos = T,test.use = 'MAST')
#若时间太长可直接载入dif
dif<-readRDS('dif.rds')
sig.dif<-dif%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
genes<-unique(sig.dif$gene)

#根据基因筛选出所有数据
data<-t(FetchData(st,vars = genes))
data<-expm1(data)
#计算每个cluster每个基因的平均表达量和表达基因的细胞比例
#按照细胞类型拆分数据
data.list<-lapply(unique(st$celltype),FUN = function(x){
  cellid<-colnames(st)[which(st$celltype==x)]
  data.cluster<-data[,cellid]
  return(data.cluster)
})
names(data.list)<-unique(st$celltype)
#计算均值和表达细胞比例
data.sum<-list()
for(i in 1:length(data.list)){
  means<-log1p(rowMeans(data.list[[i]])) #计算均值
  pro<-apply(data.list[[i]],1,FUN = function(x){
    pro<-sum(x>0.5)/length(x) #表达量大于0.5的细胞视为表达基因
    return(pro)
  })
  gene<-factor(rownames(data.list[[i]]),levels = genes)  #将基因设置为因子型向量，使后续绘图顺序固定
  celltype<-rep(names(data.list)[i],nrow(data.list[[i]]))
  stat<-data.frame(gene,celltype,means,pro)
  data.sum[[i]]<-stat
}
dot.data<-do.call("rbind",data.sum)
#数据做归一化处理
data.scale<-lapply(genes,FUN = function(x){
  data1<-dot.data[which(dot.data$gene==x),]
  data1$scale<-as.vector(scale(data1$means))
  return(data1)
})
dot.data<-do.call("rbind",data.scale)
#通过细胞类型和基因确定xy轴关系，绘制初步的散点图
p1<-ggplot(dot.data,aes(x = celltype,y = gene))+geom_point()
p1

#建立散点大小和表达细胞比例和表达量之间的关系
p2<-ggplot(dot.data,aes(x = celltype,y = gene))+
  geom_point(aes(size = pro,color = scale),shape = 16)
p2

#更改颜色
p3<-p2+scale_color_gradient(high = "#E9814E",low = "#71B5B4")
p3

#设置坐标轴、图例、图表的标题
p4<-p3+labs(x = "Cell Type",y = "Genes",size = "Proportion",color = "Expression",title = "Top 5 of DE")
p4

#自定义图表主题
mytheme<-theme_bw()+theme(plot.title = element_text(size = 18,hjust = 0.5,face = "bold"),
                          axis.text.x = element_text(size = 12,color = "black"),
                          axis.text.y = element_text(size = 12,color = "black"),
                          axis.title = element_text(size = 16,color = "black",face = "bold"),
                          legend.text = element_text(size = 12,color = "black"),
                          legend.title = element_text(size = 16,color = "black"))
p5<-p4+mytheme
p5

#可以更改横纵坐标
mytheme1<-theme_bw()+theme(plot.title = element_text(size = 18,hjust = 0.5,face = "bold"),
                          axis.text.x = element_text(size = 12,color = "black",angle = 45,vjust = 1,hjust = 1),
                          axis.text.y = element_text(size = 12,color = "black"),
                          axis.title = element_text(size = 16,color = "black",face = "bold"),
                          legend.text = element_text(size = 12,color = "black"),
                          legend.title = element_text(size = 16,color = "black"))
p6<-p4+coord_flip()+mytheme1
p6

#保存图片
ggsave(filename = "dotplot.pdf",plot = p5)
