#加载相关R包
library(Seurat)
library(dplyr)
library(ggplot2)
library(plotly)
library(plot3D)

#设置工作目录
#读取之前已经完成分析的单细胞转录组数据
load('ST_anno.Rda')
dim(st)

#1.个性化UMAP图绘制
#从seurat对象中提取坐标信息
st<-RunUMAP(st,dims = 1:30,label = T)
UMAP.coor<-st@reductions$umap@cell.embeddings
UMAP.coor<-as.data.frame(UMAP.coor)

#从seurat对象中提取细胞注释信息
celltype<-as.factor(st$celltype)
UMAP.coor<-cbind(UMAP.coor,celltype)

#从seurat对象中提取基因表达量信息
Col1a2<-FetchData(st,vars = "Col1a2")
UMAP.coor<-cbind(UMAP.coor,Col1a2)

#从seurat对象中提取Oligo细胞marker基因平均表达量
oligo.marker <- FetchData(st,vars = c('Olig1','Olig2','Mbp'))
oligo.marker<-expm1(oligo.marker)
oligo.marker<-log1p(rowMeans(oligo.marker))
UMAP.coor<-cbind(UMAP.coor,oligo.marker)

#使用ggplot2来完成个性化UMAP图的绘制
p1<-ggplot(data = UMAP.coor,aes(x = umap_1,y = umap_2))+
  geom_point(aes(color = celltype))
p1
#在图中增加标签信息
p2<-ggplot(data = UMAP.coor,aes(x = umap_1,y = umap_2,label = celltype))+
  geom_point(aes(color = celltype))+
  geom_text(color = 'black',size = 3)
p2

label.coor<-UMAP.coor%>%group_by(celltype)%>%summarise(umap_1=median(umap_1),umap_2=median(umap_2))
p2<-ggplot(data = UMAP.coor,aes(x = umap_1,y = umap_2,label = celltype))+
  geom_point(aes(color = celltype))+
  geom_text(data = label.coor,color = 'black',size = 5)
p2

#自定义颜色
color<-c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
         "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
         "#756455")
p3<-p2+scale_color_manual(values = color)
p3

#自定义主题
mytheme<-theme_classic()+
  theme(panel.background = element_rect(fill = 'white',color = 'white'),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(color = 'black',size = 10),
        axis.title = element_text(color = 'black',size = 12),
        plot.title = element_text(size = 18,hjust = 0.5)
        )
p4<-p3+ggtitle('Cell Type')+mytheme
p4

#2.3D UMAP图绘制
rm(list = ls())
load('ST_anno.Rda')
#利用seurat进行三维降维，并从seurat对象中提取坐标信息
st<-RunUMAP(st,dims = 1:30,label = T,n.components = 3L,reduction.name = 'umap',reduction.key = 'umap_')
UMAP3D.coor<-st@reductions$umap@cell.embeddings
UMAP3D.coor<-as.data.frame(UMAP3D.coor)

#从seurat对象中提取细胞注释信息
celltype<-as.factor(st$celltype)
UMAP3D.coor<-cbind(UMAP3D.coor,celltype)

#从seurat对象中提取基因表达量信息
Col1a2<-FetchData(st,vars = 'Col1a2')
UMAP3D.coor<-cbind(UMAP3D.coor,Col1a2)

#从seurat对象中提取Oligo细胞marker基因平均表达量
oligo.marker <- FetchData(st,vars = c('Olig1','Olig2','Mbp'))
oligo.marker<-expm1(oligo.marker)
oligo.marker<-log1p(rowMeans(oligo.marker))
UMAP3D.coor<-cbind(UMAP3D.coor,oligo.marker)

#利用这个表格绘制图形
#方法一：使用plotly包进行绘制
p1<-plot_ly(UMAP3D.coor,
            x = ~umap_1,y = ~umap_2,z = ~umap_3,
            type = 'scatter3d',
            mode = 'markers',
            marker = list(size=2))
p1
p2<-plot_ly(UMAP3D.coor,
            x = ~umap_1,y = ~umap_2,z = ~umap_3,
            type = 'scatter3d',
            mode = 'markers',
            marker = list(size=2),
            color = celltype,
            text = celltype,
            hoverinfo = 'text')
p2
#绘制某一个基因的表达量分布图
p3<-plot_ly(UMAP3D.coor,
            x = ~umap_1,y = ~umap_2,z = ~umap_3,
            type = 'scatter3d',
            mode = 'markers',
            marker = list(size=2,width=1),
            color = ~Col1a2,
            text = celltype,
            hoverinfo = 'text')
p3
#绘制Oligo细胞marker基因平均表达量分布图
p4<-plot_ly(UMAP3D.coor,
            x=~umap_1,y=~umap_2,z=~umap_3,
            type = 'scatter3d',
            mode='markers',
            marker=list(size=2,width=1),
            color = ~oligo.marker,
            text=celltype,
            hoverinfo='text')
p4

#方法二：使用plot3D进行绘制
UMAP3D.coor$celltype1<-as.numeric(UMAP3D.coor$celltype)
scatter3D(x = UMAP3D.coor$umap_1,y = UMAP3D.coor$umap_2,z = UMAP3D.coor$umap_3,
          colvar = UMAP3D.coor$celltype1)
with(UMAP3D.coor,scatter3D(x = umap_1,y = umap_2,z = umap_3,colvar = celltype1))
#自定义颜色
color<-c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
         "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
         "#756455")
#修改图形颜色、坐标轴名称、图例、角度调整
with(UMAP3D.coor,scatter3D(x = umap_1,y = umap_2,z = umap_3,
                         colvar = celltype1,col = color,colkey = F,
                         pch = 20,bty = 'b2',cex = 0.8,
                         theta = 25,phi = 35,
                         xlab = 'UMAP_1',ylab = 'UMAP_2',zlab = 'UMAP_3'))
legend('right',title = 'celltype',legend = unique(UMAP3D.coor$celltype),
       pch = 20,col = color,bty = 'n',title.adj = 0)

#绘制某一个基因的表达量分布图
color1<-colorRampPalette(c("#4476B6", "#FFFEC2", "#FD9B4A"))(50)
with(UMAP3D.coor,scatter3D(x = umap_1,y = umap_2,z = umap_3,
                         colvar = Col1a2,col = color1,clab = 'Col1a2',#colkey = F,
                         pch = 20,bty = 'b2',cex = 0.8,
                         theta = 25,phi = 35,
                         xlab = 'UMAP_1',ylab = 'UMAP_2',zlab = 'UMAP_3'))

#绘制oligo细胞marker基因平均表达量分布图
with(UMAP3D.coor,scatter3D(x = umap_1,y = umap_2,z = umap_3,
                         colvar = oligo.marker,col = color1,clab = 'Oligo',
                         pch = 20,bty = 'b2',cex = 0.8,
                         theta = 40,phi = 50,
                         xlab = 'UMAP_1',ylab = 'UMAP_2',zlab = 'UMAP_3'))

#保存图片
pdf('3DUMAP.pdf')
with(UMAP3D.coor,scatter3D(x = umap_1,y = umap_2,z = umap_3,
                           colvar = oligo.marker,col = color1,clab = 'Oligo',
                           pch = 20,bty = 'b2',cex = 0.8,
                           theta = 25,phi = 35,
                           xlab = 'UMAP_1',ylab = 'UMAP_2',zlab = 'UMAP_3'))
dev.off()
