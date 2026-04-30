#安装相关R包
# options(repos=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/",BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor"),timeout=1800)
# install.packages("Seurat")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("hdf5r")
# install.packages("future")
# install.packages("patchwork")
# setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
# install.packages(c("BPCells", "presto", "glmGamPoi"))
# BiocManager::install("glmGamPoi")
# remotes::install_github("bnprks/BPCells/r", mirror = "https://ghfast.top/")
# remotes::install_github("immunogenomics/presto", mirror = "https://ghfast.top/")

#加载相关R包
library(Seurat)
library(BPCells)
library(presto)
library(glmGamPoi)
library(dplyr)
library(ggplot2)
library(hdf5r)
library(patchwork)
library(future)


# plan("multicore", workers = 12)
plan("multisession", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2)  # 10 GiB
#plan("sequential")
#########################
#转换工作目录!切记！！！#
#########################
dir()
setwd("/data/h007/workspace/SpatialTranscriptomics/2.Demo/01.seurat/")
project.name <- "ST"
output_dir <- "./output/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}



#### 1. 数据导入及创建seurat对象
#读取10X ST数据(hdf5)
anterior1 <- Load10X_Spatial(data.dir = "./anterior",slice = "anterior")
posterior1 <- Load10X_Spatial(data.dir = "./posterior/",slice = "posterior")
dim(anterior1)
dim(posterior1)
# 读取ST的三元组文件
# 构建10X ST数据的读取函数，输入参数包括三元组文件路径、切片图像路径、assay名称、切片名称和项目名称
read10X_Spatial <- function (data.dir="",image.dir="",assay = "Spatial", slice = "slice1",project="CreateSeuratObject")
  {data <- Read10X(data.dir = data.dir)  #读取三元组文件
  object <- CreateSeuratObject(counts = data, assay = assay,project = project)  #创建seurat对象
  image <- Read10X_Image(image.dir = image.dir)  #读取切片图形
  image <- image[Cells(x = object)] # 提取切片图形中与seurat对象中spot对应的部分 
  DefaultAssay(object = image) <- assay # 将切片图形的默认assay设置为与seurat对象相同
  object[[slice]] <- image # 将切片图形添加到seurat对象中，切片名称为输入参数slice
  return(object)
}

anterior <- read10X_Spatial(data.dir = "./anterior/filtered_feature_bc_matrix/",
                            image.dir = "./anterior/spatial/",slice = "anterior",project = "anterior")
posterior <- read10X_Spatial(data.dir = "./posterior/filtered_feature_bc_matrix/",
                             image.dir = "./posterior/spatial/",slice = "posterior",project = "posterior")

#合并样品
st <- merge(anterior,posterior,add.cell.ids = c("ante","post"))
#多样品合并
#st <- merge(a,y=c(b,c,d),add.cell.ids = c("a","b","c","d"),project = "four_merged")
rm(anterior1,posterior1)
dim(st)
# [1] 32285  6050
# 32285个基因，一共6050个spot


#### 2. 数据标准化
#对数据进行归一化处理

# #方法一：log均一化，LogNormalize 的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 ) 
# st <- NormalizeData(st,normalization = "LogNormalize")
# #鉴定细胞间表达量高变的基因（feature selection），用于下游分析，PCA
# #这一步的目的是鉴定出spot与spot之间表达量相差很大的基因，用于后续鉴定细胞类型，
# #我们使用默认参数，即“vst”方法选取2000个高变基因。
# st <- FindVariableFeatures(st,selection.method = "vst",verbose = F)
# head(VariableFeatures(st), 10) #查看前10个高变基因
# #PCA分析数据准备，使用ScaleData()进行数据归一化；使用全部的基因进行归一化（5min）
# st <- ScaleData(st,features = rownames(st))
# #线性降维（PCA）,默认用高变基因集，但也可通过 features 参数自己指定
# #同时也可对主成分数量进行定义 
# st <- RunPCA(st,npcs = 50)
# #确定数据集的分群个数 
# #Jackstraw 置换检验算法；重复取样（原数据的 1%），重跑PCA,鉴定p-value较小的PC；计算‘null distribution’(即零假设成立时)时的基因 scores;
# #ST数据背景噪音较小，可以检验更多的主成分以提高分群准确性
# st <- JackStraw(st,dims = 50) #5min
# st <- ScoreJackStraw(st,dims = 1:50)

# pjs <- JackStrawPlot(st,dims = 1:50)
# pdf(paste0(output_dir,project.name,"_JackStrawPlot.pdf"),width = 10,height = 6)
# pjs
# dev.off()


# pelb<- ElbowPlot(st,ndims = 50) #碎石图
# pdf(paste0(output_dir,project.name,"_ElbowPlot.pdf"),width = 10,height = 6)
# pelb
# dev.off()


#方法二：SCTransform方法
#SCTranform的归一化方式无法执行Jackstraw置换检验
st <- SCTransform(st,assay = "Spatial")  #12min   #加载glmGamPoi后2min
st <- RunPCA(st)
pelb <- ElbowPlot(st,ndims = 50) #碎石图
pdf(paste0(output_dir,project.name,"_01_ElbowPlot.pdf"),width = 10,height = 6)
pelb
dev.off()

# 查看当前seurat对象中包含的assay
st@assays
# 查看当前激活的assay
DefaultAssay(st)
# 修改默认assay
DefaultAssay(st) <- "SCT"
# DefaultAssay(st) <- "Spatial"

#### 3. spot聚类和分群
#基于PCA对Spot进行聚类
st <- FindNeighbors(st,dims = 1:20)
# plan("multicore")有bug，改为plan("multisession")，但运行时间较长
st <- FindClusters(st,resolution = 0.4)


#### 4. 分群结果及基因表达量可视化
#使用组织切片图片来呈现分群结果
#col <- c("#00CED1","#32CD32","#FFA500","#FF8C00","#FF4500","#FA8072","#FF0000","#8B0000","#FF1493","#FF00FF","#FF00FF","#800080","#9400D3","#8A2BE2","#0000FF")
#col1 <- c("#005EEC","#009BB1","#50D24A","#A0CD1E","#D0CA0F","#FFC700","#FFA101","#FF7B01","#FF4601","#FE1100","#FF4867","#FF7497","#EF3AA2","#DE00AD","#DE005D","#8A2BE2","#0000FF")
p1 <- SpatialDimPlot(st,label = T,label.size = 3)
p1
pdf(paste0(output_dir,project.name,"_02_SpatialDimPlot.pdf"),width = 10,height = 6)
p1
dev.off()
#除了组织切片，数据也可以使用tSNE非线性降维进行可视化
st <- RunTSNE(st,dims = 1:30,label=T)
st <- RunUMAP(st, dims = 1:30, label = T)
#可以保存spot的坐标信息，用于导入loupe中
# data.tsne <- st@reductions[["tsne"]]@cell.embeddings
# data.umap <- st@reductions[["umap"]]@cell.embeddings
# data.cluster <- st[["seurat_clusters"]]
# write.csv(data.umap,paste0(output_dir,project.name,"_data.umap.csv"),row.names = T,quote = T)
# write.csv(data.tsne,paste0(output_dir,project.name,"_data.tsne.csv"),row.names = T,quote = T)
# write.csv(data.cluster,paste0(output_dir,project.name,"_data.cluster.csv"),row.names = T,quote = T)
p2 <- DimPlot(st,reduction = "tsne",pt.size = 1.5,label = T)
p2
p3 <- DimPlot(st, reduction = "umap",pt.size = 1.5,label = T)
p3
p2+p3
pdf(paste0(output_dir,project.name,"_03_tsne_umap.pdf"),width = 10,height = 5)
p2+p3
dev.off()

#亚群和样本的细胞分布
# 更改样本名称
# a <- c(rep(0,ncol(st)))
# for (i in 1:ncol(st)) {
#  a[i] <- strsplit(colnames(st),"_")[[i]][1]
# }
# st@meta.data[["orig.ident"]] <- a
pdf(paste0(output_dir,project.name,"_04_umap_cluster_sample.pdf"),width = 10,height = 5)
DimPlot(st,reduction = "umap",group.by = c("seurat_clusters","orig.ident"),pt.size = 1.2)
dev.off()

pdf(paste0(output_dir,project.name,"_05_tsne_cluster_sample.pdf"),width = 10,height = 5)
DimPlot(st,reduction = "tsne",group.by = c("seurat_clusters","orig.ident"),pt.size = 1.2)
dev.off()

######堆叠图绘制
#提取细胞在细胞亚群和样本的分布数据
stat <- table(st$orig.ident,st$seurat_clusters)
stat <- as.data.frame(stat)
colnames(stat) <- c('sample','seurat_clusters','Freq')

# 也可以直接tidyr包中的count函数来统计组合频数，代码更简洁
# stat <- st[[]] |>
#   count(orig.ident, seurat_clusters) |>      # 直接统计组合频数
#   rename(sample = orig.ident, Freq = n)
# 绘制基础的数量堆叠图形
p1<-ggplot(data = stat,aes(x = sample,y = Freq, fill = seurat_clusters))+
  geom_bar(stat = 'identity',position = 'fill')
# position='stack'是数量堆叠图,positon='fill'是比例堆叠图
p1

#更改配色、坐标轴名称、柱子的大小
color<-c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
         "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
         "#756455","#005EEC","#009BB1")
p2<-ggplot(data = stat,aes(x = sample,y = Freq, fill = seurat_clusters))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.5)+
  labs(x = "Sample",y = "Cell Number")+
  scale_fill_manual(values = color)+
  guides(fill = guide_legend(ncol = 1,bycol = T,override.aes = list(size = 5)))
p2

#自定义主题
mytheme <- theme(panel.background = element_rect(fill = 'white',color = 'black'),
                  panel.grid = element_line(color = 'lightgrey'),
                  axis.title.y = element_text(face = 'bold',color = 'black',size = 14),
                  axis.title.x = element_text(face = 'bold',color = 'black',size = 14),
                  axis.text.y = element_text(colour = 'black',size = 12),
                  axis.text.x = element_text(color = 'black',size = 12,angle = 45,vjust = 0.6),
                  legend.title = element_blank(),
                  legend.text = element_text(color = 'black',size = 14))
p3<-p2+mytheme
p3

#保存图片
pdf(paste0(output_dir,project.name,"_06_cluster_sample_barplot.pdf"),width = 10,height = 6)
p3
dev.off()

###############################################################################

#可以使用交互界面同时呈现两种降维方式，以完成对应性查看
LinkedDimPlot(st,reduction = "umap",image = "anterior")
LinkedDimPlot(st,reduction = "tsne",image = "posterior")
LinkedFeaturePlot(st, feature = "Stx1a",image = "anterior")
LinkedFeaturePlot(st, feature = "Car8",image = "posterior",reduction = "tsne")
#高亮cluster 1的spot
SpatialDimPlot(st,cells.highlight = CellsByIdentities(st,idents = "1"))
#如果要高亮多个cluster，需要将facet.highlight参数改为TRUE，否则会报错
SpatialDimPlot(st,cells.highlight = CellsByIdentities(st,idents = c("0","4","16")),
               images = "anterior",facet.highlight = T)
SpatialDimPlot(st,cells.highlight = CellsByIdentities(st,idents = c("0","4","16")),
               images = "posterior",facet.highlight = T)

#### 5. 差异基因分析和图形可视化
##（1）亚群间差异基因分析
# st <- JoinLayers(st,assay = "Spatial")  ##合并矩阵
# 准备 SCT 差异分析（统一多个模型）
st <- PrepSCTFindMarkers(st)
save(st,file = "ST.Rda")
load("ST.Rda")
dif <- FindAllMarkers(st,assay = "SCT",logfc.threshold = 0.6,min.pct = 0.4,only.pos = T)

# logfc.threshold定义上调倍数阈值，min.pct定义基因至少在细胞亚群中多少细胞中表达，only.pos确定只筛选上调基因
write.table(dif,paste0(output_dir,"diff.txt"),sep = "\t",quote = F,row.names = T,col.names = T)
dif <- read.table(paste0(output_dir,"diff.txt"),header = T,row.names = 1,sep = "\t")
# 通过dplyr包来完成对top基因的筛选，其中|>与linux中的管道符“|”的功能相同
sig.dif <- dif |> 
  group_by(cluster) |> 
  top_n(n = 5,wt = avg_log2FC)
write.table(sig.dif,paste0(output_dir,"diff.sig.txt"),sep = "\t",quote = F,row.names = T,col.names = T)

# 差异基因
genes <- unique(sig.dif$gene)
length(genes);genes

# 差异基因的图形可视化
# 气泡热图（一张大图）
DotPlot(st,features =genes[26:45]) + theme(axis.text.x = element_text(angle = 45,hjust = 1))
# 黄紫热图（一张大图）
DoHeatmap(st,features = genes,angle=45,assay = "SCT")
pdf(paste0(output_dir,project.name,"_06_doheatmap.pdf"),width = 15,height = 9)
DoHeatmap(st,features = genes,angle=45,assay = "SCT")
dev.off()
# umap热图和tSNE热图（每个基因一个小图）
FeaturePlot(st, features = genes[16:27],ncol = 3) #基因表达umap映射图
# 小提琴表达量图（每个基因一个小图）
VlnPlot(st, features = genes[75:79],pt.size=0) #小提琴图
# 组织映射图（每个基因一个小图）
SpatialFeaturePlot(st,features = genes[34],pt.size=1.5,alpha=c(0.1,1),interactive = F) #组织映射图
# 组织映射图（每个基因一个小图，交互式）
SpatialPlot(st,features = genes[34],pt.size=1.5,alpha=c(0.1,1),images = "posterior",interactive = T) #组织映射图

##（2）两个亚群之间的差异基因分析
# 第三和第十五亚群之间的差异基因分析
group.cluster.dif <- FindMarkers(st,assay = "SCT",ident.1 = "3",ident.2 = "15",
                       logfc.threshold = 0.6,min.pct = 0.4)
# 取绝对值最大的前10个差异基因
group.cluster.sig.dif <- group.cluster.dif |> 
  top_n(n = 10,wt = abs(avg_log2FC)) |> 
  mutate(diff.pct = pct.1 - pct.2)
# 取正调控top10和负调控top10
group.cluster.sig.dif1 <- rbind(group.cluster.dif |> top_n(n = 10,wt = avg_log2FC),
                              group.cluster.dif |> top_n(n = -10,wt = avg_log2FC))
write.table(group.cluster.sig.dif,paste0(output_dir,"group.cluster.sig.dif.txt"),sep = "\t",quote = F,row.names = T,col.names = T)
write.table(group.cluster.sig.dif1,paste0(output_dir,"group.cluster.sig.diftop10.txt"),sep = "\t",quote = F,row.names = T,col.names = T)

#（3）不同样本之间的亚群差异基因分析
# a <- c(rep(0,ncol(st)))
# for (i in 1:ncol(st)){
#  a[i] <- paste(strsplit(colnames(st),"_")[[i]][1],as.character(st@meta.data$seurat_clusters[i]),sep = "_")}
# st[["group.cluster"]] <- a

# 也可以直接在meta.data中创建一个新的列来存储样本和亚群的组合信息，代码更简洁
st[["group.cluster"]] <- paste(st$orig.ident,Idents(st),sep = '_')
# 两个样本之间组织比较
group.sample.dif <- FindMarkers(
    st,
    assay = "SCT",
    group.by = "group.cluster",
    ident.1 = "anterior_1",
    ident.2 = "posterior_1",
    logfc.threshold = 0.6,
    min.pct = 0.4
)
group.sample.sig.dif <- group.sample.dif |> top_n(n = 10, wt = abs(avg_log2FC)) |> mutate(diff.pct = pct.1 - pct.2)
group.sample.sig.dif1 <- rbind(group.sample.dif |> top_n(n = 10, wt = avg_log2FC), group.sample.dif |> top_n(n = -10, wt = avg_log2FC))
write.table(group.sample.sig.dif, paste0(output_dir,"group.sample.sig.dif.txt"), sep = "\t", quote = F, row.names = T, col.names = T)
write.table(group.sample.sig.dif1, paste0(output_dir,"group.sample.sig.dif1.txt"), sep = "\t", quote = F, row.names = T, col.names = T)

save(st,file = paste0(output_dir,"ST.Rda"))

#################################################################################
#局部区域筛选及分析
rm(list = ls())
gc()
load(paste0(output_dir,"ST.Rda"))

##6. 局部区域筛选步骤
#（1）情况一：提取亚群区域
C1_C3 <- subset(st,idents = c("1","3"))
SpatialDimPlot(C1_C3, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)

#（2）情况二：提取矩形区域
#根据样本拆分seurat对象
anterior <- SplitObject(st,split.by = "orig.ident")[[1]]

#将图形坐标提取到meat.data中，方便局部区域筛选
# VisiumV2版本的坐标提取方法
# 方法一：使用 GetTissueCoordinates()
coords <- GetTissueCoordinates(anterior)
anterior[["x"]] <- coords[["x"]]
anterior[["y"]] <- coords[["y"]]
# 方法二：直接从图像对象中提取坐标
anterior[["x"]] <- anterior@images$anterior@boundaries$centroids@coords[,1]
anterior[["y"]] <- anterior@images$anterior@boundaries$centroids@coords[,2]
# 方法三：Extract and attach spatial coordinates from anterior1 image
centroids <- anterior@images$anterior@boundaries$centroids
coords <- setNames(as.data.frame(centroids@coords), c("x", "y"))
rownames(coords) <- centroids@cells
anterior$x <- coords[colnames(anterior), "x"]
anterior$y <- coords[colnames(anterior), "y"]

# VisiumV1版本的坐标提取方法
anterior[["x"]] <- anterior@images[["anterior"]]@coordinates[["imagerow"]]
anterior[["y"]] <- anterior@images[["anterior"]]@coordinates[["imagecol"]]

#提取方法
cortex <- subset(anterior,x > 5000 & y < 9000) # 根据坐标提取行大于5000且列小于9000的区域
cortex <- subset(cortex,x < 9000 & y > 5000) # 根据坐标提取上一步提取的行小于9000且列大于5000的局部区域
#从全局查看提取结果
cortex@images[["posterior"]] <- NULL # 删除posterior图像
SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 2, label.size = 3)# crop参数设置为TRUE可以只显示提取的局部区域，设置为FALSE可以同时显示全局图像和提取的局部区域

#数据归一化
cortex <- NormalizeData(cortex,normalization = "LogNormalize")
cortex <- FindVariableFeatures(cortex,selection.method = "vst",verbose = F)
cortex <- ScaleData(cortex,features = rownames(cortex))
#线性降维（PCA）,默认用高变基因集，但也可通过 features 参数自己指定
#同时也可对主成分数量进行定义 
cortex <- RunPCA(cortex,npcs = 50)
ElbowPlot(cortex)

# 第二种数据归一化的方法用 SCTransform 
cortex <- SCTransform(cortex,assay = "Spatial")  #1min
cortex <- RunPCA(cortex)
#基于PCA对Spot进行聚类
cortex <- FindNeighbors(cortex,dims = 1:10)
cortex <- FindClusters(cortex,resolution = 0.4)

ElbowPlot(cortex)

#查看局部区域的分群结果
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)

pdf(paste0(output_dir,project.name,"_07_cortex_spatial_plots.pdf"),width = 10,height = 6)
p1+p2
dev.off()

#使用Tsne和umap图展示局部区域的细胞分群结果
cortex <- RunTSNE(cortex,dims = 1:30,label=T)
p3 <- DimPlot(cortex,reduction = "tsne",pt.size = 1.5,label = T)

cortex <- RunUMAP(cortex, dims = 1:30, label = T)
p4 <- DimPlot(cortex, reduction = "umap",pt.size = 1.5,label = T)


pdf(paste0(output_dir,project.name,"_08_cortex_tsne_umap.pdf"),width = 10,height = 5)
p3+p4
dev.off()

# save(cortex,file = paste0(output_dir,"cortex.Rda"))

# 提取0 1亚群的spot进行可视化
cortex_C1_C2 <- subset(cortex,idents = c("0","1"))
SpatialDimPlot(cortex_C1_C2, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
