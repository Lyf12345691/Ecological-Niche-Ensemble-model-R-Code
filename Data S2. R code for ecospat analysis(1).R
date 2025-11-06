setwd("D:/eco")
#加载包
library(ENMeval)
library(raster)
library(dplyr)
library(tidyverse)
library(MASS)
library(rrcov)
library(ggplot2)
library(cowplot)
library(plyr)
library(Cairo)
library(showtext)
library(spocc)
library(raster)
library(maptools)
library(rgeos)
library(dismo)
library(ecospat)
library(xlsx)
library(raster)
library(ade4)
library(ape)
library(biomod2)
#1.提取背景点

envs <- raster::stack("D:/bg.asc")
bg <- randomPoints(envs, n=10000) 
bg <- as.data.frame(bg)
write.csv( bg,"D:/bg.csv")

#2.提取环境数据点插值

#加载分布点数据

occs <- read.csv("bg.csv")

occs<-occs[,c(2,3)]  %>%  as.matrix()
head(occs)
#加载current环境数据
envs.files <- list.files(path=paste("D:/asc/current", sep=''),  
                         pattern='asc', full.names=TRUE)
envs <- raster::stack(envs.files)
#提取发生点的变量值
icurrent <- raster::extract(envs, occs)
write.csv( icurrent,"D:/eco/icurrent.csv")
#加载2050s环境数据
envs.files <- list.files(path=paste("D:/asc/50126", sep=''),  
                         pattern='asc', full.names=TRUE)
envs <- raster::stack(envs.files)
#提取发生点的变量值
i50126 <- raster::extract(envs, occs)
write.csv(i50126,"D:/eco/i50126.csv")


envs.files <- list.files(path=paste("D:/asc/50585", sep=''),  
                         pattern='asc', full.names=TRUE)
envs <- raster::stack(envs.files)
#提取发生点的变量值
i50585 <- raster::extract(envs, occs)
write.csv(i50585,"D:/eco/i50585.csv")

#加载2090s环境数据
envs.files <- list.files(path=paste("D:/asc/90126", sep=''),  
                         pattern='asc', full.names=TRUE)
envs <- raster::stack(envs.files)
#提取发生点的变量值
i90126 <- raster::extract(envs, occs)
write.csv( i90126,"D:/eco/i90126.csv")


envs.files <- list.files(path=paste("D:/asc/90585", sep=''),  
                         pattern='asc', full.names=TRUE)
envs <- raster::stack(envs.files)
#提取发生点的变量值
i90585 <- raster::extract(envs, occs)
write.csv( i90585,"D:/eco/i90585.csv")
#联合分布点及点插值数据
head(icurrent)
#调整列名
colnames(icurrent) <- c("bio12", "bio14", "bio3", "bio7","projcurrent")
colnames(i50126) <- c("bio12", "bio14", "bio3",  "proj50126")
colnames(i50585) <- c("bio12", "bio14", "bio3", "bio7", "proj50585")
colnames(i90126) <- c("bio12", "bio14", "bio3", "bio7", "proj90126")
colnames(i90585) <- c("bio12", "bio14", "bio3", "bio7", "proj90585")

#联合分布点及点插值数据
ecospatcurrent <- cbind(occs,icurrent)
ecospat50126 <- cbind(occs,i50126)
ecospat50585 <- cbind(occs,i50585)
ecospat90126 <- cbind(occs,i90126)
ecospat90585 <- cbind(occs,i90585)

#提取二值化数据点插值

envs.files <- list.files(path=paste("D:/binoutasc", sep=''),  
                         pattern='asc', full.names=TRUE)
envs <- raster::stack(envs.files)

binenvs <- raster::extract(envs, occs)

head(binenvs)
#联合分布点及点插值数据
ecospatcurrent <- cbind(occs,icurrent,binenvs[,7])
ecospat50126 <- cbind(occs,i50126,binenvs[,1])
ecospat50585 <- cbind(occs,i50585,binenvs[,2])
ecospat90126 <- cbind(occs,i90126,binenvs[,3])
ecospat90585 <- cbind(occs,i90585,binenvs[,4])
head(ecospatcurrent)
#调整列名
colnames(ecospatcurrent) <- c("x", "y", "bio12", "bio14", "bio3", "bio7","projcurrent","species_occ")
colnames(ecospat50126) <- c("x", "y", "bio12", "bio14", "bio3", "bio7","proj50126","species_occ")
colnames(ecospat50585) <- c("x", "y", "bio12", "bio14", "bio3", "bio7", "proj50585","species_occ")
colnames(ecospat90126) <- c("x", "y", "bio12", "bio14", "bio3", "bio7", "proj90126","species_occ")
colnames(ecospat90585) <- c("x", "y", "bio12", "bio14", "bio3", "bio7", "proj90585","species_occ")


#删除有空值的行
#删除有空值的行
ecospatcurrent <- ecospatcurrent[complete.cases(ecospatcurrent), ]
ecospat50126 <- ecospat50126[complete.cases(ecospat50126), ]

ecospat50585 <- ecospat50585[complete.cases(ecospat50585), ]
ecospat90126 <- ecospat90126[complete.cases(ecospat90126), ]

ecospat90585 <- ecospat90585[complete.cases(ecospat90585), ]


#写出数据
write.csv(ecospatcurrent,"D:/eco/ecospatcurrent.csv")
write.csv(ecospat50126,"D:/eco/ecospat50126.csv")
write.csv(ecospat50585,"D:/eco/ecospat50585.csv")
write.csv(ecospat90126,"D:/eco/ecospat90126.csv")
write.csv(ecospat90585,"D:/eco/ecospat90585.csv")


head(ecospatcurrent)
head(ecospat50126)
 #50126
pca.env <- dudi.pca(rbind(ecospatcurrent,ecospat50126)[,3:6],scannf=F,nf=2)

#Plot Variables Contribution

ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)      #相关圆表示原始预测值对PCA轴的贡献。

pdf("PCA50126.pdf", width = 6, height = 6)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) 
dev.off()

#[,3:6]为环境因子所在例

#预测轴上的分数

# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for 本地分布
scores.sp.ecospatcurrent <- suprow(pca.env,ecospatcurrent[which(ecospatcurrent[,8]==1),3:6])$li

#3:13  为选择的环境因子所在列
#(nat[,10]==1) 为nat第11列，存在/不存在1/0，==1表示选择存在

# PCA scores for 物种入侵分布
scores.sp.ecospat50126 <- suprow(pca.env,ecospat50126[which(ecospat50126[,8]==1),3:6])$li

# PCA scores for 整个本地研究区域
scores.clim.ecospatcurrent <- suprow(pca.env,ecospatcurrent[,3:6])$li

# PCA scores for 整个物种入侵研究区域
scores.clim.ecospat50126 <- suprow(pca.env,ecospat50126[,3:6])$li

#计算出现密度 Grid

# gridding 本地生态位
grid.clim.ecospatcurrent <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                                  glob1=scores.clim.ecospatcurrent,
                                                  sp=scores.sp.ecospatcurrent, R=100,
                                                  th.sp=0)

#R为 网格的分辨率。

# gridding 入侵生态位
grid.clim.ecospat50126 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                                glob1=scores.clim.ecospat50126,
                                                sp=scores.sp.ecospat50126, R=100,
                                                th.sp=0)

#计算生态位重叠
D.overlap <- ecospat.niche.overlap (grid.clim.ecospatcurrent, grid.clim.ecospat50126, cor = TRUE)$D
D.overlap
write.csv(D.overlap,"D:/eco/D.overlap50126.csv")

# D    Schoener's D   可改为I

#建议使用至少1000个重复进行等效性检验。作为一个例子，我们使用了rep=10，以减少计算时间。

eq.test <- ecospat.niche.equivalency.test(grid.clim.ecospatcurrent, grid.clim.ecospat50126,
                                          rep=1000)
#Plot等效性检验
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
pdf("equivalency50126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
dev.off()

#生态位相似性检验
#建议使用至少1000个重复进行相似性检验

sim.test <- ecospat.niche.similarity.test(grid.clim.ecospatcurrent, grid.clim.ecospat50126,
                                          rep=1000)
#Plot相似性检验
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
pdf("Similarity50126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
dev.off()

#画出生态位扩展、稳定性和未填充的地块相似性检验
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
pdf("expansion_Similarity50126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
dev.off()

ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
pdf("stability_Similarity50126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
dev.off()

ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
pdf("unfilling_Similarity50126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
dev.off()


#####在模拟气候中划定生态位类别和量化生态位动态，最主要

niche.dyn <- ecospat.niche.dyn.index (grid.clim.ecospatcurrent, grid.clim.ecospat50126, intersection = 0.1)

#intersection  用于消除边缘气候的环境密度分位数。如果intersection=NA，对整个环境范围进行分析（本地和入侵）。如果intersection=0，则在本地范围和入侵范围之间的中间部分执行分析。如果intersection=0.05，则分析在本地和入侵的第5个分位数的交点处执行环境密度。

###可视化生态位类别、生态位动态和范围之间的气候类比

ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat50126, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
ecospat.shift.centroids(scores.sp.ecospatcurrent, scores.sp.ecospat50126, scores.clim.ecospatcurrent, scores.clim.ecospat50126)

pdf("Niche Overlap50126.pdf", width = 6, height = 6)
ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat50126, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
ecospat.shift.centroids(scores.sp.ecospatcurrent, scores.sp.ecospat50126, scores.clim.ecospatcurrent, scores.clim.ecospat50126)
dev.off()

# 计算生态位动态，量化增加、减少和维持不变的部分
niche.dyn <- ecospat.niche.dyn.index(grid.clim.ecospatcurrent, grid.clim.ecospat50126, intersection = 0.1)

# 提取动态信息：expansion (增加), stability (维持不变), unfilling (减少)
expansion <- niche.dyn$dynamic.index.w['expansion']  # 生态位扩展
stability <- niche.dyn$dynamic.index.w['stability']  # 生态位维持不变
unfilling <- niche.dyn$dynamic.index.w['unfilling']  # 生态位未填充

# 打印结果
print(paste("生态位增加 (Expansion):", expansion))
print(paste("生态位维持不变 (Stability):", stability))
print(paste("生态位减少 (Unfilling):", unfilling))

# 保存结果到CSV文件
results <- data.frame(
  Metric = c("Expansion", "Stability", "Unfilling"),
  Value = c(expansion, stability, unfilling)
)

write.csv(results, "D:/eco/Niche_Dynamics_50126.csv", row.names = FALSE)

# 可视化生态位扩展、维持不变和未填充的部分
pdf("Niche_Dynamics_50126.pdf", width = 6, height = 6)
ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat50126, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
dev.off()


#quant=0.25 用于划分边缘气候的环境密度分位数
#选择要绘制的密度：如果interest=1，绘制本机密度，如果interest=2，绘制侵入密度图。

#3.生态位动态可视化


#50585

pca.env <- dudi.pca(rbind(ecospatcurrent,ecospat50585)[,3:6],scannf=F,nf=2)

#Plot Variables Contribution

ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)      #相关圆表示原始预测值对PCA轴的贡献。
pdf("PCA50585.pdf", width = 6, height = 6)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) 
dev.off()

#[,3:6]为环境因子所在例

#预测轴上的分数

# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for 本地分布
scores.sp.ecospatcurrent <- suprow(pca.env,ecospatcurrent[which(ecospatcurrent[,8]==1),3:6])$li

#3:6  为选择的环境因子所在列
#(nat[,10]==1) 为nat第11列，存在/不存在1/0，==1表示选择存在

# PCA scores for 物种入侵分布
scores.sp.ecospat50585 <- suprow(pca.env,ecospat50585[which(ecospat50585[,8]==1),3:6])$li

# PCA scores for 整个本地研究区域
scores.clim.ecospatcurrent <- suprow(pca.env,ecospatcurrent[,3:6])$li

# PCA scores for 整个物种入侵研究区域
scores.clim.ecospat50585 <- suprow(pca.env,ecospat50585[,3:6])$li

#计算出现密度 Grid

# gridding 本地生态位
grid.clim.ecospatcurrent <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                                  glob1=scores.clim.ecospatcurrent,
                                                  sp=scores.sp.ecospatcurrent, R=100,
                                                  th.sp=0)

#R为 网格的分辨率。

# gridding 入侵生态位
grid.clim.ecospat50585 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                                glob1=scores.clim.ecospat50585,
                                                sp=scores.sp.ecospat50585, R=100,
                                                th.sp=0)

#计算生态位重叠
D.overlap <- ecospat.niche.overlap (grid.clim.ecospatcurrent, grid.clim.ecospat50585, cor = TRUE)$D
D.overlap
write.csv(D.overlap,"D:/JEM/baomao/eco/D.overlap50585.csv")

# D    Schoener's D   可改为I

#建议使用至少1000个重复进行等效性检验。作为一个例子，我们使用了rep=10，以减少计算时间。

eq.test <- ecospat.niche.equivalency.test(grid.clim.ecospatcurrent, grid.clim.ecospat50585,
                                          rep=1000)
#Plot等效性检验
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
pdf("equivalency50585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
dev.off()

#生态位相似性检验
#建议使用至少1000个重复进行相似性检验

sim.test <- ecospat.niche.similarity.test(grid.clim.ecospatcurrent, grid.clim.ecospat50585,
                                          rep=1000)
#Plot相似性检验
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
pdf("Similarity50585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
dev.off()

#画出生态位扩展、稳定性和未填充的地块相似性检验
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
pdf("expansion_Similarity50585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
dev.off()

ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
pdf("stability_Similarity50585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
dev.off()

ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
pdf("unfilling_Similarity50585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
dev.off()



#####在模拟气候中划定生态位类别和量化生态位动态，最主要

niche.dyn <- ecospat.niche.dyn.index (grid.clim.ecospatcurrent, grid.clim.ecospat50585, intersection = 0.1)

#intersection  用于消除边缘气候的环境密度分位数。如果intersection=NA，对整个环境范围进行分析（本地和入侵）。如果intersection=0，则在本地范围和入侵范围之间的中间部分执行分析。如果intersection=0.05，则分析在本地和入侵的第5个分位数的交点处执行环境密度。

###可视化生态位类别、生态位动态和范围之间的气候类比

ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat50585, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
ecospat.shift.centroids(scores.sp.ecospatcurrent, scores.sp.ecospat50585, scores.clim.ecospatcurrent, scores.clim.ecospat50585)
pdf("Niche Overlap50585.pdf", width = 6, height = 6)
ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat50585, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
ecospat.shift.centroids(scores.sp.ecospatcurrent, scores.sp.ecospat50585, scores.clim.ecospatcurrent, scores.clim.ecospat50585)
dev.off()

# 计算生态位动态，量化增加、减少和维持不变的部分
niche.dyn <- ecospat.niche.dyn.index(grid.clim.ecospatcurrent, grid.clim.ecospat50585, intersection = 0.1)

# 提取动态信息：expansion (增加), stability (维持不变), unfilling (减少)
expansion <- niche.dyn$dynamic.index.w['expansion']  # 生态位扩展
stability <- niche.dyn$dynamic.index.w['stability']  # 生态位维持不变
unfilling <- niche.dyn$dynamic.index.w['unfilling']  # 生态位未填充

# 打印结果
print(paste("生态位增加 (Expansion):", expansion))
print(paste("生态位维持不变 (Stability):", stability))
print(paste("生态位减少 (Unfilling):", unfilling))

# 保存结果到CSV文件
results <- data.frame(
  Metric = c("Expansion", "Stability", "Unfilling"),
  Value = c(expansion, stability, unfilling)
)

write.csv(results, "D:/eco/Niche_Dynamics_50585.csv", row.names = FALSE)

# 可视化生态位扩展、维持不变和未填充的部分
pdf("Niche_Dynamics_50585.pdf", width = 6, height = 6)
ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat50585, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
dev.off()

#90126
pca.env <- dudi.pca(rbind(ecospatcurrent,ecospat90126)[,3:6],scannf=F,nf=2)

#Plot Variables Contribution

ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)      #相关圆表示原始预测值对PCA轴的贡献。
pdf("PCA90126.pdf", width = 6, height = 6)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) 
dev.off()

#[,3:8]为环境因子所在例

#预测轴上的分数

# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for 本地分布
scores.sp.ecospatcurrent <- suprow(pca.env,ecospatcurrent[which(ecospatcurrent[,8]==1),3:6])$li

#3:8  为选择的环境因子所在列
#(nat[,10]==1) 为nat第11列，存在/不存在1/0，==1表示选择存在

# PCA scores for 物种入侵分布
scores.sp.ecospat90126 <- suprow(pca.env,ecospat90126[which(ecospat90126[,8]==1),3:6])$li

# PCA scores for 整个本地研究区域
scores.clim.ecospatcurrent <- suprow(pca.env,ecospatcurrent[,3:6])$li

# PCA scores for 整个物种入侵研究区域
scores.clim.ecospat90126 <- suprow(pca.env,ecospat90126[,3:6])$li

#计算出现密度 Grid

# gridding 本地生态位
grid.clim.ecospatcurrent <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                                  glob1=scores.clim.ecospatcurrent,
                                                  sp=scores.sp.ecospatcurrent, R=100,
                                                  th.sp=0)

#R为 网格的分辨率。

# gridding 入侵生态位
grid.clim.ecospat90126 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                                glob1=scores.clim.ecospat90126,
                                                sp=scores.sp.ecospat90126, R=100,
                                                th.sp=0)

#计算生态位重叠
D.overlap <- ecospat.niche.overlap (grid.clim.ecospatcurrent, grid.clim.ecospat90126, cor = TRUE)$D
D.overlap
write.csv(D.overlap,"D:/eco/D.overlap90126.csv")
# D    Schoener's D   可改为I

#建议使用至少1000个重复进行等效性检验。作为一个例子，我们使用了rep=10，以减少计算时间。

eq.test <- ecospat.niche.equivalency.test(grid.clim.ecospatcurrent, grid.clim.ecospat90126,
                                          rep=1000)
#Plot等效性检验
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
pdf("equivalency90126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
dev.off()

#生态位相似性检验
#建议使用至少1000个重复进行相似性检验

sim.test <- ecospat.niche.similarity.test(grid.clim.ecospatcurrent, grid.clim.ecospat90126,
                                          rep=1000)
#Plot相似性检验
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
pdf("Similarity90126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
dev.off()
#画出生态位扩展、稳定性和未填充的地块相似性检验
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
pdf("expansion_Similarity90126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
dev.off()

ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
pdf("stability_Similarity90126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
dev.off()

ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
pdf("unfilling_Similarity90126.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
dev.off()



#####在模拟气候中划定生态位类别和量化生态位动态，最主要

niche.dyn <- ecospat.niche.dyn.index (grid.clim.ecospatcurrent, grid.clim.ecospat90126, intersection = 0.1)

#intersection  用于消除边缘气候的环境密度分位数。如果intersection=NA，对整个环境范围进行分析（本地和入侵）。如果intersection=0，则在本地范围和入侵范围之间的中间部分执行分析。如果intersection=0.05，则分析在本地和入侵的第5个分位数的交点处执行环境密度。

###可视化生态位类别、生态位动态和范围之间的气候类比

ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat90126, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
ecospat.shift.centroids(scores.sp.ecospatcurrent, scores.sp.ecospat90126, scores.clim.ecospatcurrent, scores.clim.ecospat90126)
pdf("Niche Overlap90126.pdf", width = 6, height = 6)
ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat90126, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
ecospat.shift.centroids(scores.sp.ecospatcurrent, scores.sp.ecospat90126, scores.clim.ecospatcurrent, scores.clim.ecospat90126)
dev.off()
#quant=0.25 用于划分边缘气候的环境密度分位数
#选择要绘制的密度：如果interest=1，绘制本机密度，如果interest=2，绘制侵入密度图。

# 计算生态位动态，量化增加、减少和维持不变的部分
niche.dyn <- ecospat.niche.dyn.index(grid.clim.ecospatcurrent, grid.clim.ecospat90126, intersection = 0.1)

# 提取动态信息：expansion (增加), stability (维持不变), unfilling (减少)
expansion <- niche.dyn$dynamic.index.w['expansion']  # 生态位扩展
stability <- niche.dyn$dynamic.index.w['stability']  # 生态位维持不变
unfilling <- niche.dyn$dynamic.index.w['unfilling']  # 生态位未填充

# 打印结果
print(paste("生态位增加 (Expansion):", expansion))
print(paste("生态位维持不变 (Stability):", stability))
print(paste("生态位减少 (Unfilling):", unfilling))

# 保存结果到CSV文件
results <- data.frame(
  Metric = c("Expansion", "Stability", "Unfilling"),
  Value = c(expansion, stability, unfilling)
)

write.csv(results, "D:/eco/Niche_Dynamics_90126.csv", row.names = FALSE)


# 可视化生态位扩展、维持不变和未填充的部分
pdf("Niche_Dynamics_90126.pdf", width = 6, height = 6)
ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat90126, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
dev.off()


#90585

pca.env <- dudi.pca(rbind(ecospatcurrent,ecospat90585)[,3:6],scannf=F,nf=2)

#Plot Variables Contribution

ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)      #相关圆表示原始预测值对PCA轴的贡献。
pdf("PCA90585.pdf", width = 6, height = 6)
ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) 
dev.off()
#[,3:8]为环境因子所在例

#预测轴上的分数

# PCA scores for the whole study area
scores.globclim <- pca.env$li

# PCA scores for 本地分布
scores.sp.ecospatcurrent <- suprow(pca.env,ecospatcurrent[which(ecospatcurrent[,8]==1),3:6])$li

#3:8  为选择的环境因子所在列
#(nat[,10]==1) 为nat第11列，存在/不存在1/0，==1表示选择存在

# PCA scores for 物种入侵分布
scores.sp.ecospat90585 <- suprow(pca.env,ecospat90585[which(ecospat90585[,8]==1),3:6])$li

# PCA scores for 整个本地研究区域
scores.clim.ecospatcurrent <- suprow(pca.env,ecospatcurrent[,3:6])$li

# PCA scores for 整个物种入侵研究区域
scores.clim.ecospat90585 <- suprow(pca.env,ecospat90585[,3:6])$li

#计算出现密度 Grid

# gridding 本地生态位
grid.clim.ecospatcurrent <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                                  glob1=scores.clim.ecospatcurrent,
                                                  sp=scores.sp.ecospatcurrent, R=100,
                                                  th.sp=0)

#R为 网格的分辨率。

# gridding 入侵生态位
grid.clim.ecospat90585 <- ecospat.grid.clim.dyn(glob=scores.globclim,
                                                glob1=scores.clim.ecospat90585,
                                                sp=scores.sp.ecospat90585, R=100,
                                                th.sp=0)

#计算生态位重叠
D.overlap <- ecospat.niche.overlap (grid.clim.ecospatcurrent, grid.clim.ecospat90585, cor = TRUE)$D
D.overlap
write.csv(D.overlap,"D:/JEM/baomao/eco/D.overlap90585.csv")

# D    Schoener's D   可改为I

#建议使用至少1000个重复进行等效性检验。作为一个例子，我们使用了rep=10，以减少计算时间。

eq.test <- ecospat.niche.equivalency.test(grid.clim.ecospatcurrent, grid.clim.ecospat90585,
                                          rep=1000)
#Plot等效性检验
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
pdf("equivalency90585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
dev.off()

#生态位相似性检验
sim.test <- ecospat.niche.similarity.test(grid.clim.ecospatcurrent, grid.clim.ecospat90585,
                                          rep=1000)
#Plot相似性检验
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
pdf("Similarity90585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
dev.off()

#画出生态位扩展、稳定性和未填充的地块相似性检验
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
pdf("expansion_Similarity90585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "expansion", "Similarity")
dev.off()

ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
pdf("stability_Similarity90585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "stability", "Similarity")
dev.off()

ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
pdf("unfilling_Similarity90585.pdf", width = 6, height = 6)
ecospat.plot.overlap.test(sim.test, "unfilling", "Similarity")
dev.off()


#####在模拟气候中划定生态位类别和量化生态位动态，最主要

niche.dyn <- ecospat.niche.dyn.index (grid.clim.ecospatcurrent, grid.clim.ecospat90585, intersection = 0.1)

#intersection  用于消除边缘气候的环境密度分位数。如果intersection=NA，对整个环境范围进行分析（本地和入侵）。如果intersection=0，则在本地范围和入侵范围之间的中间部分执行分析。如果intersection=0.05，则分析在本地和入侵的第5个分位数的交点处执行环境密度。

###可视化生态位类别、生态位动态和范围之间的气候类比

ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat90585, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
ecospat.shift.centroids(scores.sp.ecospatcurrent, scores.sp.ecospat90585, scores.clim.ecospatcurrent, scores.clim.ecospat90585)
pdf("Niche Overlap90585.pdf", width = 6, height = 6)
ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat90585, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
ecospat.shift.centroids(scores.sp.ecospatcurrent, scores.sp.ecospat90585, scores.clim.ecospatcurrent, scores.clim.ecospat90585)
dev.off()
#quant=0.25 用于划分边缘气候的环境密度分位数
#选择要绘制的密度：如果interest=1，绘制本机密度，如果interest=2，绘制侵入密度图。
# 计算生态位动态，量化增加、减少和维持不变的部分
niche.dyn <- ecospat.niche.dyn.index(grid.clim.ecospatcurrent, grid.clim.ecospat90585, intersection = 0.1)

# 提取动态信息：expansion (增加), stability (维持不变), unfilling (减少)
expansion <- niche.dyn$dynamic.index.w['expansion']  # 生态位扩展
stability <- niche.dyn$dynamic.index.w['stability']  # 生态位维持不变
unfilling <- niche.dyn$dynamic.index.w['unfilling']  # 生态位未填充

# 打印结果
print(paste("生态位增加 (Expansion):", expansion))
print(paste("生态位维持不变 (Stability):", stability))
print(paste("生态位减少 (Unfilling):", unfilling))

# 保存结果到CSV文件
results <- data.frame(
  Metric = c("Expansion", "Stability", "Unfilling"),
  Value = c(expansion, stability, unfilling)
)

write.csv(results, "D:/eco/Niche_Dynamics_90585.csv", row.names = FALSE)

# 可视化生态位扩展、维持不变和未填充的部分
pdf("Niche_Dynamics_90585.pdf", width = 6, height = 6)
ecospat.plot.niche.dyn(grid.clim.ecospatcurrent, grid.clim.ecospat90585, quant=0.25, interest=1,
                       title= "Niche Overlap", name.axis1="PC1",
                       name.axis2="PC2")
dev.off()








































