
setwd("G:/")
if (!require(tidyverse)) install.packages("tidyverse")
library(biomod2)
library(terra)
library(tidyterra)
library(ggtext)
library(export)
# 加载分布点数据
a_occ <- read.csv("G:/spcies.csv")
head(a_occ)
tail(a_occ)
# 选取研究物种名称
myRespName <- "a"

# 获取存在/不存在数据
myResp <- as.numeric(a_occ[, myRespName])

# 获取经纬度数据
myRespXY <- a_occ[, c("lon", "lat")]

# 加载环境数据
bioclim_ZA_sub <-
  raster::stack(
    c(
      bio4  = "G:/bio4.asc",
      bio5  = "G:/bio5.asc",
      bio12  = "G:bio12.asc",
      bio14  = "G:/bio14.asc",
      bio15 = "G:/bio15.asc",
      water = "G:/elev.asc"
      
    )
  )
myExpl <- rast(bioclim_ZA_sub)

# 数据格式化（选一种即可）
# random生成伪不存在点
myBiomodData.r <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = myRespName,
                                       PA.nb.rep = 2,
                                       PA.nb.absences = 1000,
                                       PA.strategy = 'random')
myBiomodData.r
summary(myBiomodData.r)
plot(myBiomodData.r)

# 运行单一模型

myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData.r,
                                    modeling.id = 'AllModels',
                                    models = c("GLM", "GBM", "GAM",            
                                               "CTA", "ANN", "SRE", "FDA", 
                                               "MARS", "RF",
                                               "MAXNET","MAXENT","XGBOOST"),
                                    CV.strategy = 'random',
                                    CV.nb.rep =1,
                                    CV.perc = 0.8,
                                    OPT.strategy = 'bigboss',
                                    metric.eval = c('TSS','ROC'),
                                    var.import = 2,
                                    seed.val = 10)
myBiomodModelOut

get_evaluations(myBiomodModelOut)

try(write.csv(
  x = get_evaluations(myBiomodModelOut),
  file = paste0(myRespName, "/", myRespName, "_models_evaluations.csv")
))

get_variables_importance(myBiomodModelOut)

try(write.csv(
  x = get_variables_importance(myBiomodModelOut),
  file = paste0(myRespName, "/", myRespName, "_variables_importance.csv")
))



# 各种权重
# 第一个图
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')

p1 <- bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
graph2ppt(p1$plot, file=paste0(myRespName,"_res_plot"), vector.graphic=TRUE, width=9, 
          aspectr=sqrt(2), append = TRUE) 
# 第二个
bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')

p2 <- bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
graph2ppt(p2$plot, file=paste0(myRespName,"_res_plot"), vector.graphic=TRUE, width=9, 
          aspectr=sqrt(2), append = TRUE)
# 第三个
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))

p3 <- bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
graph2ppt(p3$plot, file=paste0(myRespName,"_res_plot"), vector.graphic=TRUE, width=9, 
          aspectr=sqrt(2), append = TRUE)
# 第四个
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))

p4 <- bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
graph2ppt(p4$plot, file=paste0(myRespName,"_res_plot"), vector.graphic=TRUE, width=9, 
          aspectr=sqrt(2), append = TRUE)
# 第五个
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))

p5 <- bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
graph2ppt(p5$plot, file=paste0(myRespName,"_res_plot"), vector.graphic=TRUE, width=9, 
          aspectr=sqrt(2), append = TRUE)
# 第六个
bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))

p6 <- bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))
graph2ppt(p6$plot, file=paste0(myRespName,"_res_plot"), vector.graphic=TRUE, width=9, 
          aspectr=sqrt(2), append = TRUE)


# 响应曲线
# 第七个
mods <- get_built_models(myBiomodModelOut, run = 'RUN1')
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = mods,
                      fixed.var = 'median')

p7 <- bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                            models.chosen = mods,
                            fixed.var = 'median')
graph2ppt(p7$plot, file=paste0(myRespName,"_res_plot"), vector.graphic=TRUE, width=9, 
          aspectr=sqrt(2), append = TRUE)

bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = mods,
                      fixed.var = 'min')
mods <- get_built_models(myBiomodModelOut, full.name = 'a_PA2_RUN1_GLM')
bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
                      models.chosen = mods,
                      fixed.var = 'median',
                      do.bivariate = TRUE)




install.packages("xgboost")
install.packages("biomod2")
library("xgboost")
library("biomod2")
any(is.null(myBiomodData.r))

# 将单一模型投影出来
# file.proj <- paste0(myRespName, "/proj_Current/", myRespName, ".Current.projection.out")
# if (file.exists(file.proj)) {
#   myBiomodProj <- get(load(file.proj))
# } else {
#   myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
#                                     proj.name = 'Currentone',
#                                     new.env = myExpl,
#                                     models.chosen = 'all')
# }
# myBiomodProj
# plot(myBiomodProj)
# 写出单一模型
file.proj <- paste0(myRespName, "/proj_Current/", myRespName, ".Current.projection.out")
if (file.exists(file.proj)) {
  myBiomodProj <- get(load(file.proj))
} else {
  myBiomodProj <- BIOMOD_Projection(
    bm.mod = myBiomodModelOut,
    proj.name = "Currentone",
    new.env = myExpl,
    models.chosen = "all",
    metric.binary = "all",  
    metric.filter = "all",   
    build.clamping.mask = TRUE,   
    output.format = ".tif",
    do.stack = FALSE
  )
}
myBiomodProj
# plot(myBiomodProj)

# 模型集成
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c("EMmean", "EMcv", "EMci", "EMmedian", "EMca", "EMwmean"),
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.8),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      seed.val = 42)
myBiomodEM


get_evaluations(myBiomodEM)
try(write.csv(
  x = get_evaluations(myBiomodEM),
  file = paste0(myRespName, "/", myRespName, "_EMmodels_evaluations.csv")
))
get_variables_importance(myBiomodEM)
try(write.csv(
  x = get_variables_importance(myBiomodEM),
  file = paste0(myRespName, "/", myRespName, "_EMvariables_importance.csv")
))
# 权重可视化，曲线
bm_PlotEvalMean(bm.out = myBiomodEM, dataset = 'validation')
bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'algo', 'merged.by.PA'))
bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('algo', 'expl.var', 'merged.by.PA'))
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM),
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM),
                      fixed.var = 'min')
bm_PlotResponseCurves(bm.out = myBiomodEM, 
                      models.chosen = get_built_models(myBiomodEM, algo = 'EMmean'),
                      fixed.var = 'median',
                      do.bivariate = TRUE)





# 投影集成模型结果
myBiomodEMProj <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  bm.proj = myBiomodProj,
  models.chosen = "all",
  metric.filter = "all",
  output.format = ".tif",
  do.stack = FALSE
)

myBiomodEMProj <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  proj.name = "CurrentEM",
  new.env = myExpl,
  models.chosen = "all",
  metric.filter = "all",
  output.format = ".tif",
  do.stack = FALSE
)


myBiomodEMProj
plot(myBiomodEMProj)




























































