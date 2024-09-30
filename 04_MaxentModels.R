################################
####      MAXENT MODELS     ####
################################

library(raster)
library(dismo)
library(rJava)
library(usdm)

jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

##UNITEMPORAL MAXENT MODELS-----------------------------------------------------------------

me <- list()

for(i in spe.names)
{
  print(i)
  
  me[[i]] <- dismo::maxent(x=clim[[vars[[i]]]], p=spe.filter[[i]], a=rp10, removeDuplicates=TRUE,
                           args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1",
                                  "replicates=1", "replicatetype=crossvalidate", "threads=2"))
}


##EXPORT MAXENT SUITABILITY MAPS------------------------------------------------

for(i in 1:length(spe.names))
{
  print(paste0("********************* ", spe.names[i], " *********************"))
  ###For current condition###
  print("Current conditions")
  
  mod <- dismo::predict(object = me[[spe.names[i]]], x = clim)
  
  plot(mod, main=paste0(spe.names[i], " | present"))
  writeRaster(mod, filename = paste0("./ME/models_current/", gen.names[i], "_", spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  ###For LGM condition###
  print("LGM conditions")
  
  mod <- dismo::predict(object = me[[spe.names[i]]], x = clim.mpi)
  
  plot(mod, main=paste0(spe.names[i], " | LGM - MPI"))
  writeRaster(mod, filename = paste0("./ME/models_LGM/", gen.names[i], "_", spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- dismo::predict(object = me[[spe.names[i]]], x = clim.ccsm4)
  
  plot(mod, main=paste0(spe.names[i], " | LGM - CCSM4"))
  writeRaster(mod, filename = paste0("./ME/models_LGM_CCSM4/", gen.names[i], "_", spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- dismo::predict(object = me[[spe.names[i]]], x = clim.miroc)
  
  plot(mod, main=paste0(spe.names[i], " | LGM - MIROC"))
  writeRaster(mod, filename = paste0("./ME/models_LGM_MIROC/", gen.names[i], "_", spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
}

###MULTITEMPORAL MODELS---------------------------------------------------------

me.multi <- list()
me.multi.ccsm4 <- list()
me.multi.miroc <- list()

for(q in test.spe.names)
{
  print(q)
  
  occ <- c(rep(1, nrow(spe.filter[[q]])+nrow(fossils.sel[[q]])),
           rep(0, nrow(rp10)))
  
  ##MPI
  env.tab <- as.data.frame(rbind(extract(clim, spe.filter[[q]]),
                                 extract(clim.mpi, fossils.sel[[q]]),
                                 extract(clim, rp10)))
  me.multi[[q]]   <- dismo::maxent(x=env.tab[, vars.test[[q]]], p=occ, removeDuplicates=TRUE,
                                   args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1",
                                          "replicates=1", "replicatetype=crossvalidate", "threads=2"))
  ##CCSM4
  env.tab <- as.data.frame(rbind(extract(clim, spe.filter[[q]]),
                                 extract(clim.ccsm4, fossils.sel[[q]]),
                                 extract(clim, rp10)))
  me.multi.ccsm4[[q]]   <- dismo::maxent(x=env.tab[, vars.test[[q]]], p=occ, removeDuplicates=TRUE,
                                         args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1",
                                                "replicates=1", "replicatetype=crossvalidate", "threads=2"))
  
  ##MIROC
  env.tab <- as.data.frame(rbind(extract(clim, spe.filter[[q]]),
                                 extract(clim.miroc, fossils.sel[[q]]),
                                 extract(clim, rp10)))
  me.multi.miroc[[q]]   <- dismo::maxent(x=env.tab[, vars.test[[q]]], p=occ, removeDuplicates=TRUE,
                                         args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1",
                                                "replicates=1", "replicatetype=crossvalidate", "threads=2"))
  
}

###EXPORT MULTITEMPORAL MODELS--------------------------------------------------

for(i in 1:length(test.spe.names))
{
  print(paste0("********************* ", test.spe.names[i], " *********************"))
  ###For current condition###
  print("Current conditions")
  
  mod <- dismo::predict(object = me.multi[[test.spe.names[i]]], x = clim)

  plot(mod, main=paste0(test.spe.names[i], " | present"))
  writeRaster(mod, filename = paste0("./ME/models.multitemp_current/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- dismo::predict(object = me.multi.ccsm4[[test.spe.names[i]]], x = clim)
  
  plot(mod, main=paste0(test.spe.names[i], " | present"))
  writeRaster(mod, filename = paste0("./ME/models.multitemp_current_CCSM4/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- dismo::predict(object = me.multi.miroc[[test.spe.names[i]]], x = clim)
  
  plot(mod, main=paste0(test.spe.names[i], " | present"))
  writeRaster(mod, filename = paste0("./ME/models.multitemp_current_MIROC/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  ###For LGM conditions###
  print("LGM conditions")
  
  mod <- dismo::predict(object = me.multi[[test.spe.names[i]]], x = clim.mpi)
  
  plot(mod, main=paste0(test.spe.names[i], " | LGM - MPI"))
  writeRaster(mod, filename = paste0("./ME/models.multitemp_LGM/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- dismo::predict(object = me.multi.ccsm4[[test.spe.names[i]]], x = clim.ccsm4)
  
  plot(mod, main=paste0(test.spe.names[i], " | LGM - CCSM4"))
  writeRaster(mod, filename = paste0("./ME/models.multitemp_LGM_CCSM4/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- dismo::predict(object = me.multi.miroc[[test.spe.names[i]]], x = clim.miroc)
  
  plot(mod, main=paste0(test.spe.names[i], " | LGM - MIROC"))
  writeRaster(mod, filename = paste0("./ME/models.multitemp_LGM_MIROC/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
}

###MODERN SPECIES RANGES--------------------------------------------------------
##single-temporal models--------------------------------------------------------
models <- list.files("./ME/models_current/", full.names = T)
names(models) <- spe.names

for(i in 1:length(spe.names))
{
  q <- spe.names[i]
  print(q)
  
  mod <- raster(models[q])
  plot(mod, main = q)
  
  mod <- mod * frame30 #cut at 30°N
  
  mod <- mod * lake_current #remove lakes
  mod <- mod * glacier_current #remove areas under ice sheet
  
  thresh <- me[[q]]@results["Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold", 1]
  
  mod[mod >= thresh] <- 1
  mod[mod != 1] <- 0
  plot(mod, main = q)
  
  writeRaster(mod, filename = paste0("./ME/models_current_Balance/", gen.names[i], "_", spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
}

##Multitemporal models----------------------------------------------------------
models <- list.files("./ME/models.multitemp_current/", full.names = T)
names(models) <- test.spe.names

models.CCSM4 <- list.files("./ME/models.multitemp_current_CCSM4/", full.names = T)
names(models.CCSM4) <- test.spe.names

models.MIROC <- list.files("./ME/models.multitemp_current_MIROC/", full.names = T)
names(models.MIROC) <- test.spe.names

for(i in 1:length(test.spe.names))
{
  q <- test.spe.names[i]
  print(q)
  
  rs <- stack(c(models[q], models.CCSM4[q], models.MIROC[q]))
  mod <- mean(rs)
  plot(mod, main = q)
  
  writeRaster(mod, filename = paste0("./ME/models.multitemp_currentavg/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- mod * frame30 #cut at 30°N
  
  mod <- mod * lake_current #remove lakes
  mod <- mod * glacier_current #remove areas under ice sheet
  
  thresh <- mean(c(me.multi[[q]]@results["Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold", 1],
                   me.multi.ccsm4[[q]]@results["Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold", 1],
                   me.multi.miroc[[q]]@results["Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold", 1]))
  
  mod[mod >= thresh] <- 1
  mod[mod != 1] <- 0
  plot(mod, main = q)
  
  writeRaster(mod, filename = paste0("./ME/models.multitemp_currentavg_Balance/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
}

###AVERAGE LGM PREDICTIONS------------------------------------------------------
##Unitemporal models------------------------------------------------------------
models <- list.files("./ME/models_LGM/", full.names = T)
names(models) <- spe.names

models.CCSM4 <- list.files("./ME/models_LGM_CCSM4/", full.names = T)
names(models.CCSM4) <- spe.names

models.MIROC <- list.files("./ME/models_LGM_MIROC/", full.names = T)
names(models.MIROC) <- spe.names

for(i in 1:length(spe.names))
{
  q <- spe.names[i]
  print(q)
  
  rs <- stack(c(models[q], models.CCSM4[q], models.MIROC[q]))
  mod <- mean(rs)
  plot(mod, main = q)
  
  writeRaster(mod, filename = paste0("./ME/models_LGMavg/", gen.names[i], "_", spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- mod * frame30 #cut at 30°N
  
  mod <- mod * lake_lgm #remove lakes
  mod <- mod * glacier_lgm #remove areas under ice sheet
  
  thresh <- me[[q]]@results["Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold", 1]
  
  mod[mod >= thresh] <- 1
  mod[mod != 1] <- 0
  plot(mod, main = q)
  
  writeRaster(mod, filename = paste0("./ME/models_LGMavg_Balance/", gen.names[i], "_", spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
}

##multi-temporal models---------------------------------------------------------
models <- list.files("./ME/models.multitemp_LGM/", full.names = T)
names(models) <- test.spe.names

models.CCSM4 <- list.files("./ME/models.multitemp_LGM_CCSM4/", full.names = T)
names(models.CCSM4) <- test.spe.names

models.MIROC <- list.files("./ME/models.multitemp_LGM_MIROC/", full.names = T)
names(models.MIROC) <- test.spe.names

for(i in 1:length(test.spe.names))
{
  q <- test.spe.names[i]
  print(q)
  
  rs <- stack(c(models[q], models.CCSM4[q], models.MIROC[q]))
  mod <- mean(rs)
  plot(mod, main = q)
  
  writeRaster(mod, filename = paste0("./ME/models.multitemp_LGMavg/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
  mod <- mod * frame30 #cut at 30°N
  
  mod <- mod * lake_lgm #remove lakes
  mod <- mod * glacier_lgm #remove areas under ice sheet
  
  thresh <- mean(c(me.multi[[q]]@results["Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold", 1],
                   me.multi.ccsm4[[q]]@results["Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold", 1],
                   me.multi.miroc[[q]]@results["Balance.training.omission..predicted.area.and.threshold.value.Cloglog.threshold", 1]))
  
  mod[mod >= thresh] <- 1
  mod[mod != 1] <- 0
  plot(mod, main = q)
  
  writeRaster(mod, filename = paste0("./ME/models.multitemp_LGMavg_Balance/", test.gen.names[i], "_", test.spe.names[i], ".tif", sep=""), format="GTiff", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6"))
  
}
