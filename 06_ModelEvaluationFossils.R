##################################################
####      MODEL EVALUATION USING FOSSILS      ####
##################################################

###Boyce index for randomly sampled fossil records

##Multitemporal models----------------------------------------------------------
models.LGM <- list.files(path= paste("./ME/models.multitemp_LGMavg/", sep=""), pattern='tif', full.names=TRUE )
names(models.LGM) <- test.spe.names

set.seed(1234)
boyce.bg <- randomPoints(clim[[1]], 100000) ##for boyce100.ME.tab.multi

boyce100.ME.tab.multi <- list()

for(q in seq(1, length(test.spe.names)))
{
  print(test.spe.names[q])
  
  boyce100.ME.tab.multi[[q]] <- vector("numeric")
  
  if(test.spe.names[q] %in% test.spe.names[c(1,3,4,10,12,13,14,18)])
  {
    f <- fossils.range[fossils.range$Species2 == test.spe.names[q], c("POINT_X", "POINT_Y")]
  }
  else
  {
    f <- fossils[fossils$SPECIES2 == test.spe.names[q], c("X", "Y")]
  }
  
  me.pred <- raster(models.LGM[test.spe.names[q]])
  fit <- raster::extract(me.pred, boyce.bg)
  
  for(i in seq(1, length(seeds)))
  {
    print(i)
    
    set.seed(seeds[i])
    sel <- f[sample(nrow(f), size=100),]
    obs <- raster::extract(me.pred, sel)
    
    boyce100.ME.tab.multi[[q]][i] <- ecospat::ecospat.boyce(fit=c(fit[!is.na(fit)], obs[!is.na(obs)]), obs=obs[!is.na(obs)], PEplot=F, res=100)$cor
  }
}
names(boyce100.ME.tab.multi) <- test.spe.names

boyce100.ME.tab.multi <- do.call(rbind, boyce100.ME.tab.multi)
rowMeans(boyce100.ME.tab.multi)

write.table(cbind(rowMeans(boyce100.ME.tab.multi, na.rm=T),
                  apply(boyce100.ME.tab.multi, 1, sd, na.rm=T)), file="boyce100.ME.tab.multi_2024-02.txt", sep="\t", dec=".")


##Unitemporal models------------------------------------------------------------

models.LGM <- list.files(path= paste("./ME/models_LGMavg/", sep=""), pattern='tif', full.names=TRUE )
names(models.LGM) <- spe.names
models.LGM <- models.LGM[test.spe.names]

set.seed(1234)
boyce.bg <- randomPoints(clim.mpi[[1]], 100000)

boyce100.ME.tab <- list()

for(q in seq(1, length(test.spe.names)))
{
  print(test.spe.names[q])
  
  boyce100.ME.tab[[q]] <- vector("numeric")
  
  if(test.spe.names[q] %in% test.spe.names[c(1,3,4,10,12,13,14,18)])
  {
    f <- fossils.range[fossils.range$Species2 == test.spe.names[q], c("POINT_X", "POINT_Y")]
  }
  else
  {
    f <- fossils[fossils$SPECIES2 == test.spe.names[q], c("X", "Y")]
  }
  
  me.pred <- raster(models.LGM[test.spe.names[q]])
  fit <- raster::extract(me.pred, boyce.bg)
  
  for(i in seq(1, length(seeds)))
  {
    print(i)
    
    set.seed(seeds[i])
    sel <- f[sample(nrow(f), size=100),]
    obs <- raster::extract(me.pred, sel)
    
    boyce100.ME.tab[[q]][i] <- ecospat::ecospat.boyce(fit=c(fit[!is.na(fit)], obs[!is.na(obs)]), obs=obs[!is.na(obs)], PEplot=F, res=100)$cor
  }
}
names(boyce100.ME.tab) <- test.spe.names

boyce100.ME.tab <- do.call(rbind, boyce100.ME.tab)
rowMeans(boyce100.ME.tab, na.rm=T)

write.table(cbind(rowMeans(boyce100.ME.tab, na.rm=T),
                  apply(boyce100.ME.tab, 1, sd, na.rm=T)), file="boyce100.ME.tab_2024-02.txt", sep="\t", dec=".")

cbind(rowMeans(boyce100.ME.tab, na.rm=T),
      rowMeans(boyce100.ME.tab.multi, na.rm=T))

