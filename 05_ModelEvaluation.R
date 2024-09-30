####################################
####      MODEL EVALUATION      ####
####################################

library(xlsx)
library(abind)

##SINGLE TEMPORAL MODELS------------------------------------------------------

set.seed(1234)
boyce.rp <- randomPoints(mask = conv_hull, n = 50000)

boyce.bg <- raster::extract(clim, boyce.rp)

cv.maxent <- list()

for(q in spe.names)
{
  print(q)
  
  if(nrow(spe.filter[[q]]) < 20){
    k <- 2
  }
  else{
    k <- 5
  }
  
  cv.maxent[[q]] <- cv.me(p.env = raster::extract(clim, spe.filter[[q]]), 
                          a.env = raster::extract(clim, rp10), 
                          p = rep(1, nrow(spe.filter[[q]])), 
                          a = rep(0, nrow(rp10)), k = k, 
                          v = vars[[q]], 
                          boyce.bg = boyce.bg)
  
}

cv.maxent
cv.maxent <- do.call(rbind.data.frame, cv.maxent)

###MULTITEMPORAL MODELS---------------------------------------------------------

set.seed(1234)
boyce.rp <- randomPoints(mask = conv_hull, n = 100000)

###MPI-ESM----------------------------------------------------------------------
boyce.bg <- rbind(raster::extract(clim, boyce.rp[1:50000,]),
                  raster::extract(clim.mpi, boyce.rp[50001:100000,]))

cv.maxent.multi <- list()

for(q in test.spe.names)
{
  print(q)
  
  if(nrow(spe.filter[[q]]) < 20){
    k <- 2
  }
  else{
    k <- 5
  }
  
  cv.maxent.multi[[q]] <- cv.me(p.env = rbind(raster::extract(clim, spe.filter[[q]]),
                                              raster::extract(clim.mpi, fossils.sel[[q]])), 
                                a.env = raster::extract(clim, rp10), 
                                p = rep(1, nrow(spe.filter[[q]])*2), 
                                a = rep(0, nrow(rp10)), k = k, 
                                v = vars.test[[q]], 
                                boyce.bg = boyce.bg)
  
}

cv.maxent.multi
cv.maxent.multi <- do.call(rbind.data.frame, cv.maxent.multi)

###CCSM4------------------------------------------------------------------------
boyce.bg <- rbind(raster::extract(clim, boyce.rp[1:50000,]),
                  raster::extract(clim.ccsm4, boyce.rp[50001:100000,]))

cv.maxent.multi.ccsm4 <- list()

for(q in test.spe.names)
{
  print(q)
  
  if(nrow(spe.filter[[q]]) < 20){
    k <- 2
  }
  else{
    k <- 5
  }
  
  cv.maxent.multi.ccsm4[[q]] <- cv.me(p.env = rbind(raster::extract(clim, spe.filter[[q]]),
                                                    raster::extract(clim.ccsm4, fossils.sel[[q]])), 
                                      a.env = raster::extract(clim, rp10), 
                                      p = rep(1, nrow(spe.filter[[q]])*2), 
                                      a = rep(0, nrow(rp10)), k = k, 
                                      v = vars.test[[q]], 
                                      boyce.bg = boyce.bg)
  
}

cv.maxent.multi.ccsm4
cv.maxent.multi.ccsm4 <- do.call(rbind.data.frame, cv.maxent.multi.ccsm4)

###MIROC------------------------------------------------------------------------
boyce.bg <- rbind(raster::extract(clim, boyce.rp[1:50000,]),
                  raster::extract(clim.miroc, boyce.rp[50001:100000,]))

cv.maxent.multi.miroc <- list()

for(q in test.spe.names)
{
  print(q)
  
  if(nrow(spe.filter[[q]]) < 20){
    k <- 2
  }
  else{
    k <- 5
  }
  
  cv.maxent.multi.miroc[[q]] <- cv.me(p.env = rbind(raster::extract(clim, spe.filter[[q]]),
                                                    raster::extract(clim.miroc, fossils.sel[[q]])), 
                                      a.env = raster::extract(clim, rp10), 
                                      p = rep(1, nrow(spe.filter[[q]])*2), 
                                      a = rep(0, nrow(rp10)), k = k, 
                                      v = vars.test[[q]], 
                                      boyce.bg = boyce.bg)
  
}

cv.maxent.multi.miroc
cv.maxent.multi.miroc <- do.call(rbind.data.frame, cv.maxent.multi.miroc)


###EXPORT CROSS-VALIDATION TABLES-----------------------------------------------

write.xlsx(cv.maxent, file="cross-validation_maxent_2024-02.xlsx", sheetName = "Maxent")

temp_array <- abind(cv.maxent.multi, cv.maxent.multi.ccsm4, cv.maxent.multi.miroc, along=3)
res <- apply(temp_array, 1:2, mean)

write.xlsx(res, file="cross-validation_maxent.multi_2024-02.xlsx", sheetName = "Maxent")

