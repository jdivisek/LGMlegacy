######################################################
####      OCCUPANCY OF BIOGEOGRAPHIC REGIONS      ####
######################################################

library(raster)
library(rgdal)
library(sp)
library(xlsx)
library(viridis)
library(landscapemetrics)

reg <- readOGR("shapefile")
plot(reg, col=heat.colors(9)[as.numeric(reg$RegionID)])

###MODERN DISTRIBUTION----------------------------------------------------------
# with multitemporal models
models <- list.files(path= paste("./ME/models_current_Balance/", sep=""), pattern='tif', full.names=T )
names(models) <- spe.names
models[names(models) %in% test.spe.names] <- gsub("models_current_Balance", "models.multitemp_currentavg_Balance", models[names(models) %in% test.spe.names])
models <- models[-21]

range.stats <- list(total.area=list(), max.patch.area=list(), total.area.perc=list(), 
                    max.patch.area.perc=list())

for(q in spe.names[-21])#
{
  print(q)
  
  me.pred <- raster(models[q])
  plot(me.pred)
  plot(reg, add=T)
  
  tab <- sample_lsm(me.pred, y = reg, what = "lsm_p_area", verbose = F, all_classes =T, directions = 8)
  tab$value <- tab$value * 0.01 ##to km2
  
  region.area <- aggregate(value ~ plot_id, data=tab, FUN=sum)$value
  
  ##total area
  ag <- aggregate(value ~ plot_id + class, data=tab, FUN=sum, drop=F) 
  ag$value[is.na(ag$value)] <- 0
  
  range.stats$total.area[[q]] <- ag[ag$class == 1, "value"]
  
  ##max patch ares
  ag <- aggregate(value ~ plot_id + class, data=tab, FUN=max, drop=F) 
  ag$value[is.na(ag$value)] <- 0
  
  range.stats$max.patch.area[[q]] <- ag[ag$class == 1, "value"]
  
  ##% of total area
  range.stats$total.area.perc[[q]] <- round((range.stats$total.area[[q]] / region.area) *100,1)
  
  ##% of max patch
  range.stats$max.patch.area.perc[[q]] <- round((range.stats$max.patch.area[[q]] / region.area) *100,1)
}

range.stats$total.area <- as.data.frame(do.call(rbind, range.stats$total.area))
colnames(range.stats$total.area) <- 1:9
range.stats$max.patch.area <- as.data.frame(do.call(rbind, range.stats$max.patch.area))
colnames(range.stats$max.patch.area) <- 1:9
range.stats$total.area.perc <- as.data.frame(do.call(rbind, range.stats$total.area.perc))
colnames(range.stats$total.area.perc) <- 1:9
range.stats$max.patch.area.perc <- as.data.frame(do.call(rbind, range.stats$max.patch.area.perc))
colnames(range.stats$max.patch.area.perc) <- 1:9

range.stats.withFossils <- range.stats

##write excel file
write.xlsx(range.stats.withFossils$total.area, file="Occupancy-withFossils_modern_2024-02.xlsx", sheetName = "total.area")
write.xlsx(range.stats.withFossils$max.patch.area, file="Occupancy-withFossils_modern_2024-02.xlsx", sheetName = "max.patch.area", append = T)
write.xlsx(range.stats.withFossils$total.area.perc, file="Occupancy-withFossils_modern_2024-02.xlsx", sheetName = "total.area%", append = T)
write.xlsx(range.stats.withFossils$max.patch.area.perc, file="Occupancy-withFossils_modern_2024-02.xlsx", sheetName = "max.patch.area%", append = T)

###LGM DISTRIBUTION-------------------------------------------------------------
# with multitemporal models
models.LGM <- list.files(path= paste("./ME/models_LGMavg_Balance/", sep=""), pattern='tif', full.names=T )
names(models.LGM) <- spe.names
models.LGM[names(models.LGM) %in% test.spe.names] <- gsub("models", "models.multitemp", models.LGM[names(models.LGM) %in% test.spe.names])
models.LGM <- models.LGM[-21]

range.stats <- list(total.area=list(), max.patch.area=list(), total.area.perc=list(), 
                    max.patch.area.perc=list())

for(q in spe.names[-21])#
{
  print(q)
  
  me.pred <- raster(models.LGM[q])
  plot(me.pred)
  plot(reg, add=T)
  
  tab <- sample_lsm(me.pred, y = reg, what = "lsm_p_area", verbose = F, all_classes =T, directions = 8)
  tab$value <- tab$value * 0.01 ##to km2
  
  region.area <- aggregate(value ~ plot_id, data=tab, FUN=sum)$value
  
  ##total area
  ag <- aggregate(value ~ plot_id + class, data=tab, FUN=sum, drop=F) 
  ag$value[is.na(ag$value)] <- 0
  
  range.stats$total.area[[q]] <- ag[ag$class == 1, "value"]
  
  ##max patch ares
  ag <- aggregate(value ~ plot_id + class, data=tab, FUN=max, drop=F) 
  ag$value[is.na(ag$value)] <- 0
  
  range.stats$max.patch.area[[q]] <- ag[ag$class == 1, "value"]
  
  ##% of total area
  range.stats$total.area.perc[[q]] <- round((range.stats$total.area[[q]] / region.area) *100,1)
  
  ##% of max patch
  range.stats$max.patch.area.perc[[q]] <- round((range.stats$max.patch.area[[q]] / region.area) *100,1)
}

range.stats$total.area <- as.data.frame(do.call(rbind, range.stats$total.area))
colnames(range.stats$total.area) <- 1:9
range.stats$max.patch.area <- as.data.frame(do.call(rbind, range.stats$max.patch.area))
colnames(range.stats$max.patch.area) <- 1:9
range.stats$total.area.perc <- as.data.frame(do.call(rbind, range.stats$total.area.perc))
colnames(range.stats$total.area.perc) <- 1:9
range.stats$max.patch.area.perc <- as.data.frame(do.call(rbind, range.stats$max.patch.area.perc))
colnames(range.stats$max.patch.area.perc) <- 1:9

range.stats.LGM.withFossils <- range.stats

##write excel file
write.xlsx(range.stats.LGM.withFossils$total.area, file="Occupancy-withFossils_LGM_2024-02.xlsx", sheetName = "total.area")
write.xlsx(range.stats.LGM.withFossils$max.patch.area, file="Occupancy-withFossils_LGM_2024-02.xlsx", sheetName = "max.patch.area", append = T)
write.xlsx(range.stats.LGM.withFossils$total.area.perc, file="Occupancy-withFossils_LGM_2024-02.xlsx", sheetName = "total.area%", append = T)
write.xlsx(range.stats.LGM.withFossils$max.patch.area.perc, file="Occupancy-withFossils_LGM_2024-02.xlsx", sheetName = "max.patch.area%", append = T)
