########################################################
####      DISPLACEMENT OF MODERN AND LGM RANGES     ####
########################################################

library(raster)
library(rgdal)
library(reshape2)
library(xlsx)
library(landscapemetrics)
library(geosphere)

grat.p <- read.delim("points_along_meridians", header=T, row.names = 1)
head(grat.p)

###With multitemporal models
models <- list.files(path= paste("./ME/models_current_Balance/", sep=""), pattern='tif', full.names=T )
names(models) <- spe.names
models[names(models) %in% test.spe.names] <- gsub("models_current_Balance", "models.multitemp_currentavg_Balance", models[names(models) %in% test.spe.names])
models <- models[-21]

models.LGM <- list.files(path= paste("./ME/models_LGMavg_Balance/", sep=""), pattern='tif', full.names=T )
names(models.LGM) <- spe.names
models.LGM[names(models.LGM) %in% test.spe.names] <- gsub("models", "models.multitemp", models.LGM[names(models.LGM) %in% test.spe.names])
models.LGM <- models.LGM[-21]

range45.dists <- list()

for(q in spe.names[-21])
{
  print(q)
  
  mod <- raster(models[q])
  mod.lgm <- raster(models.LGM[q])
  
  mod <- getPatch(mod, min.size = 763) #763 ~ 45000 km2
  mod <- sum(mod)
  
  mod.lgm <- getPatch(mod.lgm, min.size = 763) #763 ~ 45000 km2
  mod.lgm <- sum(mod.lgm)
  
  # plot(mod)
  # plot(mod.lgm)
  
  p.pres <- extract(mod, grat.p[, c("X", "Y")])
  p.pres[is.na(p.pres)] <- 0
  p.lgm <- extract(mod.lgm, grat.p[, c("X", "Y")])
  p.lgm[is.na(p.lgm)] <- 0
  
  p.pres <- grat.p[p.pres == 1,]
  p.pres$period <- "present"
  p.lgm <- grat.p[p.lgm == 1,]
  p.lgm$period <- "lgm"
  
  p.pres <- split(p.pres, as.factor(p.pres$grid_code))
  p.lgm <- split(p.lgm, as.factor(p.lgm$grid_code))
  
  p.pres <- do.call(rbind.data.frame, lapply(p.pres, FUN=sel.min))
  p.lgm <- do.call(rbind.data.frame, lapply(p.lgm, FUN=sel.max))
  
  # points(p.lgm[, c("X", "Y")], pch=16, cex=0.5, col="red")
  # points(p.pres[, c("X", "Y")], pch=16, cex=0.5, col="blue")
  
  p <- rbind(p.pres, p.lgm)
  # head(p)
  
  p <- split(p, as.factor(p$grid_code))
  # p <- lapply(p, FUN=function(x){ifelse(nrow(x)>1, return(x), return(NULL))})
  
  ##calculate distnaces in km
  p <- lapply(p, FUN=function(x){if(nrow(x)>1){x$Distance <- distGeo(x[, c("Lon", "Lat")])[1]/1000; return(x)}})
  
  p <- do.call(rbind.data.frame, p)
  # head(p)
  
  p <- p[order(p$grid_code, p$period), ]
  
  p.pres <- p[p$period == "present",]
  p.lgm <- p[p$period == "lgm",]
  # dim(p.pres)
  # dim(p.lgm)
  # head(p.pres)
  
  ##positive difference = the southern edge of the modern range is farther north than the northern edge of the LGM range = the ranges don't overlap
  ##negative difference = negative distance = the southern edge of the modern range is further south than the northern edge of the LGM range = the ranges overlap
  
  res <- p.pres[,c(1:7, 9)]
  res$Diff <- p.pres$Lat - p.lgm$Lat
  # head(res)
  
  res$Distance[res$Diff < 0] <- res$Distance[res$Diff < 0]*-1
  
  res$Sp <- q
  
  boxplot(Distance ~ Name, data=res, name=q)
  title(q)
  
  range45.dists[[q]] <- res
}

names(range45.dists)

range45.dists.mean1 <- lapply(range45.dists, FUN= function(x){aggregate(x$Distance, by=list(x$Name, x$Sp), mean)})
range45.dists.mean1 <- do.call(rbind.data.frame, range45.dists.mean1)
range45.dists.mean1 <- dcast(range45.dists.mean1, Group.2 ~ Group.1, value.var="x", fill=NA)
rownames(range45.dists.mean1) <- range45.dists.mean1$Group.2
range45.dists.mean1 <- range45.dists.mean1[,-1]

range45.dists.sd1 <- lapply(range45.dists, FUN= function(x){aggregate(x$Distance, by=list(x$Name, x$Sp), sd)})
range45.dists.sd1 <- do.call(rbind.data.frame, range45.dists.sd1)
range45.dists.sd1 <- dcast(range45.dists.sd1, Group.2 ~ Group.1, value.var="x", fill=NA)
rownames(range45.dists.sd1) <- range45.dists.sd1$Group.2
range45.dists.sd1 <- range45.dists.sd1[,-1]

range45.dists.mean2 <- lapply(range45.dists, FUN= function(x){aggregate(x$Distance, by=list(x$RegionID, x$Sp), mean)})
range45.dists.mean2 <- do.call(rbind.data.frame, range45.dists.mean2)
range45.dists.mean2 <- dcast(range45.dists.mean2, Group.2 ~ Group.1, value.var="x", fill=NA)
rownames(range45.dists.mean2) <- range45.dists.mean2$Group.2
range45.dists.mean2 <- range45.dists.mean2[,-1]

range45.dists.sd2 <- lapply(range45.dists, FUN= function(x){aggregate(x$Distance, by=list(x$RegionID, x$Sp), sd)})
range45.dists.sd2 <- do.call(rbind.data.frame, range45.dists.sd2)
range45.dists.sd2 <- dcast(range45.dists.sd2, Group.2 ~ Group.1, value.var="x", fill=NA)
rownames(range45.dists.sd2) <- range45.dists.sd2$Group.2
range45.dists.sd2 <- range45.dists.sd2[,-1]

write.xlsx(range45.dists.mean1, file="Range_distances_withFossils_2024-03.xlsx", sheetName = "BigRegions_mean_km")
write.xlsx(range45.dists.sd1, file="Range_distances_withFossils_2024-03.xlsx", sheetName = "BigRegions_sd_km", append = T)
write.xlsx(range45.dists.mean2, file="Range_distances_withFossils_2024-03.xlsx", sheetName = "SmallRegions_mean_km", append = T)
write.xlsx(range45.dists.sd2, file="Range_distances_withFossils_2024-03.xlsx", sheetName = "SmallRegions_sd_km", append = T)

