####################################################
####      OVERLAP OF MODERN AND LGM RANGES      ####
####################################################

library(raster)
library(landscapemetrics)
library(xlsx)
library(rgdal)
library(reshape2)

reg <- readOGR("shapefile")
plot(reg, col=heat.colors(9)[as.numeric(reg$RegionID)])

### With multitemporal models
models <- list.files(path= paste("./ME/models_current_Balance/", sep=""), pattern='tif', full.names=T )
names(models) <- spe.names
models[names(models) %in% test.spe.names] <- gsub("models_current_Balance", "models.multitemp_currentavg_Balance", models[names(models) %in% test.spe.names])
models <- models[-21]

models.LGM <- list.files(path= paste("./ME/models_LGMavg_Balance/", sep=""), pattern='tif', full.names=T )
names(models.LGM) <- spe.names
models.LGM[names(models.LGM) %in% test.spe.names] <- gsub("models", "models.multitemp", models.LGM[names(models.LGM) %in% test.spe.names])
models.LGM <- models.LGM[-21]

tab <- list()

for(q in spe.names[-21])
{
  print(q)
  
  mod <- raster(models[q])
  mod.lgm <- raster(models.LGM[q])
  
  mod[mod == 1] <- 2 ##assign 2 to modern range
  mod[is.na(mod)] <- 0 ##assign 0 to NA values (LGM range currently below seal level would be lost otherwise)
  # plot(mod)
  
  overlap <- mod.lgm+mod
  plot(overlap, main=q, maxpixels=3000000)
  #0 = LGM landscape
  #1 = only LGM range
  #2= only modern range
  #3= overlap of LGM and modern ranges
  
  stat <- sample_lsm(overlap, y = reg, what = "lsm_c_ca", verbose = F, all_classes =T, directions = 8)
  stat$value <- stat$value * 0.01 ##to km2
  
  tab[[q]] <- aggregate(value ~ class + plot_id, data=stat, FUN = sum)
  tab[[q]]$species <- q
}

tab <- do.call(rbind.data.frame, tab)
head(tab)

range.overlap <- list()

range.overlap$LGM.only <- dcast(tab[tab$class == "1", ], species ~ plot_id, fun.aggregate=sum, fill=0, value.var="value")
range.overlap$Modern.only <- dcast(tab[tab$class == "2", ], species ~ plot_id, fun.aggregate=sum, fill=0, value.var="value")
range.overlap$Overlap <- dcast(tab[tab$class == "3", ], species ~ plot_id, fun.aggregate=sum, fill=0, value.var="value")

rownames(range.overlap$Overlap) <- range.overlap$Overlap[,1]
range.overlap$Overlap <- range.overlap$Overlap[spe.names[-21],]

rownames(range.overlap$LGM.only) <- range.overlap$LGM.only[,1]
range.overlap$LGM.only <- range.overlap$LGM.only[spe.names[-21],]

rownames(range.overlap$Modern.only) <- range.overlap$Modern.only[,1]
range.overlap$Modern.only <- range.overlap$Modern.only[spe.names[-21],]

write.xlsx(range.overlap$Overlap, file="Range_overlap-withFossils_2024-02.xlsx", sheetName = "Overlap")
write.xlsx(range.overlap$LGM.only, file="Range_overlap-withFossils_2024-02.xlsx", sheetName = "LGM.only", append = T)
write.xlsx(range.overlap$Modern.only, file="Range_overlap-withFossils_2024-02.xlsx", sheetName = "Modern.only", append = T)
