######################################################
####      FILTERING OF SPECIES OCURRENCE DATA     ####
######################################################

### SPECIES SPECIFIC STANDARDIZATION AND FILTERING-----------------------------

##extract climate
spe.clim <- raster::extract(clim, spe[, c("X", "Y")])

#remove points with NAs for climate
spe <- spe[complete.cases(spe.clim),]
spe.clim <- spe.clim[complete.cases(spe.clim),]

###filtering cycle----------------------------------------------------------
spe.filter <- list()
filtering.info <- list(stats = list(), used.axes = list())

for(i in seq(1,length(spe.names)))
{
  print(paste0("********************* ", spe.names[i] ," *********************"))
  
  sel <- spe.clim[spe$SPECIES == spe.names[i],]#
  PCA <- rda(sel, scale=TRUE)
  Naxes <- which((cumsum(eigenvals(PCA))/PCA$tot.chi)*100 > 90)[1]
  # Naxes <- 3
  
  if(nrow(sel) > 80)
  {
    Result <- envDist.filter(env = scores(PCA, display = "sites", choices = 1:Naxes), 
                             coord = spe[spe$SPECIES == spe.names[i], c("X", "Y")], 
                             thresh = 0.1, 
                             min.size = 80, 
                             Seed = 1234,
                             OriginalPoints = FALSE, 
                             RemovedPoints = FALSE,
                             SupplInfo = TRUE)
    
    spe.filter[[i]] <- Result$filtered$coord
    filtering.info$stats[[i]] <- unlist(Result$suppl)
    filtering.info$used.axes[[i]] <- c((cumsum(eigenvals(PCA))/PCA$tot.chi)*100)[1:Naxes]
  }
  else
  {
    spe.filter[[i]] <- spe[spe$SPECIES == spe.names[i], c("X", "Y")]
    filtering.info$stats[[i]] <- c(Original = nrow(sel), Filtered = nrow(sel), OmittedData = 0)
    filtering.info$used.axes[[i]] <- c((cumsum(eigenvals(PCA))/PCA$tot.chi)*100)[1:Naxes]
  }
  
  
}
names(spe.filter) <- spe.names
names(filtering.info$stats) <- spe.names
names(filtering.info$used.axes) <- spe.names

filtering.info$stats
filtering.info$used.axes

###SELECT EQUAL NO OF FOSSILS---------------------------------------------------
fossils.sel <- list()

for(i in seq(1, length(test.spe.names)))
{
  set.seed(1234)
  
  sel <- fossils.range[fossils.range$Species2 == test.spe.names[i], c("POINT_X", "POINT_Y")]
  fossils.sel[[test.spe.names[i]]] <- sel[sample(nrow(sel), size=filtering.info$stats[[test.spe.names[i]]]["Filtered"]),]
  
  if(test.spe.names[i] == "fulvus") ##1/2 from NA, 1/2 from EU
  {
    sel <- fossils.range[fossils.range$Species2 == test.spe.names[i], c("POINT_X", "POINT_Y")]
    
    sel1 <- sel[sel$POINT_X > -2e+06,]
    sel2 <- sel[sel$POINT_X < -2e+06,]
    
    fossils.sel[[test.spe.names[i]]] <- rbind(sel1[sample(nrow(sel1), size=floor(filtering.info$stats[[test.spe.names[i]]]["Filtered"]/2)),],
                                              sel2[sample(nrow(sel2), size=floor(filtering.info$stats[[test.spe.names[i]]]["Filtered"]/2)+1),])
    rm(sel1, sel2)
    
  }
}
names(fossils.sel) <- test.spe.names

lapply(fossils.sel, FUN=nrow)