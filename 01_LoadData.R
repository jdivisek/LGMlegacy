############################
####      LOAD DATA     ####
############################

library(raster)

##ASSEMBLE CLIMATE DATA---------------------------------------------------------
# get WorldClim grids
wc.paths <- list.files(path="folder", pattern='tif', full.names=TRUE )

# get WorldClim grids
envirem.paths <- list.files(path="folder", pattern='tif', full.names=TRUE)

clim <- stack(c(wc.paths,envirem.paths))
names(clim)
plot(clim[[1:4]])

###IMPORT PRUNED SPECIES OCCURRENCE DATA----------------------------------------
spe <- read.delim("path", header=T)
head(spe)

#merge species
spe[which(spe$SPECIES == "coloradensis" | spe$SPECIES == "cristata" | spe$SPECIES == "pisewensis"), "SPECIES"] <- "cristata agg."
spe[which(spe$SPECIES == "hebes pithodes"), "SPECIES"] <- "hebes"

table(spe$SPECIES)

spe <- spe[order(spe$GENUS, spe$SPECIES), ]

spe.names <- unique(spe$SPECIES)
gen.names <- rep(unique(spe$GENUS), c(4,11,30))

cbind(gen.names, spe.names)

table(spe$SPECIES)

##IMPORT FOSSILS FROM EXTIMATED LGM RANGES--------------------------------------
###European species are represented by the intersection of loess and LGM ranges estimated by Michal
###NA species are represented by the intersection of loess and LGM ranges estimated by Jeff

fossils.range <- read.delim("path", header=T)
head(fossils.range)

test.spe.names <- unique(fossils.range$SPECIES2)

test.gen.names <- rep(unique(fossils.range$GENUS), c(1,7,10))

cbind(test.gen.names, test.spe.names)

table(fossils.range$SPECIES2)

###IMPORT LGM CLIMATE DATA---------------------------------------------------------------
# MPI-ESM simulations
wc.mpi.paths <- list.files(path="folder", pattern='tif', full.names=TRUE )

envirem.mpi.paths <- list.files(path="folder", pattern='tif', full.names=TRUE )

clim.mpi <- stack(c(wc.mpi.paths,envirem.mpi.paths))
plot(clim.mpi[[1:4]])

# CCSM4 simulations
wc.ccsm4.paths <- list.files(path="folder", pattern='tif', full.names=TRUE )

envirem.ccsm4.paths <- list.files(path="folder", pattern='tif', full.names=TRUE )

clim.ccsm4 <- stack(c(wc.ccsm4.paths,envirem.ccsm4.paths))
plot(clim.ccsm4[[1:4]])

# MIROC simulations
wc.miroc.paths <- list.files(path="folder", pattern='tif', full.names=TRUE )

envirem.miroc.paths <- list.files(path="folder", pattern='tif', full.names=TRUE )

clim.miroc <- stack(c(wc.miroc.paths,envirem.miroc.paths))
plot(clim.miroc[[1:4]])



