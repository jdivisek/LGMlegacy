######################################
####      VARIABLE SELECTION      ####
######################################

library(raster)
library(dismo)
library(rJava)
library(usdm)
library(stringr)

##SELECT VARIABLES FOR MODELLING------------------------------------------------
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

# import convex hull for sampling of backgroud points
conv_hull <- raster("path")
plot(conv_hull)

###sample background points
set.seed(1234)
rp10 <- randomPoints(mask = conv_hull, n = 10000)
colnames(rp10) <- c("X", "Y")
points(rp10, cex=0.5, pch=16)

#create constant for variable selection
r <- clim[[1]]
values(r)[!is.na(values(r))] <- 1
plot(r)
names(r) <- "constant"

###VARIABLE SELECTION----------------------------------------------

vars <- list()

for(q in seq(1, length(spe.names)))
{
  print(spe.names[q])
  
  var.imp <- vector("numeric")
  
  for(i in seq(1, length(names(clim))))
  {
    print(names(clim)[i])
    
    env.try <- stack(clim[[i]], r)
    
    mod <- dismo::maxent(x=env.try, p=spe.filter[[q]], a=rp10, removeDuplicates=TRUE,
                         args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1",
                                "replicates=5", "replicatetype=crossvalidate", "threads=2")) #with cross-validation
    
    var.imp[i] <- mod@results[str_detect(rownames(mod@results), "Test.AUC"), ncol(mod@results)]
  }
  
  names(var.imp) <- names(clim)
  var.imp <- sort(var.imp, decreasing = T)
  
  #extract env variables and create correlation table
  env.tab <- as.data.frame(raster::extract(clim, rbind(rp10, spe.filter[[q]])))
  
  #select variables using VIF - best selection method according fossils test
  sel <- var.sel.vif(imp = var.imp, var.tab = env.tab, vif.threshold=10)
  vars[[q]] <- sel[[1]]
  attributes(vars[[q]]) <- list(var.imp=sel[[2]])
  
  #check variable inflation factor
  print(vif(env.tab[, vars[[q]]], maxobservations=20000))
}

names(vars) <- spe.names

###SELECT VARIABLES FOR SPECIS WITH FOSSILS-------------------------------------

vars.test <- list()

for(q in test.spe.names)
{
  print(q)
  
  var.imp <- vector("numeric")
  
  env.tab <- as.data.frame(rbind(raster::extract(clim, spe.filter[[q]]),
                                 raster::extract(clim.mpi, fossils.sel[[q]]),
                                 raster::extract(clim, rp10)))
  
  occ <- c(rep(1, nrow(spe.filter[[q]])+nrow(fossils.sel[[q]])),
           rep(0, nrow(rp10)))
  
  
  for(i in seq(1, ncol(env.tab)))
  {
    print(colnames(env.tab)[i])
    
    env.try <- as.data.frame(cbind(env.tab[, i], rep(1, nrow(env.tab))))
    colnames(env.try) <- c(colnames(env.tab)[i],"constant")
    
    mod <- dismo::maxent(x=env.try, p=occ, removeDuplicates=TRUE,
                         args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1",
                                "replicates=5", "replicatetype=crossvalidate", "threads=2")) #with crossvalidation
    
    var.imp[i] <- mod@results[str_detect(rownames(mod@results), "Test.AUC"), ncol(mod@results)]
  }
  
  names(var.imp) <- colnames(env.tab)
  var.imp <- sort(var.imp, decreasing = T)
  
  sel <- var.sel.vif(imp = var.imp, var.tab = env.tab, vif.threshold=10)
  vars.test[[q]] <- sel[[1]] 
  attributes(vars.test[[q]]) <- list(var.imp=sel[[2]])
  
  #check variable inflation factor
  print(vif(env.tab[, vars.test[[q]]], maxobservations=20000))
}
names(vars.test) <- test.spe.names

