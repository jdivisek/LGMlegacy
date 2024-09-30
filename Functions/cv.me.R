###CROSS-VALIDATION FOR MAXENT MODELS

cv.me <-  function(p.env, a.env, p, a, k = 5, v, boyce.bg)
{
  require(dismo)
  require(raster)
  require(ecospat)
  
  p <- as.data.frame(cbind(p, p.env[, v]))
  a <- as.data.frame(cbind(a, a.env[, v]))
  boyce.bg <- boyce.bg[,v]
  
  ##create bins
  set.seed(12)
  p.folds <- kfold(p, k=k)
  
  auc.train <- vector("numeric")
  auc.test <- vector("numeric")
  overfit <- vector("numeric")
  tss.train <- vector("numeric")
  tss.test <- vector("numeric")
  boyce.train <- vector("numeric")
  boyce.test <- vector("numeric")
  
  for(i in 1:k)
  {
    print(paste("Cross-validation:", i ,"out of", k, sep= " "))
    
    train.p <- as.data.frame(p[p.folds != i, ]) ##select training presences from modern occurence data
    test.p <- as.data.frame(p[p.folds == i, ]) ##select testing presences
    
    mod <- dismo::maxent(x = rbind(train.p[,-1], a[,-1]), p = c(train.p[,1], a[,1]), 
                         removeDuplicates = TRUE, args=c("addsamplestobackground=true", "autofeature=true", "betamultiplier=1",
                                                         "replicates=1", "threads=1"))
    
    e.train <- evaluate(p = train.p[,-1], a = a[,-1], model = mod, tr = seq(0.001, 1, 0.001))
    e.test <- evaluate(p = test.p[,-1], a = a[,-1], model = mod, tr = seq(0.001, 1, 0.001))
    
    auc.train[i] <- e.train@auc
    auc.test[i] <- e.test@auc
    overfit[i] <- e.train@auc - e.test@auc
    tss.train[i] <- max(e.train@TPR - e.train@FPR)
    tss.test[i] <- max(e.test@TPR - e.test@FPR)
    
    pr <- dismo::predict(mod, boyce.bg)
    pr <- pr[!is.na(pr)]
    boyce.train[i] <- ecospat.boyce(pr, dismo::predict(mod, train.p[,-1]),
                                    PEplot=F, res=100)$cor
    boyce.test[i] <- ecospat.boyce(pr, dismo::predict(mod, test.p[,-1]),
                                   PEplot=F, res=100)$cor
    
  }
  
  res <-   list(auc.train = mean(auc.train),
                auc.train.sd = sd(auc.train),
                auc.test = mean(auc.test),
                auc.test.sd = sd(auc.test),
                overfit = mean(overfit),
                overfit.sd = sd(overfit),
                tss.train = mean(tss.train),
                tss.train.sd = sd(tss.train),
                tss.test = mean(tss.test),
                tss.test.sd = sd(tss.test),
                boyce.train = mean(boyce.train),
                boyce.train.sd = sd(boyce.train),
                boyce.test = mean(boyce.test),
                boyce.test.sd = sd(boyce.test))
  
  return(do.call(cbind.data.frame, res))
  
}
