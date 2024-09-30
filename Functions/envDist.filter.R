####ENVIRONMENTAL FILTERING OF SPECIES OCURRENCE RECORDS

envDist.filter <- function(env, coord, thresh, min.size = 80, Seed, OriginalPoints = TRUE, RemovedPoints = TRUE, SupplInfo = TRUE)
{
  if(nrow(env) != nrow(coord))
  {
    stop("Error: 'env' and 'coord' must have the same number of rows")
  }
  
  set.seed(Seed)
  rand <- sample(nrow(coord))
  
  coord <- as.data.frame(coord[rand, ])
  env <- as.data.frame(env[rand, ])
  
  env.dist <- dist(env, method="euclidean") #calculate environmental distance matrix
  
  env.dist <- as.matrix(env.dist) ##make matrix of distances
  diag(env.dist) <- NA
  d <- colMeans(env.dist, na.rm = TRUE) #calculate average distnace of points from each other point
  
  env.dist <- env.dist[order(d), order(d)] #order distance matrix
  env <- env[order(d), ]
  coord <- coord[order(d), ]
  
  original.env <- env
  original.coord <- coord
  
  if(RemovedPoints)
  {
    removed.coord <- list()
    removed.env <- list()
  }
  
  step <- 1
  repeat
  {
    if(any(env.dist < thresh & dim(env.dist)[1] > min.size, na.rm = T) == FALSE) #test criterion
    {
      break
    }
    
    most.similar <- apply(env.dist, 2, FUN=min, na.rm=TRUE) #find most similar pairs
    to.rm <- which.min(most.similar)
    
    if(RemovedPoints)
    {
      removed.coord[[i]] <- coord[to.rm,]
      removed.env[[i]] <- env[to.rm,]
    }
    
    coord <- coord[-c(to.rm),] #update coordinates
    env <- env[-c(to.rm),] #update envi. data
    
    env.dist <- env.dist[-c(to.rm), -c(to.rm)] #update environmental matrix
    print(paste0(step, ". plot removed"))
    step <- step+1
  }
  
  Result <- list()
  Result$filtered <- list(coord = coord, env = env)
  
  if(OriginalPoints)
  {
    Result$original <- list(coord = original.coord, env = original.env)
  }
  
  if(RemovedPoints)
  {
    Result$removed <- list(coord = do.call(rbind.data.frame, removed.coord),
                           env = do.call(rbind.data.frame, removed.env))
  }
  
  if(SupplInfo)
  {
    Result$suppl <- list(Original = nrow(original.coord), 
                         Filtered = nrow(coord),
                         OmittedData = nrow(original.coord) - nrow(coord))
  }
  
  return(Result)
}
