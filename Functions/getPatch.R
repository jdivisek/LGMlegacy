##FIND SPATIALLY CONTINUOUS PATCHES OF THE RANGE

getPatch <- function(pa, min.size, directions=8)
{
  require(raster)
  
  if(any(values(pa) > 0, na.rm = T))
  {
    #find spatial clusters
    refs <- raster::clump(pa, directions=directions, gaps=FALSE)
    
    #get No. cells per cluster
    a <- table(values(refs))
    
    for(i in 1:length(a)) {if(a[i] < min.size){refs[refs == i] <- NA}}
    
    ##assign new values to single refugia
    a <- unique(values(refs))[-1]
    
    if(length(a) > 0)
    {
      refs <- reclassify(refs, rcl = matrix(c(a, seq(1, length(a))), ncol=2, byrow = F))
      
      ##assign 0 to NA values
      refs[is.na(refs)] <- 0
      
      #make raterstack with single refugia (origins)
      rs <- stack()
      
      for(i in unique(values(refs))[-1])
      {
        rs <- addLayer(rs, refs)
        
        rs[[i]][rs[[i]] != i] <- 0
        rs[[i]][rs[[i]] == i] <- 1
      }
      names(rs) <- paste0("Ref", 1:maxValue(refs))
      
      return(rs)
    }
    else
    {
      return(NULL)
    } 
  }
  else
  {
    return(NULL)
  }
}