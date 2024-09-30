###SELECTION OF UNCORRELATED VARIABLES USING VIF

var.sel.vif <- function(imp, var.tab, vif.threshold=10, maxobservations=20000)
{
  require(usdm)
  
  var.tab <- var.tab[, names(imp)]#sort table
  
  sel <- 1
  
  for(i in seq(2, ncol(var.tab)))
  {
    test <- vif(var.tab[, c(sel, i)], maxobservations=maxobservations)$VIF
    
    if(sum(test > vif.threshold) == 0)
    {
      sel <- c(sel, i)
    }
  }
  
  return(list(colnames(var.tab)[sel], sel))
}