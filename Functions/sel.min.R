###SELECT SOUTHERN RANGE LIMIT

sel.min <- function(x)
{
  if(is.data.frame(x)) {x[which.min(x$Lat),]}
  else{x}
}