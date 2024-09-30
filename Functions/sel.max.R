###SELECT NORTHERN RANGE LIMIT

sel.max <- function(x)
{
  if(is.data.frame(x)) {x[which.max(x$Lat),]}
  else{x}
}