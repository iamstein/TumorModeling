#for getting the sourced file name whether run by R CMD BATCH, sourced from RStudio, or sourced from another function
get.dirs = function(sys.calls,dirs) {
  ind = grep('.R',sys.calls)
  if (length(ind)==0) {  #function was likely called from R CMD BATCH, need a different way to get name  
    a = commandArgs()
    dirs$Rname = a[grep("\\.R$",a)]
  } else {
    ind = ind[length(ind)]
    dirs$Rname   = basename( as.character(sys.calls[[ind]])[2] )
  }
  dirs$outpref = substr(dirs$Rname,1,7)
  return(dirs)
}