#get filename of current script being executed
#input should be sys.calls()
get.dirs = function(sys.calls,dirs) {
  if (is.null(sys.calls)) {
    dirs$Rname = '_';
  } else {
      sys.calls = sys.calls[!grepl("\n",sys.calls)] #don't fully understand why this is sometimes needed, but it is
      ind = grep('.R',sys.calls)
    if (length(ind)==0) {  #function was likely called from R CMD BATCH, need a different way to get name because sys.calls won't work
      a = commandArgs()
      dirs$Rname = a[grep("\\.R$",a)]
    } else {
      ind = ind[length(ind)]
      dirs$Rname   = basename( as.character(sys.calls[[ind]])[2] )
    }
  }
  dirs$outpref = substr(dirs$Rname,1,gregexpr("_",dirs$Rname)[[1]][1])
  dirs$outpref.full = paste0(dirs$top,dirs$results,dirs$outpref)
  if (grepl("_\\d\\d\\d.R$",dirs$Rname)) #when I call as specific monolix/nonmem run, I end the script with _###.R.  
    dirs$run = paste0("run",substr(dirs$Rname,nchar(dirs$Rname)-4,nchar(dirs$Rname)-2))
  else
    dirs$run = NA
  return(dirs)
}