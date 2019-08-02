#Use to find possible MS2 from a RawSpec File. Returns a table of scans with file names.
# Rt is a range in seconds

LOB_findMS2 <- function(rawSpec,mz,rt,ppm){

  mzrange <- mz*(0.000001*ppm)
  mzlow <- (mz-mzrange)
  mzhigh <- (mz+mzrange)

  rthigh<-rt[[2]]
  rtlow<-rt[[1]]

  ms1mz <- as.data.frame(MSnbase::precursorMz(rawSpec))
  ms1rt <- as.data.frame(rtime(rawSpec))
  colnames(ms1mz) <- "precursorMz"
  colnames(ms1rt) <- "rtime"

  ms2candid <- subset.data.frame(x = ms1mz,subset = precursorMz>=mzlow & precursorMz<=mzhigh)

  ms2candid$retention <- ms1rt[rownames(ms2candid),]

  ms2matchs <- subset.data.frame(ms2candid, subset = retention>=rtlow&retention<=rthigh)

  ms2matchs$file <- rep(0,nrow(ms2matchs))

  if(nrow(ms2matchs)>0){
    for (j in 1:nrow(ms2matchs)){
      ms2matchs[j,"file"] <- rawSpec@featureData@data[row.names(ms2matchs[j,]),"fileIdx"]
      ms2matchs[j,"file"]<- sampleNames(rawSpec)[as.numeric(ms2matchs[j,"file"])]
    }
  }
  return(ms2matchs)
}

#Example - find MS2 of DNPPE
#DNPPE_ms2 <- LOB_findMS2(rawSpec = rawSpec,mz = 875.550487,rt= c(840,1020),ppm = 2.5)



