#Use to find possible MS2 from a RawSpec File. Returns a table of scans with file names.
# Rt is a range in seconds

LOB_findMS2 <- function(rawSpec,data=NULL,mz,rt,rtspan=175,ppm=2.5){

  if (is.null(data)) {
   run <- data.frame(mz,rt)
   compound_name <- as.character(mz)
  }else{
   run <- data.frame(data$LOBdbase_mz,data$peakgroup_rt)
   compound_name <- as.character(data[,"compound_name"])
  }

  ms2store <- list()

  for (i in 1:nrow(run)) {
    cat("\r")
    flush.console()
    cat("Searching for compounds spectra",i,"of",nrow(run),"...")

  mz <- run[i,1]
  rt <- run[i,2]

  mzrange <- mz*(0.000001*ppm)
  mzlow <- (mz-mzrange)
  mzhigh <- (mz+mzrange)

  rthigh<-rt+rtspan
  rtlow<-rt-rtspan

  ms1mz <- as.data.frame(MSnbase::precursorMz(rawSpec))
  ms1rt <- as.data.frame(xcms::rtime(rawSpec))
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
  if (nrow(ms2matchs)==0) {
    ms2store[i] <- "No ms2 spectra found."
  }else{
  ms2store[[i]] <- ms2matchs
  }
  names(ms2store)[i] <- compound_name[i]
  }
  return(ms2store)
}



