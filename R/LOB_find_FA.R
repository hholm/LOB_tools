LOB_find_FA <- function(rawSpec, data = NULL, mz, rt, rtspan = 175, ppm = 100){
  
# Find MS2 spectra scans for lipids
  scans  <- LOB_findMS2(rawSpec = rawSpec,
                        data = data,
                        mz = mz,
                        rt = rt,
                        rtspan = rtspan,
                        ppm = ppm)
  
# Find the file with the most scans for each lipid
  most <- lapply(done,
                 function(x){names(which(table(x$file) == max(table(x$file))))})

# plot ms1 chromatogram of lipid
  for (i in 1:length(done)) {
  
  if(!is.null(data)){
  mz <- data[i,"LOBdbase_mz"]
  rt <- data[i,"peakgroup_rt"]
  }
    
  mzrange <- mz * (0.000001 * ppm)
  mzlow <- (mz - mzrange)
  mzhigh <- (mz + mzrange)
  
test  <- rawSpec_ms2 %>%
    filterFile(file = most[[i]][1]) %>%
    filterRt(rt = c(rt-rtspan,rt+rtspan)) %>%
    filterMz(mz = c(mzlow,mzhigh),msLevel = 1)

plotpeak <- filterMsLevel(test, msLevel = 1)
plot(featureData(plotpeak)@data$retentionTime,featureData(plotpeak)@data$basePeakIntensity)

plot(featureData(test))
%>%
    plot(type = "XIC")
  
  }
}