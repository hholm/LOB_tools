LOB_adduct_summary<- function(peakdata){
  
  # check inputs class
  if (!class(peakdata) %in% c("data.frame","LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }
  
  
  #Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    if(c("C1","C2a","C2b") %in% colnames(peakdata)){
      stop("Peakdata does not have adduct diagnostic columns C1, C2a, and C2b needed to summerize adduct higherarchy.")
    }
    flagged_set <- peakdata
  }else{
    flagged_set <- LOBSTAHS::peakdata(peakdata)
  }
  
  #Add adduct summary
  flagged_set[which(flagged_set$C2a == 1),"adduct_summary"] <- "Completely_satisfied"
  flagged_set[which(flagged_set$C2b == 1),"adduct_summary"] <- "Somewhat_satisfied"
  flagged_set[which(flagged_set$C1 == 1),"adduct_summary"] <- "No_other_adducts"
  
  if (is.data.frame(peakdata)) {
    return(flagged_set)
  }else{
    peakdata@peakdata <- flagged_set
    return(peakdata)
  }
}