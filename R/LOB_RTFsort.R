# LOB_RTFsort

LOB_RTFsort <- function(LOBSet, RT_Factor_Dbase, choose_class = FALSE, plot_data = FALSE, save_plots = FALSE, data_title){

  library(tidyr)
  library(dplyr)
  
  peakdata <- peakdata(LOBSet)
  
  ### Check Inputs ###
  
  if (is.null(RT_Factor_Dbase$Mean_DNPPE_Factor)) {
    
    stop("Looks like input 'RT_Factor_Dbase' is not the right one.\n",
         "Make sure you're using an up to date database with",
         "a Mean DNPPE Factor for each compound.")
    
  }
  
  if (class(LOBSet) != "LOBSet") {
    
    stop("Input LOBSet is not a S4 LOBSet object.")
    
  }
  
  # make sure peakgroup rt is numeric
  original_data$peakgroup_rt <- as.numeric(original_data$peakgroup_rt)
  
  # Extract correct DNPPE retention time
  DNPPE_RT <- as.numeric(original_data$peakgroup_rt[which(grepl("DNPPE", original_data$compound_name))])
  
  if(length(DNPPE_RT) > 1){
    stop("Looks like there's more than one DNPPE peak.",
         "Please resolve in input LOBpeaklist and try again.")
  }
  
}