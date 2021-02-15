LOB_adduct_summary<- function(peakdata,choose_class = NULL, plot_data = FALSE,include_ox = FALSE){

  # check inputs class
  if (!class(peakdata) %in% c("data.frame","LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }


  #Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    if(!all(c("C1","C2a","C2b") %in% colnames(peakdata))){
      stop("Peakdata does not have adduct diagnostic columns C1, C2a, and C2b needed to summerize adduct higherarchy.")
    }
    og <- peakdata
    flagged_set <- peakdata
  }else{
    flagged_set <- LOBSTAHS::peakdata(peakdata)
    og <- flagged_set
  }

  if(!is.null(choose_class)){
    flagged_set <- subset(flagged_set, species == choose_class)
  }

  #Add adduct summary
  flagged_set[which(flagged_set$C2a == 1),"adduct_summary"] <- "Completely_satisfied"
  flagged_set[which(flagged_set$C2b == 1),"adduct_summary"] <- "Somewhat_satisfied"
  flagged_set[which(flagged_set$C1 == 1),"adduct_summary"] <- "No_other_adducts"

  #
  if(plot_data == TRUE) {
    cat("\n")
    names <- unique(flagged_set$species)
    if (include_ox == FALSE) {
      plot <- subset(flagged_set,degree_oxidation == 0)
    }else{
      plot <- flagged_set
    }
    for (z in 1:length(names)) {
    print(ggplot(data = subset(plot,species == names[z]))+
            geom_point(aes(x = peakgroup_rt, y = LOBdbase_mz, color =  adduct_summary), size = 3) +
            scale_color_manual(values=c("Somewhat_satisfied"="#fdbd4c","Completely_satisfied"="#9bc53d","No_other_adducts"="#5bc0eb")) +
            geom_text(aes(x = peakgroup_rt, y = LOBdbase_mz, label = paste0(FA_total_no_C, ":", FA_total_no_DB), hjust = 1, vjust = 2, color = adduct_summary)) +
            ggtitle(paste0("Adduct Summary for ", names[z])) +
            xlab("Peak Group Retention Time (sec)") +
            ylab("Peak Group m/z"))
      cat("\r")
      flush.console()
      cat("Generating plot", z, "of", length(names), "...")
    }
  }

  for (a in 1:length(unique(og$match_ID))) {
    og[which(og$match_ID == flagged_set[a,"match_ID"]), "adduct_summary"] <- flagged_set[a,"adduct_summary"]
  }

  if (is.data.frame(peakdata)) {
    return(og)
  }else{
    peakdata@peakdata <- og
    return(peakdata)
  }
}
