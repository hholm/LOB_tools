LOB_plotdata <- function(peakdata,plot_type,choose_class = NULL) {

  ### Check Inputs ###

  if (!class(peakdata) %in% c("data.frame", "LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }
  
  if(!plot_type %in% c("RTF")){
    stop(
      "Input 'plot_type' must be a 'character' and contain one or more of the following:\n",
      "RTF"
    )
  }
  
  # Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    LOBpeaklist <- peakdata
  } else {
    LOBpeaklist <- LOBSTAHS::peakdata(peakdata)
  }

  # plot RTF data
  if (plot_type %in% "RTF") {
    
    if(is.null(LOBpeaklist$Flag_RF)){"No RTF screened data to plot."}
    
    # subset to classes based on inputs
    if (is.null(choose_class)) {
      # get lipid classes for which there are rtf DB entries
      lipid_classes <- unique(LOBpeaklist[which(!is.na(LOBpeaklist$Flag_RF)), "species"])
      if(length(lipid_classes)==0){stop("No classes detected to have RTF screening data. Specify a lipid class to force plotting.")}
    } else {
      lipid_classes <- choose_class
    }

    # and plot each
    for (i in 1:length(lipid_classes)) {
      #insert a class check
      print(ggplot(data = subset(LOBpeaklist, species == lipid_classes[i] & degree_oxidation == 0)) +
              geom_point(aes(x = peakgroup_rt, y = LOBdbase_mz, color = Flag), size = 3) +
              scale_color_manual(values = c("Red" = "#e55934", "ms2v" = "#fdbd4c", "5%_rtv" = "#9bc53d", "10%_rtv" = "#0f5f20", "Double_Peak?" = "#5bc0eb", "Unknown" = "#ff99ff")) +
              geom_point(aes(x = DBase_std_RF * Std_rt_RF, y = LOBdbase_mz, color = Flag), shape = 3) +
              geom_errorbarh(aes(xmax = DBase_std_RF * Std_rt_RF * 1.1, xmin = DBase_std_RF * Std_rt_RF * 0.9, height = 0.2, y = LOBdbase_mz, color = Flag)) +
              geom_text(aes(x = peakgroup_rt, y = LOBdbase_mz, label = paste0(FA_total_no_C, ":", FA_total_no_DB), hjust = 1, vjust = 2, color = Flag)) +
              ggtitle(paste0("RtF Screening for ", lipid_classes[i])) +
              xlab("Peak Group Retention Time (sec)") +
              ylab("Peak Group m/z"))

      cat("\r")
      flush.console()
      cat("Generating plot", i, "of", length(lipid_classes), "...")
    }
    cat("Done!")
  }
}
