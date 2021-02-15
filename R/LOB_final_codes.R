LOB_final_codes <- function(peakdata){

  # check inputs class
  if (!class(peakdata) %in% c("data.frame","LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }

  #Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    flagged_set <- peakdata
  }else{
    flagged_set <- LOBSTAHS::peakdata(peakdata)
  }

  #library(tidyverse)
  flagged_set <- flagged_set %>% dplyr::mutate(code = "Unknown")

  # find xcms_peakgroups that have been multiply assigned
  split_by_peakgroup <- split(flagged_set,
                       duplicated(flagged_set$xcms_peakgroup) | duplicated(flagged_set$xcms_peakgroup, fromLast = TRUE))
  unique_peakgroups <- split_by_peakgroup[["FALSE"]]
  duplicate_peakgroups <- split_by_peakgroup[["TRUE"]]

  # split duplicate assignments by whether they've been identified by retention time factor
  separated_duplicates <- split(duplicate_peakgroups, duplicate_peakgroups$Flag == "ms2v" |
                       duplicate_peakgroups$Flag == "5%_rtv" )

  # separate
  confirmed <- separated_duplicates[["TRUE"]]
  unlikely <- separated_duplicates[["FALSE"]]

  if(length(confirmed) > 0){#give a confirmed code
    for (k in 1:length(confirmed$match_ID)){
      confirmed_row <- as.numeric(which(grepl(paste0("^", confirmed$match_ID[k], "$"), flagged_set$match_ID)))
      flagged_set$code[confirmed_row] <- "RTF_Confirmed"
    }

    # split the ones that are still duplicated into a potential isomer list
    double_positives <- split(confirmed, duplicated(confirmed$xcms_peakgroup) | duplicated(confirmed$xcms_peakgroup, fromLast = TRUE))[["TRUE"]]
    double_positives$code <- rep("Probable Isomer")

    # get all of the xcms peakgroups for the confident rtf confirmations
    peakgroup_list <- as.data.frame(confirmed$xcms_peakgroup)

    # for each of the "unlikely" peaks, if you find an xcms peakgroup match in the rtf confirmed peakgroup list,
    # you can confidently say if was a "false assignment"
    cat("\n")
    for (i in 1:length(unlikely$match_ID)){
      if(grepl(unlikely$xcms_peakgroup[i], peakgroup_list) == TRUE){
        disconfirmed_row <- as.numeric(which(grepl(paste0("^", unlikely$match_ID[i], "$"), flagged_set$match_ID)))
        flagged_set$code[disconfirmed_row] <- "False_Assignment"
      }else{unique_peakgroups <- rbind(unique_peakgroups, unlikely[i, ])}
      cat("\r")
      flush.console()
      cat("Identifying duplicate assignments.","Compound",i,"of",length(unlikely$match_ID),"...")
    }
    cat("Done! (Warnings indicate there were multiple 'TRUE' hits when looking for duplicates.)")
  }

  # assigning final codes
  for (m in 1:length(unique_peakgroups$match_ID)){
      if(grepl("ms2v|5%_rtv", unique_peakgroups$Flag[m]) == TRUE){
        unique_peakgroups$code[m] <- "RTF_Confirmed"
      }
  }

  for (m in 1:length(unique_peakgroups$match_ID)){
    if(grepl("10%_rtv|Double_Peak?", unique_peakgroups$Flag[m]) == TRUE){
      unique_peakgroups$code[m] <- "Double Check"
    }
  }

  for (m in 1:length(unique_peakgroups$match_ID)){
    if(grepl("Yes", unique_peakgroups$lpSolve[m]) == TRUE & grepl("Unknown", unique_peakgroups$Flag[m]) == TRUE){
      unique_peakgroups$code[m] <- "LP_Solve_Confirmed"
    }
  }

  for (m in 1:length(unique_peakgroups$match_ID)){
    if(grepl("Red", unique_peakgroups$Flag[m]) == TRUE){
      unique_peakgroups$code[m] <- "RTF_Failure"
    }
  }

  for (m in 1:length(unique_peakgroups$match_ID)){
    if(grepl("Unknown", unique_peakgroups$Flag[m]) == TRUE & grepl("No", unique_peakgroups$lpSolve[m]) == TRUE){
      unique_peakgroups$code[m] <- "LP_Solve_Failure"
    }
  }

  for (m in 1:length(unique_peakgroups$match_ID)){
    if(grepl("Unknown", unique_peakgroups$Flag[m]) == TRUE & grepl("Maybe", unique_peakgroups$lpSolve[m]) == TRUE){
      unique_peakgroups$code[m] <- "LP_Solve_Maybe"
    }
  }

  cat("\n")
  for (j in 1:length(unique_peakgroups$match_ID)){
    match_row <- as.numeric(which(grepl(paste0("^", unique_peakgroups$match_ID[j], "$"), flagged_set$match_ID)))
    flagged_set$code[match_row] = as.character(unique_peakgroups$code[j])
    cat("\r")
    flush.console()
    cat("Recombining unique peakgroups.","Compound",j,"of",length(unique_peakgroups$match_ID),"...")
  }
  cat("Done!")
  cat("\n")
  cat("Summarizing adduct info...")
  flagged_set <- LOB_adduct_summary(flagged_set)
  cat("Done!")

  #Adding a column to indicate whether isomers have been resolved
  # flagged_set$resolved <- rep(NA,length(flagged_set$code))
  # i <- NULL
  # for (i in 1:length(unique(flagged_set$xcms_peakgroup))) {
  #   run <- flagged_set[which(flagged_set$xcms_peakgroup == unique(flagged_set$xcms_peakgroup)[i]),]
  #   if(nrow(run)=1){
  #     flagged_set[which(flagged_set$xcms_peakgroup == unique(flagged_set$xcms_peakgroup)[i]),"resolved"] <- TRUE
  #   }else{
  #     if("run$code")
  #   }
  # }

  if (is.data.frame(peakdata)) {
    return(flagged_set)
  }else{
    peakdata@peakdata <- flagged_set
    return(peakdata)
  }

}

