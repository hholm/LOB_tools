eval_dupls <- function(flagged_set){
  
  library(tidyverse)
  flagged_set <- flagged %>% mutate(code = "Unknown")
  
  # find xcms_peakgroups that have been multiply assigned
  split_by_peakgroup <- split(flagged_set, 
                       duplicated(flagged_set$xcms_peakgroup) | duplicated(flagged_set$xcms_peakgroup, fromLast = TRUE))
  unique_peakgroups <- split_by_peakgroup[["FALSE"]]
  duplicate_peakgroups <- split_by_peakgroup[["TRUE"]]
  
  # split duplicate assignments by whether they've been identified by retention time factor
  separated_duplicates <- split(duplicate_peakgroups, duplicate_peakgroups$Flag == "ms2v" | 
                       duplicate_peakgroups$Flag == "5%_rtv" | 
                       duplicate_peakgroups$Flag == "10%_rtv")
  
  
  confirmed <- separated_duplicates[["TRUE"]]
  unlikely <- separated_duplicates[["FALSE"]]
  
  double_positives <- split(confirmed, duplicated(confirmed$xcms_peakgroup) | duplicated(confirmed$xcms_peakgroup, fromLast = TRUE))[["TRUE"]]
  double_positives$code <- rep("Confirmed")
  peakgroup_list <- unique(confirmed$xcms_peakgroup)
  
  
  cat("\n")
  for (i in 1:length(unlikely$match_ID)){
    if(grepl(unlikely$xcms_peakgroup[i], peakgroup_list) == TRUE){
      flagged_set$code[flagged_set$match_ID == unlikely$match_ID[i]] <- "False_Assignment"
    }
    cat("\r")
    flush.console()
    cat("Identifying duplicate assignments.","Compound",i,"of",length(unlikely$match_ID),"...")
  }
  cat("Done! (Warnings indicate there were multiple 'TRUE' hits when looking for duplicates.")
  
  for (k in 1:length(confirmed$match_ID)){
    flagged_set$code[flagged_set$match_ID == confirmed$match_ID[k]] <- "Yes"
  }
  
  
  for (m in 1:length(unique_peakgroups$match_ID)){
      if(grepl(c("ms2v", "5%_rtv", "10%_rtv"), unique_peakgroups$Flag[m]) == TRUE | grepl("Yes", unique_peakgroups$lpSolve[m]) == TRUE){
        unique_peakgroups$code[m] <- "Confirmed"
      }else if(grepl("Red", unique_peakgroups$Flag[m]) == TRUE | grepl("No", unique_peakgroups$lpSolve[m]) == TRUE){
        unique_peakgroups$code[m] <- "False_Assignment"
      }else if(grepl("Maybe", unique_peakgroups$lpSolve[m]) == TRUE | grepl("Double_Pek?", unique_peakgroups$Flag[m]) == TRUE){
        unique_peakgroups$code[m] <- "Possible Double Peak"
      }
    }
  
}
