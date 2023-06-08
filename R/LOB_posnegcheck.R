# This is a program to cross check LOBSets in positive against negative modes
# to tag peaks identified within in both modes

LOB_posnegcheck <- function(pos_coded_LOBpeaklist,
                             neg_coded_LOBpeaklist,
                             class = NULL,
                             ignore_class = NULL,
                             rt_window = NULL){


  cat("Cross checking positive vs. negative ion modes for matching peakgroups.")
  cat("\n")

  # can't both specify classes and ignore classes
  if(!is.null(class) & !is.null(ignore_class)){
    stop("Get outta here! You can't specify 'class' and 'ignore_class' at the same time, silly!")
  }

  # cut down to a class or list of classes if specified
  if(!is.null(class)){
    pos_peaklist <- subset(pos_coded_LOBpeaklist, (pos_coded_LOBpeaklist$species %in% class & pos_coded_LOBpeaklist$degree_oxidation == 1))
    neg_peaklist <- subset(neg_coded_LOBpeaklist, (neg_coded_LOBpeaklist$species %in% class & pos_coded_LOBpeaklist$degree_oxidation == 1))
  }else{
    pos_peaklist <- pos_coded_LOBpeaklist
    neg_peaklist <- neg_coded_LOBpeaklist
  }

  # remove ignored classes if specified
  if(!is.null(ignore_class)){
    pos_peaklist <- subset(pos_coded_LOBpeaklist, (!(pos_coded_LOBpeaklist$species %in% ignore_class) & pos_coded_LOBpeaklist$degree_oxidation == 1))
    neg_peaklist <- subset(neg_coded_LOBpeaklist, (!(neg_coded_LOBpeaklist$species %in% ignore_class) & neg_coded_LOBpeaklist$degree_oxidation == 1))
  }else{
    pos_peaklist <- pos_coded_LOBpeaklist
    neg_peaklist <- neg_coded_LOBpeaklist
  }

  #extra_pos_peaklist <- subset(pos_coded_LOBpeaklist, !(pos_coded_LOBpeaklist$match_ID %in% pos_peaklist$match_ID))

  # set auto RT window of 20 seconds unless specified
  if(!is.null(rt_window)){
    rt_win <- rt_window
  }else{
    rt_win <- 15
  }

  posneg_crosschecked <- pos_peaklist
  posneg_crosschecked$posneg_check <- "Unknown"

  for(i in 1:length(posneg_crosschecked$match_ID)){


    # i = 2275
    # test <- pos_peaklist[i,]
    #
    #
    # subset peaks in negative mode with the same compound name
    same_compound_name <- neg_peaklist[neg_peaklist$compound_name %in% pos_peaklist$compound_name[i],]

    if(length(same_compound_name$match_ID) > 0){
      # calculate acceptable rt range of positive mode rt
      rt_pos <- pos_peaklist$peakgroup_rt[i]
      max_rt <- rt_pos + rt_win
      min_rt <- rt_pos - rt_win

      #set up an empty counter to manage negative matches
      negative_matches <- 0

      # check every row with the same compound name
      for(j in length(same_compound_name$peakgroup_rt)){
        if(same_compound_name$peakgroup_rt[j] >= min_rt & same_compound_name$peakgroup_rt[j] <= max_rt){
          # add matching match_ID to the d
          negative_matches <- negative_matches + 1
        }
      }

      if(negative_matches == 0){
        posneg_crosschecked$posneg_check[i] <- "No Match"
      }else if(negative_matches == 1){
        posneg_crosschecked$posneg_check[i] <- "PosNeg Match"
      }else{
        posneg_crosschecked$posneg_check[i] <- "Check Multiple Matches"
      }
    }else{
      posneg_crosschecked$posneg_check[i] <- "No Match"
    }
    cat("\r")
    flush.console()
    cat("Checking assignment", i, "of", nrow(pos_peaklist))
  }

  #extra_pos_peaklist$posneg_check <- "Unknown"
  #recombined <- bind_rows(posneg_crosschecked, extra_pos_peaklist)

  cat("\nDone!")
  #return(recombined)
  return(posneg_crosschecked)

}


