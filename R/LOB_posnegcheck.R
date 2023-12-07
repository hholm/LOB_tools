# This is a program to cross check LOBSets in positive against negative modes
# to tag peaks identified within in both modes

LOB_posnegcheck <- function(pos_coded_LOBpeaklist,
                             neg_coded_LOBpeaklist,
                             class = NULL,
                             ignore_class = NULL,
                             rt_window = NULL,
                             neg_replacement = NULL){


  cat("Cross checking positive vs. negative ion modes for matching peakgroups.")

  # can't both specify classes and ignore classes
  if(!is.null(class) & !is.null(ignore_class)){
    stop("Get outta here! You can't specify 'class' and 'ignore_class' at the same time, silly!")
  }

  # cut down to a class or list of classes if specified
  # the "else" section of this code prioritizes unoxidized lipids,
  #         except in the case of GSL's, where our predominant known species are mono-hydroxylated
  if(!is.null(class)){
    pos_peaklist <- subset(pos_coded_LOBpeaklist, (pos_coded_LOBpeaklist$species %in% class & pos_coded_LOBpeaklist$degree_oxidation == 1))
    neg_peaklist <- subset(neg_coded_LOBpeaklist, (neg_coded_LOBpeaklist$species %in% class & pos_coded_LOBpeaklist$degree_oxidation == 1))
  }else{
    pos_peaklist <- subset(pos_coded_LOBpeaklist, (pos_coded_LOBpeaklist$degree_oxidation == 0) | (pos_coded_LOBpeaklist$species == "dLCB_GSL_No_FA_OH" & pos_coded_LOBpeaklist$degree_oxidation == 1))
    neg_peaklist <- subset(neg_coded_LOBpeaklist, (neg_coded_LOBpeaklist$degree_oxidation == 0) | (neg_coded_LOBpeaklist$species == "dLCB_GSL_No_FA_OH" & neg_coded_LOBpeaklist$degree_oxidation == 1))
  }

  # remove ignored classes if specified
  # the "else" section of this code prioritizes unoxidized lipids,
  #         except in the case of GSL's, where our predominant known species are mono-hydroxylated
  if(!is.null(ignore_class)){
    pos_peaklist <- subset(pos_coded_LOBpeaklist, (!(pos_coded_LOBpeaklist$species %in% ignore_class) & pos_coded_LOBpeaklist$degree_oxidation == 1))
    neg_peaklist <- subset(neg_coded_LOBpeaklist, (!(neg_coded_LOBpeaklist$species %in% ignore_class) & neg_coded_LOBpeaklist$degree_oxidation == 1))
  }else{
    pos_peaklist <- subset(pos_coded_LOBpeaklist, (pos_coded_LOBpeaklist$degree_oxidation == 0) | (pos_coded_LOBpeaklist$species == "dLCB_GSL_No_FA_OH" & pos_coded_LOBpeaklist$degree_oxidation == 1))
    neg_peaklist <- subset(neg_coded_LOBpeaklist, (neg_coded_LOBpeaklist$degree_oxidation == 0) | (neg_coded_LOBpeaklist$species == "dLCB_GSL_No_FA_OH" & neg_coded_LOBpeaklist$degree_oxidation == 1))
  }

  # if replacing positive lipid classes with their negative counterparts (e.g. DAG's)
  # this section sets up the necessary dataframes
  if(!is.null(neg_replacement)){

    # make dataframes of lipids that won't be crosschecked against the corresponding positive/negative sets
    # first, extra_pos, which is a subset of the full "raw" coded peaklist,
    # including all oxidized lipids, that haven't already been subset into pos_peaklist
    extra_pos_peaklist <- subset(pos_coded_LOBpeaklist, !(pos_coded_LOBpeaklist$match_ID %in% pos_peaklist$match_ID) & !(pos_coded_LOBpeaklist$species %in% neg_replacement))
    #next, do the same thing except it's only the species being included by neg_replacement (e.g. oxidized FFA's and DAG's)
    extra_neg_peaklist <- subset(neg_coded_LOBpeaklist, !(neg_coded_LOBpeaklist$match_ID %in% neg_peaklist$match_ID) & (neg_coded_LOBpeaklist$species %in% neg_replacement))

    # temporarily remove classes where the negative ion mode data will replace positive ion mode data
    # to be cross checked after the positive mode
    pos_to_replace <- subset(pos_peaklist, (!(pos_peaklist$species %in% neg_replacement)))  # positive species to be replaced
    neg_to_replace <- subset(neg_peaklist, ((neg_peaklist$species %in% neg_replacement))) # negative species to replace them

    #set up positive and negative dataframes to cross check
    posneg_crosschecked <- subset(pos_peaklist, (!(pos_peaklist$species %in% neg_replacement)))  # subset pos_peaklist down
    posneg_crosschecked$posneg_check <- "Unknown"

    negpos_crosschecked <- neg_to_replace
    negpos_crosschecked$posneg_check <- "Unknown"

    }else{

      #set up dataframes to cross check
      extra_pos_peaklist <- subset(pos_coded_LOBpeaklist, !(pos_coded_LOBpeaklist$match_ID %in% pos_peaklist$match_ID))
      posneg_crosschecked  <- pos_peaklist
      posneg_crosschecked$posneg_check <- "Unknown"

    }

  # set auto RT window of 20 seconds unless specified
  if(!is.null(rt_window)){
    rt_win <- rt_window
  }else{
    rt_win <- 20
  }


  # cross checking positive set
  for(i in 1:length(posneg_crosschecked$compound_name)){

    # subset peaks in negative mode with the same compound name
    same_compound_name <- neg_peaklist[neg_peaklist$compound_name %in% posneg_crosschecked$compound_name[i],]

    if(length(same_compound_name$match_ID) > 0){

      # calculate acceptable rt range of positive mode rt
      rt_pos <- posneg_crosschecked$peakgroup_rt[i]
      max_rt <- rt_pos + rt_win
      min_rt <- rt_pos - rt_win

      #set up an empty counter to manage negative matches
      negative_matches <- 0

      # check every row with the same compound name
      for(j in length(same_compound_name$peakgroup_rt)){
        if(same_compound_name$peakgroup_rt[j] >= min_rt & same_compound_name$peakgroup_rt[j] <= max_rt){
          # add to the counter
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
  }

  # Now using the same loop to cross check the classes where negative mode classes will replace positive mode
  if(!is.null(neg_replacement)){

    cat("\nNow cross checking negative vs. positive ion modes in selected lipid classes to replace.")

    for(i in 1:length(negpos_crosschecked$compound_name)){


      # i = 2275
      # test <- pos_peaklist[i,]
      #
      #
      # subset peaks in positive mode with the same compound name
      same_compound_name <- pos_to_replace[pos_to_replace$compound_name %in% negpos_crosschecked$compound_name[i],]

      if(length(same_compound_name$match_ID) > 0){
        # calculate acceptable rt range of negative mode rt
        rt_neg <- negpos_crosschecked$peakgroup_rt[i]
        max_rt <- rt_neg + rt_win
        min_rt <- rt_neg - rt_win

        #set up an empty counter to manage positive matches
        positive_matches <- 0

        # check every row with the same compound name
        for(j in length(same_compound_name$peakgroup_rt)){
          if(same_compound_name$peakgroup_rt[j] >= min_rt & same_compound_name$peakgroup_rt[j] <= max_rt){
            # add matching match_ID to the d
            positive_matches <- positive_matches + 1
          }
        }

        if(positive_matches == 0){
          negpos_crosschecked$posneg_check[i] <- "No Match"
        }else if(positive_matches == 1){
          negpos_crosschecked$posneg_check[i] <- "PosNeg Match"
        }else{
          negpos_crosschecked$posneg_check[i] <- "Check Multiple Matches"
        }
      }else{
        negpos_crosschecked$posneg_check[i] <- "No Match"
      }
    }
    #recombine
    extra_neg_peaklist$posneg_check <- "Unknown"
    recombined_neg <- dplyr::bind_rows(negpos_crosschecked, extra_neg_peaklist)
  }


  extra_pos_peaklist$posneg_check <- "Unknown"
  recombined_pos <- dplyr::bind_rows(posneg_crosschecked, extra_pos_peaklist)



  if(!is.null(neg_replacement)){
    recombined <- dplyr::bind_rows(recombined_pos, recombined_neg)
    return(recombined)
  }else{
    return(recombined_pos)
  }

  cat("\nDone!")

}

#######
# library(tidyverse)
#
test <- LOB_posnegcheck(pos_coded_LOBpeaklist = raw_pos_LOBset,
                          neg_coded_LOBpeaklist = raw_neg_LOBset,
                          neg_replacement = c("DAG", "FFA"))
#test <- raw_pos_LOBset #74282

test1 <- test %>% filter(species == "FFA") #269
test2 <- raw_neg_LOBset %>% filter(species == "FFA") #269

test3 <- raw_pos_LOBset %>% filter(species == "DAG") #1410
test4 <- raw_neg_LOBset %>% filter(species == "DAG") #646
test5 <- test %>% filter(species == "DAG")#646

test6 <- raw_pos_LOBset %>% filter(species != "TAG")# 19835
test7 <- raw_pos_LOBset %>% filter(species == "TAG")#2075
test8 <- raw_pos_LOBset #2075
test9 <- CoePro_Final #5996

#########################

