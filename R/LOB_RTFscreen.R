#### RT Factor Filtering Function
# Functional Draft
# Feb. 20, 2019
# input a raw lobset

# Load data and RT Factor database
# setwd("C:/Users/TSQ/Desktop/Daniel Lowenstein/Nano_Take_Two/Intermediate_Third/")
# original_data <- read.csv("Nano_Intermediate_Third_Pos_Raw_LOBSTAHS_screened_peakdata_2019-06-18T2-27-43_PM-0400.csv")
# RT_Factor_Dbase <-read.csv("C:/Users/TSQ/Desktop/Daniel Lowenstein/Older_Projects/RT_Factors/Hummel RtF Master Database - rtf_data.csv")

LOB_RTFscreen <- function(peakdata, LOBdbase, standard = "DNPPE",entry = NULL, choose_class = NULL,
                        plot_data = FALSE, save_plots = FALSE, data_title) {

  # library(tidyr)
  # library(dplyr)

  ### Check Inputs ###

  # check that lipid input is a LOBSet or a data.frame
  if (!class(peakdata) %in% c("data.frame", "LOBSet")) {
    stop("Input original_data is neither a 'data.frame' nor 'LOBSet'.")
  }

  # read peaklist data from LOBset if necessary
  if (is.data.frame(peakdata)) {
    original_data <- peakdata
  } else {
    original_data <- LOBSTAHS::peakdata(peakdata)
  }

  if (is.null(original_data$match_ID)) {
    stop(
      "Input data.frame does not contain a 'match_ID' column.",
      "Please use a data.frame generated by 'getLOBpeaklist'."
    )
  }

  if (is.null(original_data$compound_name)) {
    stop(
      "Input data.frame does not contain a 'compound_name' column.",
      "Please use a data.frame generated by 'getLOBpeaklist'."
    )
  }

  if (is.null(original_data$LOBdbase_mz)) {
    stop(
      "Input data.frame does not contain a 'LOBdbase_mz' column.",
      "Please use a data.frame generated by 'getLOBpeaklist'."
    )
  }

  if (is.null(original_data$peakgroup_rt)) {
    stop(
      "Input data.frame does not contain a 'peakgroup_rt' column.",
      "Please use a data.frame generated by 'getLOBpeaklist'."
    )
  }

  if (is.null(original_data$FA_total_no_C)) {
    stop(
      "Input data.frame does not contain a 'FA_total_no_C' column.",
      "Please use a data.frame generated by 'getLOBpeaklist'."
    )
  }

  if (is.null(original_data$FA_total_no_DB)) {
    stop(
      "Input data.frame does not contain a 'FA_total_no_DB' column.",
      "Please use a data.frame generated by 'getLOBpeaklist'."
    )
  }

  if (is.null(entry)) {
    warning(
      "No RTF entry choosen. Using average RTF for the standard. Availible entries include: ",names(LOBdbase@rtf[[standard]])
    )
  }

  #store the rtf data for usage later
  RT_Factor_Dbase <- LOBdbase@rtf[[standard]][[entry]]$data

  # make sure peakgroup rt is numeric
  original_data$peakgroup_rt <- as.numeric(original_data$peakgroup_rt)

  if (length("standard") > 1) {
    stop(
      "Input 'standard' should only be one value. Screen using annotations using specific\n",
      "standards by running the function multipule times, subsetting with 'choose_class'."
    )
  }

  # Extract correct standard retention time
  if (class(standard) == "character") {
    if (standard %in% original_data$compound_name) {
      std_rt <- as.numeric(original_data$peakgroup_rt[which(grepl(standard, original_data$compound_name))])
    } else {
      stop(
        "Input 'standard' not found in peakdata compound names.\n",
        "Check the compound_name slot of peakdata for your standard.\n",
        "Alternatively, use a numeric vector in input 'standard' to select a standard by 'match_ID'"
      )
    }
  }

  if (class(standard) == "numeric") {
    if (standard %in% original_data$match_ID) {
      std_rt <- as.numeric(original_data$peakgroup_rt[which(grepl(standard, original_data$match_ID))])
    } else {
      stop(
        "Input 'standard' not found in peakdata match IDs.\n",
        "Check the match_ID slot of peakdata for your standard.\n",
        "Alternatively, use a charactor vector in input 'standard' to select a standard by 'compound_name'"
      )
    }
  }


  if (length(std_rt) > 1) {
    stop(
      "Looks like there's more than one annotation in 'peakdata' for your standard.",
      "Please enter the match_id of the annotation you would like to use in the 'standard' input."
    )
  }

  # Add column for DNPPE factor and an empty one for flagging IF none exists
  if (is.null(original_data$Flag)){
  original_data <- original_data %>%
    dplyr::mutate(DNPPE_Factor = peakgroup_rt / std_rt, Flag = "None", DBase_std_RF = "NA") %>%
    subset(species != "NA")
}
  # isolate major intact polar lipid classes, unoxidized (need to add pigments, etc.)
  Main_Lipids <- original_data %>%
    subset(degree_oxidation == "0") %>%
    subset(lipid_class == "IP_DAG" | lipid_class == "IP_MAG" | lipid_class == "TAG" | lipid_class == "FFA")

  # if you're choosing a class
  if (!is.null(choose_class)) {
    Main_Lipids <- original_data[original_data$species %in% choose_class, ] %>%
      subset(degree_oxidation == "0")
  }

  # isolate oxidized lipids into df
  Ox_Lipids <- original_data %>%
    subset(degree_oxidation > 0)

  # flag known vs unknown verse the database
  cat("Checking for databases entries...")
  Main_Lipids$Flag <- unlist(apply(
    Main_Lipids, 1,
    function(x) {
      if (x["compound_name"] %in% LOBdbase@parent_compound_name) {
        #since LOBdbbase has entry for every adduct just get the highest one (should all be same info)
        if(!is.na(RT_Factor_Dbase$factor[which(LOBdbase@parent_compound_name == x["compound_name"] & LOBdbase@adduct_rank == 1)])){
          "Known"
        }else{
          "Unknown"
        }
      }else{
        "Not In Database"
      }
    }
  ))
  cat("Done!")

  # separate into two dfs
  Known_RtFs <- Main_Lipids %>% subset(Flag == "Known")
  Unknown_RtFs <- Main_Lipids %>% subset(Flag == "Unknown")

  # for each compound in the "Known" df, grab its corresponding row number in RT_Factor_Dbase,
  # then flag it as follows
  # first check if it's in a 10%
  cat("\n")
  for (i in 1:length(Known_RtFs$compound_name)) {
    which_row <- as.numeric(which(grepl(paste0("^", Known_RtFs$compound_name[i], "$"), LOBdbase@parent_compound_name)))
    Known_RtFs$DBase_std_RF[i] <- RT_Factor_Dbase$factor[which_row]

    if (Known_RtFs$DNPPE_Factor[i] < (RT_Factor_Dbase$factor[which_row] * 1.1) & Known_RtFs$DNPPE_Factor[i] > (RT_Factor_Dbase$factor[which_row] * 0.9)) {
      if (Known_RtFs$DNPPE_Factor[i] < (RT_Factor_Dbase$factor[which_row] * 1.05) & Known_RtFs$DNPPE_Factor[i] > (RT_Factor_Dbase$factor[which_row] * 0.95)) {
        if (!is.na(RT_Factor_Dbase$ms2v[which_row]) & RT_Factor_Dbase$ms2v[which_row] == TRUE) {
          Known_RtFs$Flag[i] <- "ms2v"
        } else {
          Known_RtFs$Flag[i] <- "5%_rtv"
        }
      } else {
        (Known_RtFs$Flag[i] <- "10%_rtv")
      }
    } else {
      Known_RtFs$Flag[i] <- "Red"
    }

    cat("\r")
    flush.console()
    cat("Comparing database RtFs to data.", "Compound", i, "of", length(Known_RtFs$compound_name), "...")
  }
  cat("Done!")

  # check for possible double peaks
  Green_Flags <- Known_RtFs %>% subset(Flag == "ms2v" | Flag == "5%_rtv" | Flag == "10%_rtv")
  split_by_compound_name <- split(
    Green_Flags,
    duplicated(Green_Flags$compound_name) | duplicated(Green_Flags$compound_name, fromLast = TRUE)
  )
  Multiple_Flags <- split_by_compound_name[["TRUE"]]
  Multiple_Flags$Flag <- rep("Double_Peak?", length(Multiple_Flags$compound_name))


  cat("\n")
  for (j in 1:length(Multiple_Flags$match_ID)) {
    match_row <- as.numeric(which(grepl(paste0("^", Multiple_Flags$match_ID[j], "$"), Known_RtFs$match_ID)))
    Known_RtFs$Flag[match_row] <- as.character(Multiple_Flags$Flag[j])

    cat("\r")
    flush.console()
    cat("Flagging possible multiple peaks.", "Compound", j, "of", length(Multiple_Flags$match_ID), "...")
  }
  cat("Done!")

  # change levels of colors so they plot in the right order
  # Combined$Flag = factor(Combined$Flag, levels = c("Red", "ms2v", "5%_rtv", "10%_rtv", "Unknown"))

  Flagged_Data <- original_data

  #create new columns only if they don't exist to preserve existing data.
  if(is.null(Flagged_Data$Flag)){ Flagged_Data$Flag <- rep("Unknown", length(Flagged_Data$Flag))}
  if(is.null(Flagged_Data$Std_name_RF)){ Flagged_Data$Std_name_RF <- NA}
  if(is.null(Flagged_Data$Std_rt_RF)){   Flagged_Data$Std_rt_RF <- NA}

  #combine data
  cat("\n")
  for (j in 1:length(Known_RtFs$match_ID)) {
    match_row <- as.numeric(which(grepl(paste0("^", Known_RtFs$match_ID[j], "$"), Flagged_Data$match_ID)))
    Flagged_Data$Flag[match_row] <- as.character(Known_RtFs$Flag[j])
    Flagged_Data$DBase_std_RF[match_row] <- Known_RtFs$DBase_std_RF[j]
    Flagged_Data$Std_name_RF[match_row] <- standard
    Flagged_Data$Std_rt_RF[match_row] <- std_rt

    cat("\r")
    flush.console()
    cat("Combining data.", "Compound", j, "of", length(Known_RtFs$match_ID), "...")
  }
  cat("Done!")

  #insure the correct classes on certain output
  Flagged_Data$Flag_RF <- factor(Flagged_Data$Flag, levels = c("Red", "ms2v", "5%_rtv", "10%_rtv", "Double_Peak?", "Unknown"))
  Flagged_Data$DBase_std_RF <- as.numeric(Flagged_Data$DBase_std_RF)


  if (plot_data == TRUE) {
    cat("\n")
    LOB_plotdata(Flagged_Data,"RTF",choose_class = choose_class)
  }

  # output to match inputed format, whether that is LOBSet of data.frame
  cat("\n\n")
  if (is.data.frame(peakdata)) {
    return(Flagged_Data)
  } else {
    peakdata@peakdata <- Flagged_Data
    return(peakdata)
  }
}
