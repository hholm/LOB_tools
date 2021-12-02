LOB_addRTF <- function(peakdata, LOBdbase, standards, entry, ms2v, only_add_ms2v = FALSE) {

  ### Check Inputs ###

  if (!class(peakdata) %in% c("data.frame", "LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }

  # add check for unique compound names

  cat("Formatting...\n")
  # Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    LOBpeaklist <- peakdata
  } else {
    LOBpeaklist <- LOBSTAHS::peakdata(peakdata)
  }

  # define addition function
  # add_rtf(x){
  #   if (x %in% LOBpeaklist$compound_name) {
  #     LOBpeaklist[which(LOBpeaklist$compound_name == x),"peakgroup_rt"]
  #   }
  # }

  # loop through every standard to create proper entries
  for (i in 1:length(standards)) {

    # check if standard is already in LOBdbase object; create it if it isnt
    if (!(standards[i] %in% names(LOBdbase@rtf))) {
      LOBdbase@rtf[[standards[i]]] <- list()
    }

    # warn if you are overwriting
    if ((entry %in% names(LOBdbase@rtf[[standards[i]]]))) {
      warning(
        "An entry titled '", entry,
        "' already exists for standard ", standard,
        ". New RTFs will overwrite old ones for this entry."
      )
    }

    # add new entry
    LOBdbase@rtf[[standards[i]]][[entry]] <-
      list(
        data = data.frame(
          factor = rep(NA, LOBdbase@num_entries),
          rt = rep(NA, LOBdbase@num_entries),
          ms2_file = rep(NA, LOBdbase@num_entries)
        ),
        std_RT = LOBpeaklist[which(LOBpeaklist$compound_name == standards[i]), "peakgroup_rt"]
      )
  }

  # vector of relivant compounds to run through
  run <- LOBpeaklist$compound_name[which(LOBpeaklist$compound_name %in% LOBdbase@parent_compound_name)]

  # calculate RTF
  for (i in 1:length(run)) {
    # feedback
    cat("\r")
    flush.console()
    cat("Calculating RTFs for compound ", i, " of ", length(run), "...")

    # which statments for subseting
    cmp_pl <- which(LOBpeaklist$compound_name == run[i])
    cmp_lb <- which(LOBdbase@parent_compound_name == run[i])

    for (j in 1:length(standards)) {
      LOBdbase@rtf[[standards[j]]][[entry]]$data[cmp_lb, "rt"] <- LOBpeaklist[cmp_pl, "peakgroup_rt"]
      LOBdbase@rtf[[standards[j]]][[entry]]$data[cmp_lb,"factor"] <- LOBdbase@rtf[[standards[j]]][[entry]]$data[cmp_lb,"rt"] / LOBdbase@rtf[[standards[j]]][[entry]]$std_RT
      LOBdbase@rtf[[standards[j]]][[entry]]$data[cmp_lb,"ms2v"] <- if(LOBpeaklist[cmp_pl,"match_ID"] %in% ms2v){TRUE}else{FALSE}
    }
  }
  cat("Done!\n")
  return(LOBdbase)
}
