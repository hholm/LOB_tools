LOB_addRTF <- function(peakdata, LOBdbase, match_ID, standards = "DNPPE", dataset, ms2v = FALSE, ms2file = NA) {

  ### Check Inputs ###

  if (!class(peakdata) %in% c("data.frame", "LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }

  # add check for unique compound names
  # Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    LOBpeaklist <- peakdata
  } else {
    LOBpeaklist <- LOBSTAHS::peakdata(peakdata)
  }

  #rt entries
  rt <- data.frame(compound_name = LOBpeaklist[which(!LOBpeaklist$compound_name %in% standards),"compound_name"],
                    value = LOBpeaklist[which(!LOBpeaklist$compound_name %in% standards),"peakgroup_rt"],
                    value_type  = "rt",
                    dataset = dataset,
                    standard = NA)

  #rtf entries
  rtfs <- data.frame()
  for (i in 1:length(standards)) {
    r <- rt
    stdrt <- LOBpeaklist[which(LOBpeaklist$compound_name == standards[i]),"peakgroup_rt"]
    r$value <- r$value/stdrt
    r$value_type  <- "rtf"
    r$standard <- standards[i]
    rtfs <- rbind(r,rtfs)
  }
  LOBdbase@rtf <- rbind(LOBdbase@rtf,rt,rtfs)
  return(LOBdbase)

}
