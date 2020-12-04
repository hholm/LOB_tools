# Use to find possible MS2 from a RawSpec File. Returns a table of scans with file names.
# Rt is a range in seconds

LOB_findMS2 <- function(rawSpec, data = NULL, mz, rt, rtspan = 175, ppm = 100) {

  # Test if the data is entered in a dataframe (data argument) or with mz, rt, ect arguments.
  # Reformate info into temp 'run' df and list of compounds to search for
  if (is.null(data)) {
    run <- data.frame(mz, rt)
    compound_name <- as.character(mz)
  } else {
    run <- data.frame(data$LOBdbase_mz, data$peakgroup_rt)
    compound_name <- as.character(data[, "compound_name"])
  }

  # Create a empty object to store ms2 scans that match our search
  ms2store <- list()

  # Run a search for each compound
  for (i in 1:nrow(run)) {
    cat("\r")
    flush.console()
    cat("Searching for compounds spectra", i, "of", nrow(run), "...")

    # calculate RT, mz, and ppm into high and low ends of the search window
    mz <- run[i, 1]
    rt <- run[i, 2]

    mzrange <- mz * (0.000001 * ppm)
    mzlow <- (mz - mzrange)
    mzhigh <- (mz + mzrange)

    rthigh <- rt + rtspan
    rtlow <- rt - rtspan

    # Extract the MS2 spectra info from the rawspec object and store it in a dataframe
    ms1mz <- as.data.frame(MSnbase::precursorMz(rawSpec))
    ms1rt <- as.data.frame(xcms::rtime(rawSpec))
    colnames(ms1mz) <- "precursorMz"
    colnames(ms1rt) <- "rtime"

    # Subset data frame to things that fall in the mz window
    ms2candid <- subset.data.frame(x = ms1mz, subset = precursorMz >= mzlow & precursorMz <= mzhigh)

    # Attach rt info
    ms2candid$retention <- ms1rt[rownames(ms2candid), ]

    # Subset again using RT window
    ms2matchs <- subset.data.frame(ms2candid, subset = retention >= rtlow & retention <= rthigh)

    # Create empty column for file names from which the ms2 came from
    ms2matchs$file <- rep(0, nrow(ms2matchs))

    # Add file names as long as the df has info in it
    if (nrow(ms2matchs) > 0) {
      for (j in 1:nrow(ms2matchs)) {
        ms2matchs[j, "file"] <- rawSpec@featureData@data[row.names(ms2matchs[j, ]), "fileIdx"]
        ms2matchs[j, "file"] <- MSnbase::sampleNames(rawSpec)[as.numeric(ms2matchs[j, "file"])]
      }
    }

    # If there were no matches return that, otherwise save df of matchs. Return list of df after each compound has been screened.
    if (nrow(ms2matchs) == 0) {
      ms2store[i] <- "No ms2 spectra found."
    } else {
      ms2store[[i]] <- ms2matchs
    }
    names(ms2store)[i] <- compound_name[i]
  }
  return(ms2store)
}



