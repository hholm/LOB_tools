# Use to find possible MS2 from a RawSpec File. Returns a table of scans with file names.
# Rt is a range in seconds

LOB_findMS2 <- function(rawSpec, peakdata = NULL, mz, rt, rtspan = 175, ppm = 100) {

  # Test if the data is entered in a dataframe (data argument) or with mz, rt, ect arguments.
  # Reformate info into temp 'run' df and list of compounds to search for
  if (is.null(peakdata)) {
    run <- data.frame(mz, rt)
    compound_name <- as.character(mz)
  } else {
    if (class(peakdata) != "data.frame") {
      stop("Input 'peakdata' must be of class 'data.frame'")
    } else {
      if (!all(c("peakgroup_rt", "LOBdbase_mz", "compound_name") %in% colnames(peakdata))) { # with three columns.
        stop("The input 'peakdata' must have columns 'peakgroup_rt', 'LOBdbase_mz', and 'compound_name'.")
      }
      run <- data.frame(peakdata$LOBdbase_mz, peakdata$peakgroup_rt)
      compound_name <- as.character(peakdata[, "compound_name"])
    }
  }

  cat("Loading data...")
  cat("\n")
  # Extract the MS2 spectra info from the rawspec object and store it in a dataframe
  ms1 <- filterMsLevel(rawSpec, msLevel = 2)@featureData@data[, c("fileIdx", "precursorMZ", "retentionTime")]
  samples <- MSnbase::sampleNames(rawSpec)
  ms1$file <- sapply(ms1$fileIdx, function(x) {
    samples[x]
  })
  ms1 <- ms1[, -1]

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

    # Subset data frame to things that fall in the mz window
    ms2candid <- subset.data.frame(x = ms1, subset = precursorMZ >= mzlow & precursorMZ <= mzhigh)

    # Subset again using RT window
    ms2candid <- subset.data.frame(ms2candid, subset = retentionTime >= rtlow & retentionTime <= rthigh)

    # If there were no matches return that, otherwise save df of matchs. Return list of df after each compound has been screened.
    if (nrow(ms2candid) == 0) {
      ms2store[i] <- "No ms2 spectra found."
    } else {
      ms2store[[i]] <- ms2candid
    }
    names(ms2store)[i] <- compound_name[i]
  }
  return(ms2store)
}
