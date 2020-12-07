LOB_plotMS2 <- function(rawSpec, peakdata = NULL, mz = NULL, rt = NULL, rtspan = 175,
                        ppm_pre = 100, ppm = 2.5, samples = NULL, plot_type = "most_scans") {

  # check inputs
  if (!is.null(peakdata)) { # if peakdata isnt NULL
    if (class(peakdata) == "LOBset") { # and it is a LOBset
      peakdata <- LOBSTAHS::peakdata(peakdata) # extract the peakdata.
    } else { # otherwise
      if (class(peakdata) != "data.frame") { # it should be a data.frame
        stop("Input 'peakdata' must be of class 'data.frame'.")
      } else {
        if (!all(c("peakgroup_rt", "LOBdbase_mz", "compound_name") %in% colnames(peakdata))) { # with three columns.
          stop("The input 'peakdata' must have columns 'peakgroup_rt', 'LOBdbase_mz', and 'compound_name'.")
        }
      }
    }
  }


  if (!is.null(peakdata) & (!is.null(mz) | !is.null(rt))) { # Peakdata overides mz and rt slots
    warning("You have provided 'peakdata' as well as 'mz' and/or 'rt' values. 'mz' and 'rt' inputs will be ignored and will be read from 'peakdata'.")
    mz <- NULL
    rt <- NULL
  }

  # Find MS2 spectra scans for lipids
  scans <- LOB_findMS2(
    rawSpec = rawSpec,
    data = peakdata,
    mz = mz,
    rt = rt,
    rtspan = rtspan,
    ppm = ppm_pre
  )

  if (plot_type == "most_scans") {
    # Find the file with the most scans for each lipid
    most <- lapply(
      scans,
      function(x) {
        if (class(x) != "data.frame") {
          "No ms2 spectra found."
        } else {
          names(which(table(x$file) == max(table(x$file))))
        }
      }
    )
  }

  if (plot_type == "highest_int") {
    for (j in 1:length(scans)) {
      if (scans[[j]] == "No ms2 spectra found.") {
        most[j] <- "No ms2 spectra found."
      } else {
        most[j] <- colnames(peakdata[, samples][, unique(scans[[2]]$file)])[which(x == max(x))]
      }
    }
  }

  # plot ms1 chromatogram of lipid
  for (i in 1:length(scans)) {
    cat("\n") # for console feedback
    flush.console()
    cat("Plotting spectra", i, "of", length(scans), "...")

    if (class(scans[[i]]) != "data.frame") { # Dont plot if no scans were found
      cat("\n")
      cat("No ms2 spectra found for mass/lipid", names(scans[i]), "... Moving to next lipid.")
    } else {
      if (!is.null(peakdata)) { # if using a peaklist, set the terms. if not ignore
        mz <- peakdata[i, "LOBdbase_mz"]
        rt <- peakdata[i, "peakgroup_rt"]
      }

      mzrange <- mz * (0.000001 * ppm) # calculate mz range for filtering chrom
      mzlow <- (mz - mzrange)
      mzhigh <- (mz + mzrange)

      plot <- xcms::filterMsLevel(
        filterMz( # filter to only to the one file with the most scans and correct mz range
          xcms::filterFile(rawSpec,
            file = most[[i]][1]
          ),
          mz = c(mzlow, mzhigh)
        ),
        msLevel = 1
      )

      # get the scans IDs from the file with the most scans
      plot_scans <- scans[[i]][which(scans[[i]]$file == most[[i]][1]), ]

      # find the name of the scan in this file that is closest too the center of the rt provided
      closest_scan <- rownames(plot_scans[which(abs(plot_scans$retention - rt) == min(abs(plot_scans$retention - rt))), ])

      # extract a chromatogram from our filtered XCMSnexp object and set non detected ion intensities to 0 for ploting
      df <- xcms::chromatogram(plot)
      df[[1]]@intensity[which(is.na(df[[1]]@intensity))] <- 0


      temp <- tempfile() # create temp file to supress plotting the spectra
      png(filename = temp)
      spec <- plot(rawSpec_ms2[[closest_scan]]) # spectra for the closest scan
      dev.off()
      unlink(temp)

      gridExtra::grid.arrange(
        ggplot() +
          geom_line(aes(
            x = df[[1]]@rtime,
            y = df[[1]]@intensity
          )) +
          xlab("Retention Time") +
          ylab("Intensity") +
          geom_vline(aes(xintercept = plot_scans$retention), color = "blue", alpha = 0.5) +
          geom_vline(aes(xintercept = c(rt + rtspan, rt - rtspan)), color = "green", alpha = 0.75) +
          geom_vline(aes(xintercept = plot_scans[which(abs(plot_scans$retention - rt) == min(abs(plot_scans$retention - rt))), "retention"]), color = "red") +
          ggtitle(as.character(paste("Lipid Name =", names(scans[i]))), subtitle = paste(" M/Z = ", mz, " File = ", most[[i]][1])),
        spec + geom_text(aes(label = round(mtc, 2), y = i), vjust = -0.5)
      )
    }
  }
  return(scans)
}
