LOB_plotPosNeg <- function(XCMSnExp_pos, XCMSnExp_neg, peakdata_pos = NULL, adduct_offset = NULL, mz = NULL, rt = NULL, rtspan = 175,
                           ppm = 2.5, file = NULL, window = 1) {

  # check window size
  if (window < 1) {
    stop("Window can not be less than 1 (Full rt window searched for scans).")
  }

  if(is.null(file)){
    stop("Please supply a vector with two filenames too plot from XCMSnExp_pos and XCMSnExp_neg respectively.")
  }

  # check for 'file' in both objects
  if (any(!file[1] %in% MSnbase::sampleNames(XCMSnExp_pos) | !file[2] %in% MSnbase::sampleNames(XCMSnExp_neg))) {
    stop("File(s) '", paste(file[which(!file %in% MSnbase::sampleNames(XCMSnExp))], collapse = ", "), "' not found in both XCMSnExp.
           Check MSnbase::sampleNames(XCMSnExp_pos) and MSnbase::sampleNames(XCMSnExp_neg) to see files in both objects.")
  }

  # check format of peakdata
  if (!is.null(peakdata_pos)) { # if peakdata isnt NULL
    if (class(peakdata_pos) == "LOBset") { # and it is a LOBset
      peakdata_pos <- LOBSTAHS::peakdata(peakdata_pos) # extract the peakdata.
    } else { # otherwise
      if (class(peakdata_pos) != "data.frame") { # it should be a data.frame
        stop("Input 'peakdata' must be of class 'data.frame'.")
      } else {
        if (!all(c("peakgroup_rt", "LOBdbase_mz", "compound_name") %in% colnames(peakdata_pos))) { # with three columns.
          stop("The input 'peakdata' must have columns 'peakgroup_rt', 'LOBdbase_mz', and 'compound_name'.")
        }
      }
    }
  }

  if (!is.null(peakdata_pos) & (!is.null(mz) | !is.null(rt))) { # Peakdata overides mz and rt slots
    warning("You have provided 'peakdata' as well as 'mz' and/or 'rt' values. 'mz' and 'rt' inputs will be ignored and will be read from 'peakdata'.")
    mz <- NULL
    rt <- NULL
  }

  if (is.null(peakdata_pos)) { # if user just supplied mz and rt
    peakdata_pos <- data.frame(LOBdbase_mz = mz, peakgroup_rt = rt, compound_name = as.character(mz))
  }


  range_calc <- function(x) { # a function tp calculate mz range for filtering chromatogram
    range <- x * (0.000001 * ppm)
    low <- (x - range)
    high <- (x + range)
    c(low, high)
  }


  # plot ms1 chromatogram of lipid data
  for (i in 1:nrow(peakdata_pos)) {
    cat("\n") # for console feedback
    flush.console()
    cat("Plotting spectra", i, "of", nrow(peakdata_pos), "...")

    # set rt and mz terms
    mz <- peakdata_pos[i, "LOBdbase_mz"]
    rt <- peakdata_pos[i, "peakgroup_rt"]

    # calculate range of both
    pos_range <- range_calc(mz)
    neg_range <- range_calc(mz + adduct_offset)

    plot_pos <- xcms::filterMsLevel( # filter to only to the one file
      xcms::filterMz( # and correct mz range at MS1
        xcms::filterFile(XCMSnExp_pos,
          file = file[1]
        ),
        mz = pos_range
      ),
      msLevel = 1
    )

    plot_neg <- xcms::filterMsLevel( # repeat for negative
      xcms::filterMz(
        xcms::filterFile(XCMSnExp_neg,
          file = file[2]
        ),
        mz = neg_range
      ),
      msLevel = 1
    )

    # extract a chromatogram from our filtered XCMSnexp objects
    df_pos <- xcms::chromatogram(plot_pos)
    df_neg <- xcms::chromatogram(plot_neg)

    #  set non detected ion intensities to 0 for plotting
    df_pos[[1]]@intensity[which(is.na(df_pos[[1]]@intensity))] <- 0
    df_neg[[1]]@intensity[which(is.na(df_neg[[1]]@intensity))] <- 0

    plot(rbind( # plot both graphs
      ggplot2::ggplotGrob(ggplot() +
          geom_line(aes(
            x = df_pos[[1]]@rtime,
            y = df_pos[[1]]@intensity
          )) +
          xlab("Retention Time") +
          ylab("Intensity") +
          xlim(rt - rtspan * window, rt + rtspan * window) +
          geom_vline(aes(xintercept = c(rt + rtspan, rt - rtspan)), color = "green", alpha = 0.75) +
          ggtitle(as.character(paste("Lipid Name =", peakdata_pos[i, "compound_name"]," Mode = Positive")),
                  subtitle = paste(" M/Z = ", mz, " File = ", file[1]," PPM =",ppm))),
      ggplot2::ggplotGrob(ggplot() +
          geom_line(aes(
            x = df_neg[[1]]@rtime,
            y = df_neg[[1]]@intensity
          )) +
          xlab("Retention Time") +
          ylab("Intensity") +
          xlim(rt - rtspan * window, rt + rtspan * window) +
          geom_vline(aes(xintercept = c(rt + rtspan, rt - rtspan)), color = "green", alpha = 0.75) +
          ggtitle(as.character(paste("Lipid Name =", peakdata_pos[i, "compound_name"]," Mode = Negative")),
                  subtitle = paste(" M/Z = ", mz, " File = ", file[2]," PPM =",ppm))
      )))
  }
}
