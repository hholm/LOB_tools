LOB_plotMS2 <- function(XCMSnExp, peakdata = NULL, plot_file = "most_scans", mz = NULL, rt = NULL, rtspan = 175,
                        ppm_pre = 100, ppm = 2.5, file = NULL, window = 1, diagnostic = NULL, diagnostic_ppm = 15, NL = NULL) {

  # check inputs
  if (window < 1) {
    stop("Window can not be less than 1 (Full rt window searched for scans).")
  }

  if (!(plot_file %in% c("most_scans", "highest_int"))) {
    cat("Input 'plot_file' not a charactor vector reading either 'most_scans' or 'highest_int'. Treating as file name.")
    file <- plot_file
  }

  if (is.null(file)) {
    if (any(!file %in% MSnbase::sampleNames(XCMSnExp))) {
      stop("File(s) '", paste(file[which(!file %in% sampleNames(XCMSnExp))], collapse = ", "), "' not found in XCMSnExp. Check sampleNames(XCMSnExp) to see files in object.")
    }
  }

  if (!is.null(peakdata)) { # if peakdata isnt NULL
    if (class(peakdata) == "LOBset") { # and it is a LOBset
      peakdata <- LOBSTAHS::peakdata(peakdata) # extract the peakdata.
    } else { # otherwise
      if (class(peakdata) != "data.frame") { # it should be a data.frame
        stop("Input 'peakdata' must be of class 'data.frame'.")
      } else {
        if (!all(c("peakgroup_rt", "LOBdbase_mz", "compound_name") %in% colnames(peakdata))) { # with at least three columns.
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


  range_calc <- function(x) {
    range <- x * (0.000001 * diagnostic_ppm) # calculate mz range for filtering chrom
    low <- (x - range)
    high <- (x + range)
    c(low, high)
  }

  if (!is.null(diagnostic)) {
    if (is.null(names(diagnostic))) {
      stop("'diagnostic' must be a named numeric.")
    } else {
      if (anyDuplicated((names(diagnostic))) != 0) {
        stop("'diagnostic' names must be unique")
      }
    }
    # calculate diagnostic mz ranges
    drange <- lapply(diagnostic, range_calc)
  }

  if (!is.null(NL)) {
    if (is.null(names(NL))) {
      stop("'NL' must be a named numeric.")
    } else {
      if (anyDuplicated((names(NL))) != 0) {
        stop("'NL' names must be unique.")
      }
    }
    # calculate NL mz ranges
    nlrange <- lapply(NL, range_calc)
  }

  # Find MS2 spectra scans for lipids
  scans <- LOB_findMS2(
    XCMSnExp = XCMSnExp,
    peakdata = peakdata,
    mz = mz,
    rt = rt,
    rtspan = rtspan,
    ppm_pre = ppm_pre
  )

  # if-then statments too select file to plot scans from
  if (!is.null(file)) { # check if user specified a file
    most <- lapply(
      scans,
      function(x) {
        if (class(x) != "data.frame") {
          "No ms2 spectra found."
        } else {
          file
        }
      }
    )
  } else {
    if (plot_file == "most_scans") {
      # Find the file with the most scans for each lipid
      most <- lapply(
        scans,
        function(x) {
          if (class(x) != "data.frame") {
            "No ms2 spectra found."
          } else {
            names(which(table(x$file) == max(table(x$file))))[1]
          }
        }
      )
    }

    if (plot_file == "highest_int") {
      # Find the file with the highest precursor int for each lipid
      most <- lapply(
        scans,
        function(x) {
          if (class(x) != "data.frame") {
            "No ms2 spectra found."
          } else {
            x[which(x$precursorIntensity == max(x$precursorIntensity)), "file"][1]
          }
        }
      )
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
      for (z in 1:length(most[[i]])) {
        if (!(most[[i]][z] %in% scans[[i]]$file)) {
          cat("\n")
          cat("No ms2 spectra found for mass/lipid", names(scans[i]), "in file specified... Moving to next lipid.")
          next
        }

        if (!is.null(peakdata)) { # if using a peaklist, set the terms. if not ignore
          mz <- peakdata[i, "LOBdbase_mz"]
          rt <- peakdata[i, "peakgroup_rt"]
        }

        mzrange <- mz * (0.000001 * ppm) # calculate mz range for filtering chrom
        mzlow <- (mz - mzrange)
        mzhigh <- (mz + mzrange)

        plot <- xcms::filterMsLevel(
          xcms::filterMz( # filter to only to the one file with the most scans and correct mz range
            xcms::filterFile(XCMSnExp,
              file = most[[i]][z]
            ),
            mz = c(mzlow, mzhigh)
          ),
          msLevel = 1
        )

        # get the scans IDs from the file with the most scans
        plot_scans <- scans[[i]][which(scans[[i]]$file == most[[i]][z]), ]

        # find the name of the scan in this file that is closest too the center of the rt provided
        closest_scan <- rownames(plot_scans[which(abs(plot_scans$retentionTime - rt) == min(abs(plot_scans$retentionTime - rt))), ])

        # extract a chromatogram from our filtered XCMSnexp object and set non detected ion intensities to 0 for ploting
        df <- xcms::chromatogram(plot)
        df[[1]]@intensity[which(is.na(df[[1]]@intensity))] <- 0


        temp <- tempfile() # create temp file to suppress plotting the spectra
        png(filename = temp)
        spec <- plot(XCMSnExp[[closest_scan]]) # spectra for the closest scan
        dev.off()
        unlink(temp) # delete file

        if (!is.null(diagnostic)) { # find fragments
          which <- lapply(drange, function(x) {
            which(spec$data$mtc > x[1] & spec$data$mtc < x[2])
          })
          diff <- sapply(unlist(which), function(x) {
            (spec$data$mtc[x] - diagnostic[names(which(which == x))][[1]]) / (diagnostic[names(which(which == x))][[1]] * 0.000001)
          })
        } else {
          diff <- NULL
        }

        if (!is.null(NL)) { # find fragments
          losses <- (mz - spec$data$mtc)
          which_nl <- lapply(nlrange, function(x) {
            which(losses > x[1] & losses < x[2])
          })
          diff_nl <- sapply(unlist(which_nl), function(x) {
            (losses[x] - NL[names(which(which_nl == x))][[1]]) / (NL[names(which(which_nl == x))][[1]] * 0.000001)
          })
          if (length(diff_nl) > 0) {
            names(diff_nl) <- paste("NL", names(diff_nl), sep = "")
          }
        } else {
          diff_nl <- NULL
        }

        splits <- c(seq(min(spec$data$mtc), max(spec$data$mtc), by = 10), max(spec$data$mtc)) # create a vector to split up the MS2 spectra by mz

        bincode <- as.character(cut(spec$data$mtc, breaks = splits, include.lowest = TRUE)) # create these bins

        i_plot <- rep(NA, length(spec$data$i)) # NA vector

        suppressWarnings( # find the highest ms2 peak every ~10 m/z and mark that for plotting with a number
          for (k in 1:length(unique(bincode))) {
            bin <- unique(bincode)[k]
            sub <- spec$data$i[which(bincode == bin)]
            i_plot[which(spec$data$i == sub[which(sub == max(sub))])] <- sub[which(sub == max(sub))]
          }
        )

        i_plot[which(i_plot == 0)] <- NA # if a bin was all 0s make sure they are NA so they dont plot

        colors <- rep("black", length(i_plot)) # make colors for labels
        if (!is.null(diagnostic)) { # make sure diagnostic fragments are plotted and in red
          i_plot[unlist(which)] <- spec$data$i[unlist(which)]
          colors[unlist(which)] <- "red"
        }
        if (!is.null(NL)) { # make sure NL fragments are plotted and in green
          i_plot[unlist(which_nl)] <- spec$data$i[unlist(which_nl)]
          colors[unlist(which_nl)] <- "green"
        }

        diff_group <- c(diff, diff_nl)
        plot(rbind( # plot both graphs
          ggplot2::ggplotGrob(ggplot() +
            geom_line(aes(
              x = df[[1]]@rtime,
              y = df[[1]]@intensity
            )) +
            xlab("Retention Time") +
            ylab("Intensity") +
            xlim(rt - rtspan * window, rt + rtspan * window) +
            geom_vline(aes(xintercept = plot_scans$retentionTime), color = "blue", alpha = 0.5) +
            geom_vline(aes(xintercept = c(rt + rtspan, rt - rtspan)), color = "green", alpha = 0.75) +
            geom_vline(aes(xintercept = plot_scans[closest_scan, "retentionTime"]), color = "red") +
            ggtitle(as.character(paste("Lipid Name =", names(scans[i]))), subtitle = paste(" M/Z = ", mz, " File = ", most[[i]][z]))),
          ggplot2::ggplotGrob(spec +
            geom_text(aes(label = round(mtc, 2), y = i_plot), color = colors, vjust = -0.5) +
            ggtitle(paste(spec$labels$title, "; Collision Energy", XCMSnExp[[closest_scan]]@collisionEnergy, "; Retention Time", XCMSnExp[[closest_scan]]@rt, "sec")) +
            if (!is.null(diagnostic) | !is.null(NL)) {
              ggplot2::annotate(
                geom = "text",
                label = paste(
                  "Fragments\n",
                  if (length(diff_group) > 0) {
                    paste(names(diff_group), "=",
                      round(unlist(diff_group), digits = 5),
                      collapse = "\n"
                    )
                  } else {
                    NULL
                  }
                ), x = Inf, y = Inf, hjust = 1, vjust = 1
              )
            } else {
              NULL
            })
        ))
      }
    }
  }
  cat("\n\n")
  return(scans)
}
