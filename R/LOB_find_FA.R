LOB_find_FA <- function(rawSpec, data = NULL, mz, rt, rtspan = 175, ppm_pre = 100, ppm = 5) {

  # Find MS2 spectra scans for lipids
  scans <- LOBtools::LOB_findMS2(
    rawSpec = rawSpec,
    data = data,
    mz = mz,
    rt = rt,
    rtspan = rtspan,
    ppm = ppm_pre
  )

  # Find the file with the most scans for each lipid
  most <- lapply(
    scans,
    function(x) { if(class(x)!="data.frame"){
      "No ms2 spectra found."
    }else{
      names(which(table(x$file) == max(table(x$file))))
    }
      }
  )

  # plot ms1 chromatogram of lipid
  for (i in 1:length(scans)) {
    cat("\n")
    flush.console()
    cat("Plotting spectra", i, "of", length(scans), "...")

    if(class(scans[[i]])!="data.frame"){
      cat("\n")
      cat("No ms2 spectra found for mass/lipid",names(scans[i]),"... Moving to next lipid.")
    }else{

    if (!is.null(data)) {
      mz <- data[i, "LOBdbase_mz"]
      rt <- data[i, "peakgroup_rt"]
    }

    mzrange <- mz * (0.000001 * ppm)
    mzlow <- (mz - mzrange)
    mzhigh <- (mz + mzrange)

    plot <- rawSpec %>%
      filterFile(file = most[[i]][1]) %>%
      # filterRt(rt = c(rt-rtspan,rt+rtspan)) %>%
      filterMz(mz = c(mzlow, mzhigh), msLevel = 1)

    plot_scans <- scans[[i]][which(scans[[i]]$file == most[[i]][1]), ]
    closest_scan <- rownames(plot_scans[which(abs(plot_scans$retention - rt) == min(abs(plot_scans$retention - rt))), ])

    df <- chromatogram(plot)

    nulldev()
    spec <- plot(rawSpec_ms2[[closest_scan]])
    dev.off()

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
        spec + geom_text(aes(label = round(mtc,2),y = i),vjust = -0.5)
    )
    }
  }
}
