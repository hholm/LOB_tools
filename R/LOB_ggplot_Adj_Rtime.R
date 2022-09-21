LOB_ggplot_Adj_Rtime <- function(object, adjustedRtime = TRUE) {
  
  #load helper function (since this isn't in xcms package anymore)
  hackyRtAdjustment <- function(x, rtraw, rtadj) {
    ## re-order everything if rtraw is not sorted; issue #146
    if (is.unsorted(rtraw)) {
      idx <- order(rtraw)
      rtraw <- rtraw[idx]
      rtadj <- rtadj[idx]
    }
    adjFun <- stepfun(rtraw[-1] - diff(rtraw) / 2, rtadj)
    res <- adjFun(x)
    ## Fix margins.
    idx_low <- which(x < rtraw[1])
    if (length(idx_low)) {
      first_adj <- idx_low[length(idx_low)] + 1
      res[idx_low] <- x[idx_low] + res[first_adj] - x[first_adj]
    }
    idx_high <- which(x > rtraw[length(rtraw)])
    if (length(idx_high)) {
      last_adj <- idx_high[1] - 1
      res[idx_high] <- x[idx_high] + res[last_adj] - x[last_adj]
    }
    if (is.null(dim(res)))
      names(res) <- names(x)
    res
  }
  
  #start process (just copied from plotAdjustedRtime)
  if (!is(object, "XCMSnExp")) {
    stop("'object' has to be an 'XCMSnExp' object.")
  }
  if (!hasAdjustedRtime(object)) {
    warning("No alignment/retention time correction results present.")
  }
  diffRt <- rtime(object, adjusted = TRUE) - rtime(object,
    adjusted = FALSE
  )
  diffRt <- split(diffRt, fromFile(object))
  xRt <- rtime(object, adjusted = adjustedRtime, bySample = TRUE)
  if (length(col) == 1) {
    col <- rep(col, length(diffRt))
  }
  if (length(lty) == 1) {
    lty <- rep(lty, length(diffRt))
  }
  if (length(lwd) == 1) {
    lwd <- rep(lwd, length(diffRt))
  }
  if (length(col) != length(diffRt)) {
    warning(
      "length of 'col' does not match the number of samples! Will ",
      "use 'col[1]' for all samples."
    )
    col <- rep(col[1], length(diffRt))
  }
  if (length(lty) != length(lty)) {
    warning(
      "length of 'lty' does not match the number of samples! Will ",
      "use 'lty[1]' for all samples."
    )
    lty <- rep(lty[1], length(diffRt))
  }
  if (length(lwd) != length(lwd)) {
    warning(
      "length of 'lwd' does not match the number of samples! Will ",
      "use 'lwd[1]' for all samples."
    )
    lwd <- rep(lwd[1], length(diffRt))
  }
  
  #restructure data for lines and save
  forplot <- data.frame(xRt)
  colnames(forplot) <- sampleNames(res)
  forplot <- tidyr::pivot_longer(forplot,cols = tidyr::everything())
  forplot[,"diffRt"] <- tidyr::pivot_longer(data.frame(diffRt),cols = tidyr::everything())[,"value"]

  ph <- processHistory(res,type = "Retention time correction")
  if (length(ph)) {
    ph <- ph[[length(ph)]]
    if (is(ph, "XProcessHistory")) {
      prm <- processParam(ph)
      if (is(prm, "PeakGroupsParam")) {
        rm(diffRt)
        rm(xRt)
        rawRt <- rtime(object, adjusted = FALSE, bySample = TRUE)
        adjRt <- rtime(object, adjusted = TRUE, bySample = TRUE)
        pkGroup <- peakGroupsMatrix(prm)
        subs <- subset(prm)
        if (!length(subs)) {
          subs <- 1:ncol(pkGroup)
        }
        rawRt <- rawRt[subs]
        adjRt <- adjRt[subs]
        pkGroupAdj <- pkGroup
        for (i in 1:ncol(pkGroup)) {
          pkGroupAdj[, i] <- hackyRtAdjustment(pkGroup[
            ,
            i
          ], rawRt[[i]], adjRt[[i]])
        }
        diffRt <- pkGroupAdj - pkGroup
        if (adjustedRtime) {
          xRt <- pkGroupAdj
        } else {
          xRt <- pkGroup
        }

        #restructure data for points
        forpoint <- tidyr::pivot_longer(data.frame(xRt),cols = tidyr::everything())
        forpoint[,"diffRt"] <- tidyr::pivot_longer(data.frame(diffRt),cols = tidyr::everything())[,"value"] 
        
        #plot
        ggplot2::ggplot() +
          ggplot2::geom_path(ggplot2::aes(forplot$value, forplot$diffRt, group = forplot$name, color = forplot$name),size = 0.6) +
          ggplot2::theme_minimal() +
          ggplot2::geom_point(ggplot2::aes(forpoint$value,forpoint$diffRt),shape = 4) +
          ggplot2::labs(x = if(adjustedRtime == TRUE){"Adjusted RT"}else{"Raw RT"}, y = "RT Difference", color = "Samples")
      }
    }
  }
}
