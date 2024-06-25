# functions for detection of retention time patterns in LC-MS data

### Define Helper Functions Common to both Methods ###

screen_lpsolveAPI <- function(X, done, presolve, max_solutions) {
  # for (k in 1:length(unique(X$species))) {

  run <- X
  return <- done
  # Binary string for each point
  Binary_String <- rep(1, nrow(run))

  # Empty String we will build a restrictions from
  Empty_String <- rep(0, nrow(run))

  # Matrix of our Exclusions
  Exclusion_Matrix <- matrix(nrow = 1, ncol = nrow(run))


  # Run a loop to find what to exclude for each point
  cat("\n")
  for (i in 1:nrow(run)) {

    # Get our row
    subject <- run[i, ]
    Exclusion_Table <- run

    # Make a table to store our exclusion info
    Exclusion_Table$Exclude <- rep(FALSE, nrow(run))

    # lets add more info about the peaks to our exclusion table based on info in it
    Exclusion_Table$C_num_diff <- Exclusion_Table$FA_total_no_C - subject$FA_total_no_C
    Exclusion_Table$DB_num_diff <- Exclusion_Table$FA_total_no_DB - subject$FA_total_no_DB
    Exclusion_Table$diff_factor <- Exclusion_Table$C_num_diff / Exclusion_Table$DB_num_diff

    # Lets sort the compounds run above and below our point in terms of rt
    lower_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt < subject$peakgroup_rt), ]
    higher_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt > subject$peakgroup_rt), ]

    # Now find ones that break the rules for lower and higher and set Exclude to TRUE in the Exclusion_Table

    # Exclude the lower with more carbon and less double bonds
    lower_match_ID <- lower_rt[lower_rt$FA_total_no_C >= subject$FA_total_no_C & lower_rt$FA_total_no_DB <= subject$FA_total_no_DB, "match_ID"]
    Exclusion_Table[which(Exclusion_Table$match_ID %in% lower_match_ID), "Exclude"] <- TRUE

    lower_diff <- lower_rt[which(lower_rt$diff_factor > 0 & lower_rt$diff_factor <= 1 & lower_rt$C_num_diff < 0), ]
    lower_diff_exclude <- lower_diff[which(lower_diff$C_num_diff == -1 & lower_diff$DB_num_diff == -1), ]
    if (nrow(lower_diff_exclude) > 0) {
      lower_diff <- lower_diff[!(lower_diff$match_ID %in% lower_diff_exclude$match_ID), ]
    }
    lower_diff_names <- row.names(lower_diff)
    Exclusion_Table[lower_diff_names, "Exclude"] <- TRUE


    # Exclude the higher
    higher_match_ID <- higher_rt[higher_rt$FA_total_no_C <= subject$FA_total_no_C & higher_rt$FA_total_no_DB >= subject$FA_total_no_DB, "match_ID"]
    Exclusion_Table[higher_match_ID, "Exclude"] <- TRUE

    higher_diff <- higher_rt[which(higher_rt$diff_factor > 0 & higher_rt$diff_factor <= 1 & higher_rt$C_num_diff > 0), ]
    higher_diff_exclude <- higher_diff[which(higher_diff$C_num_diff == 1 & higher_diff$DB_num_diff == 1), ]
    if (nrow(higher_diff_exclude) > 0) {
      higher_diff <- higher_diff[!(higher_diff$match_ID %in% higher_diff_exclude$match_ID), ]
    }
    higher_diff_names <- row.names(higher_diff)
    Exclusion_Table[higher_diff_names, "Exclude"] <- TRUE

    # Exclude the compounds with the same name
    Exclusion_Table[which(Exclusion_Table$compound_name == subject$compound_name), "Exclude"] <- TRUE

    Exclusion_String <- Empty_String

    for (j in 1:nrow(run)) {
      if (j != i) {
        if (Exclusion_Table[j, "Exclude"] == TRUE) {
          Exclusion_String <- Empty_String
          Exclusion_String[j] <- 1
          Exclusion_String[i] <- 1
          Exclusion_Matrix <- rbind(Exclusion_Matrix, Exclusion_String)
          rownames(Exclusion_Matrix) <- NULL
          Exclusion_Matrix <- unique(Exclusion_Matrix)
        }
      }
      cat("\r")
      flush.console()
      cat("Writing rules for", as.character(unique(X$species)), "compound number", i, "of", nrow(run), ". Number of Rules created:", nrow(Exclusion_Matrix) - 1, "...")
    }
  }
  cat("Done!")
  Final_Exclusion_Matrix <- Exclusion_Matrix[-1, ]

  if (NROW(Final_Exclusion_Matrix) == 0) {
    cat("\n No rules created. Compound class is to small or any lacks noise to screen. Returning all points as Yes")
    return[return$match_ID %in% run$match_ID, "lpSolve"] <- "Yes"
  } else {
    if (NCOL(Final_Exclusion_Matrix) == 1) {
      rows <- 1
      Final_Exclusion_Matrix <- t(matrix(Final_Exclusion_Matrix))
    } else {
      rows <- nrow(Final_Exclusion_Matrix)
    }


    # time to screen

    cat("\nApplying lpSolveAPI algorythm...")

    num_con <- rows
    num_points <- nrow(run)

    lpmodel <- make.lp(num_con, num_points)
    set.constr.type(lpmodel, rep("<=", num_con))
    set.rhs(lpmodel, rep(1, num_con))
    set.type(lpmodel, columns = c(1:num_points), "binary")
    lp.control(lpmodel, sense = "max", presolve = presolve)
    for (i in 1:num_points) {
      set.column(lpmodel, i, rep(1, length(which(Final_Exclusion_Matrix[, i] == 1))), which(Final_Exclusion_Matrix[, i] == 1))
    }
    set.objfn(lpmodel, run$weights)


    # find first solution
    status <- solve(lpmodel)
    cat(" First solution reached!")
    sols <- matrix(ncol = num_points)
    sols <- sols[-1, ] # create list for more solutions
    obj0 <- get.objective(lpmodel)
    counter <- 1 # construct a counter so you wont get more than 100 solutions
    cat("\n")
    # find more solutions
    while (counter < max_solutions) {
      cat("\r")
      flush.console()
      cat("Looking for more solutions... ",counter+1," of ",max_solutions,"total solutions...")
      sol <- get.variables(lpmodel)
      sols <- rbind(sols, sol)
      add.constraint(lpmodel, 2 * sol - 1, "<=", sum(sol) - 1)
      rc <- solve(lpmodel)
      counter <- counter + 1
      if (status != 0) break
      if (get.objective(lpmodel) < obj0) break
    }

    cat(" Done!")

    bar <- data.frame(colSums(sols))
    run$Type <- rep(0, nrow(run))
    run[which(bar == 0), "Type"] <- "No"
    run[which(bar != nrow(sols) & bar != 0), "Type"] <- "Maybe"
    run[which(bar == nrow(sols)), "Type"] <- "Yes"
    for (g in 1:nrow(run)) {
      return[which(return$match_ID == run[g,"match_ID"]), "lpSolve"] <- run[g,"Type"]
    }
  }
  # }
  return(return)
}

screen_ga <- function(X, done, iters, mutationChance, elitism, popSize, plot_data) {

  # two ga helper functions

  # algorythm evaluation function
  evalFunc <- function(x) {
    current_solution_points <- sum(X[x == 1,"weights"])
    current_lipids <- length(X[x == 1,"weights"])
    if (all(x == 0)) {
      return(0)
    }
    if (any(rowSums(Exclusion_Matrix[, x == 1]) == 2, na.rm = T)) {
      return(0)
    } else {
      return((-1*current_solution_points) + current_lipids)
    }
  }

  # algorythm monitoring function
  monitorFunc <- function(obj,x) {
    cat("\r")
    flush.console()
    cat("Optimizing... Number of Lipids in Best Solution:",sum(obj$population[which(obj$evaluations == min(obj$evaluations))[1],]),
        ", Iterations Remaining:",as.character(sum(is.na(obj$best))))
  }

  run <- X
  return <- done
  # Binary string for each point
  Binary_String <- rep(1, nrow(run))

  # Empty String we will build a restrictions from
  Empty_String <- rep(0, nrow(run))

  # Matrix of our Exclusions
  Exclusion_Matrix <- matrix(nrow = 1, ncol = nrow(run))


  # Run a loop to find what to exclude for each point
  cat("\n")
  for (i in 1:nrow(run)) {

    # Get our row
    subject <- run[i, ]
    Exclusion_Table <- run

    # Make a table to store our exclusion info
    Exclusion_Table$Exclude <- rep(FALSE, nrow(run))

    # lets add more info about the peaks to our exclusion table based on info in it
    Exclusion_Table$C_num_diff <- Exclusion_Table$FA_total_no_C - subject$FA_total_no_C
    Exclusion_Table$DB_num_diff <- Exclusion_Table$FA_total_no_DB - subject$FA_total_no_DB
    Exclusion_Table$diff_factor <- Exclusion_Table$C_num_diff / Exclusion_Table$DB_num_diff

    # Lets sort the compounds run above and below our point in terms of rt
    lower_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt < subject$peakgroup_rt), ]
    higher_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt > subject$peakgroup_rt), ]

    # Now find ones that break the rules for lower and higher and set Exclude to TRUE in the Exclusion_Table

    # Exclude the lower with more carbon and less double bonds
    lower_match_ID <- lower_rt[lower_rt$FA_total_no_C >= subject$FA_total_no_C & lower_rt$FA_total_no_DB <= subject$FA_total_no_DB, "match_ID"]
    Exclusion_Table[which(Exclusion_Table$match_ID %in% lower_match_ID), "Exclude"] <- TRUE

    lower_diff <- lower_rt[which(lower_rt$diff_factor > 0 & lower_rt$diff_factor <= 1 & lower_rt$C_num_diff < 0), ]
    lower_diff_exclude <- lower_diff[which(lower_diff$C_num_diff == -1 & lower_diff$DB_num_diff == -1), ]
    if (nrow(lower_diff_exclude) > 0) {
      lower_diff <- lower_diff[!(lower_diff$match_ID %in% lower_diff_exclude$match_ID), ]
    }
    lower_diff_names <- row.names(lower_diff)
    Exclusion_Table[lower_diff_names, "Exclude"] <- TRUE


    # Exclude the higher
    higher_match_ID <- higher_rt[higher_rt$FA_total_no_C <= subject$FA_total_no_C & higher_rt$FA_total_no_DB >= subject$FA_total_no_DB, "match_ID"]
    Exclusion_Table[higher_match_ID, "Exclude"] <- TRUE

    higher_diff <- higher_rt[which(higher_rt$diff_factor > 0 & higher_rt$diff_factor <= 1 & higher_rt$C_num_diff > 0), ]
    higher_diff_exclude <- higher_diff[which(higher_diff$C_num_diff == 1 & higher_diff$DB_num_diff == 1), ]
    if (nrow(higher_diff_exclude) > 0) {
      higher_diff <- higher_diff[!(higher_diff$match_ID %in% higher_diff_exclude$match_ID), ]
    }
    higher_diff_names <- row.names(higher_diff)
    Exclusion_Table[higher_diff_names, "Exclude"] <- TRUE

    # Exclude the compounds with the same name
    Exclusion_Table[which(Exclusion_Table$compound_name == subject$compound_name), "Exclude"] <- TRUE

    Exclusion_String <- Empty_String

    for (j in 1:nrow(run)) {
      if (j != i) {
        if (Exclusion_Table[j, "Exclude"] == TRUE) {
          Exclusion_String <- Empty_String
          Exclusion_String[j] <- 1
          Exclusion_String[i] <- 1
          Exclusion_Matrix <- rbind(Exclusion_Matrix, Exclusion_String)
          rownames(Exclusion_Matrix) <- NULL
          Exclusion_Matrix <- unique(Exclusion_Matrix)
        }
      }
      cat("\r")
      flush.console()
      cat("Writing rules for", as.character(unique(X$species)), "compound number", i, "of", nrow(run), ". Number of Rules created:", nrow(Exclusion_Matrix) - 1, "...")
    }
  }
  cat("Done!")
  Final_Exclusion_Matrix <- Exclusion_Matrix[-1, ]

  if (NROW(Final_Exclusion_Matrix) == 0) {
    cat("\n No rules created. Compound class is to small or any lacks noise to screen. Returning all points as Yes")
    return[return$match_ID %in% run$match_ID, "lpSolve"] <- "Yes"
  } else {
    if (NCOL(Final_Exclusion_Matrix) == 1) {
      rows <- 1
      Final_Exclusion_Matrix <- t(matrix(Final_Exclusion_Matrix))
    } else {
      rows <- nrow(Final_Exclusion_Matrix)
    }


    # time to screen

    cat("\nApplying Genetic Algorythm...")

    GAmodel <- rbga.bin(size = ncol(Exclusion_Matrix), popSize = popSize, iters = iters,
                        mutationChance = mutationChance, elitism = elitism, evalFunc = evalFunc,
                        monitorFunc = monitorFunc)

    cat(" Done!")

    if(plot_data){plot(GAmodel)}

    filter = GAmodel$evaluations == min(GAmodel$evaluations)
    bestSolution = GAmodel$population[filter, , drop= FALSE][1,]

    run[which(bestSolution == 0), "Type"] <- "No"
    run[which(bestSolution == 1), "Type"] <- "Yes"

    for (g in 1:nrow(run)) {
      return[which(return$match_ID == run[g,"match_ID"]), "lpSolve"] <- run[g,"Type"]
    }
  }
  # }
  return(return)
}

# to show plots
showplots <- function(X, extra) {
  # for (i in 1:length(unique(X$species))){

  # run <- X[X$species == unique(X$species)[i],]
  # run <-run[run$match_ID %in% LOBpeaklist$match_ID,]

  print(ggplot(X, aes(x = peakgroup_rt, y = LOBdbase_mz, color = lpSolve)) +
          scale_color_manual(values = c("Maybe" = "#fdbd4c", "No" = "#e55934", "Yes" = "#9bc53d")) +
          geom_point(size = 3) +
          geom_text(
            label = paste0(
              as.character(X$FA_total_no_C), ":",
              as.character(X$FA_total_no_DB)
            ),
            hjust = 1,
            vjust = 2,
            size = 3,
            color = "black"
          ) +
          ggtitle(paste0("Rt Pattern Screening for ", as.character(X$species), "-", extra)) +
          xlab("Peak Group Retention Time (sec)") +
          ylab("Peak Group m/z"))
}

##this helper function is used to allow the "optim" function to evaluate how well the model fits the points.
evaluate_distance = function(par, df, species){
  ##first we read in the paramaters provided by optim
  rt_30 <- par[1]
  rt_per_acyl <- par[2]
  rt_per_db <- par[3]
  rt_db <- par[4]
  ##initialize distance (the score optim is trying to minimize)
  total_distance <- 0
  ##subset the dataframe to only include the lipid species we are interested in, and no degrees of oxidation
  df <- df[which(df$species == species & df$degree_oxidation == 0),]
  ##initilize the loop
  i <- 1
  n <- 0
  ##this loop evaluates how well the model fits rtf_confirmed points and manually kept points.
  while(i <= nrow(df)){
    if(df$TestSet[i] == TRUE | df$code[i] == "hijacked"){
      n <- n + 1
      mz_distance <- df$LOBdbase_ppm_match[i]
      rt_guess <- rt_30 + rt_per_acyl*(df$FA_total_no_C[i]-30) + rt_per_db*df$FA_total_no_DB[i] + rt_db^df$FA_total_no_DB[i]
      rt_distance <- rt_guess - df$peakgroup_rt[i]
      total_distance <- sqrt(mz_distance^2 + rt_distance^2) + total_distance
    }
    i <- i + 1
  }
  return(total_distance / n)
}

##this helper function will provide the predicted retention times of every lipid in the class,
##based on the model, in a column called rt_prediction
calculate_predictions = function(par, df, species){
  ##this creates the column if there is no column called rt_prediction yet
  if(is.null(df$rt_prediction)){
    df["rt_prediction"] <- NA
  }
  ##first we read in the parameters provided by the model
  rt_30 <- par[1]
  rt_per_acyl <- par[2]
  rt_per_db <- par[3]
  rt_db <- par[4]
  ##now we use the model to calculate predicted rts for every point
  i <- 1
  while(i <= nrow(df)){
    if(df$species[i] == species){
      df$rt_prediction[i] <- rt_30 + rt_per_acyl*(df$FA_total_no_C[i]-30) + rt_per_db*df$FA_total_no_DB[i] + rt_db^df$FA_total_no_DB[i]
    }
    i <- i + 1
  }
  df
}

##This simple helper function computes how well each each point fits the final model.
calculate_scores = function(df){
  df["grid_scores"] <- sqrt((df$rt_prediction - df$peakgroup_rt)^2 + df$LOBdbase_ppm_match^2)
  df
}

##This helper function should be run once for each lipid species and will fit a model and
##calculate how well each point fits that model.
score_grid = function(data, species, params = c(800, 16, -20, 1.3), standard_devs = 2.8,
                      print_out = TRUE, hijacking = FALSE, hijacking_level = 2.8){
  cat("Performing weighting for",species)

  data <- data[which(data$species == species & data$degree_oxidation == 0),]

  # create code column for hijacking
  data["code"] <- FALSE

  #this creates the 'rt_prediction' column if there is none yet.
  if(is.null(data$rt_prediction)){
    data["rt_prediction"] <- NA
  }

  # don't attempt if less than 3 valid points
  if((nrow(data[which(data$TestSet),]) < 3)){
    data$weights <- 1
    data["predict"] <- NA
    cat("Not enough point in known set for weighting. Weighting equally.")
    return(data)
  }

  solution <- optim(par = params, fn = evaluate_distance, df = data, species = species)

  print(solution)

  data <- calculate_predictions(par = solution$par, df = data, species = species)
  data <- calculate_scores(data)

  std <- sd(data[which(data$TestSet == TRUE),"grid_scores"])
  print(std)

  if(hijacking){
    print("hijacking . . .")
    hijacked_data <- data
    hijacked_data[which(data$grid_scores <= hijacking_level*std + solution$value),"code"] <- 'hijacked'
    if(all(hijacked_data$code == data$code)){
      print("no points found for hijacking... hijacking terminated.")
    }else{
      solution <- optim(par = params, fn = evaluate_distance, df = hijacked_data, species = species)
      print(solution)

      data <- calculate_predictions(par = solution$par, df = data, species = species)
      data <- calculate_scores(data)

      std = sd(data[which(data$TestSet == TRUE),"grid_scores"])
      print(std)
    }
  }

  data$predict <- data$grid_scores
  data[which(data$predict <= standard_devs * std + solution$value),"predict"] <- 1
  data[which(data$predict > standard_devs * std + solution$value),"predict"] <- 0

  #data$weights <- data$grid_scores
  data$weights <- (standard_devs * std + solution$value)/data$grid_scores
  #data$weights <- 1/((data$grid_scores + 1))*20
  #data$weights <- (1/((data$grid_scores + 1)) - (1/((standard_devs * std + solution$value)+1)))
  #data$weights <- data$grid_scores

  if(print_out){
    pplot <- ggplot(data, aes(x=peakgroup_rt, y=peakgroup_mz, color = weights)) +
      geom_point()
    print(pplot + scale_color_gradient2(low="green", mid = "orange", high="red"))

    data[which(data$code == "Double Check"),"code"] <- "LP_Solve_Failure"
    pplot <- ggplot(data, aes(x=peakgroup_rt, y=peakgroup_mz, shape=code, color = as.character(code))) +
      geom_point()
    print(pplot)

    pplot <- ggplot(data, aes(x=peakgroup_rt, y=peakgroup_mz, color=predict)) +
      geom_point()
    print(pplot + scale_color_gradient(low="red", high="green"))
  }

  # Need Max to fix these few lines
  # data$grid_scores[which(data$grid_scores > standard_devs * 2 * std + solution$value)] <- mm[2] + 1
  # data$weights <- (mm[2] - data$grid_scores)^3

  print(ggplot(data) + geom_point(aes(weights,grid_scores)) + geom_hline(yintercept = (standard_devs * std + solution$value)))

  return(data)
}

#this was worse
LOB_trim = function(data,species){
  data <- data[which(data$species == species & data$degree_oxidation == 0),]
  d <- data[which(data$lpSolve == "Yes"),]
  for (i in 1:nrow(d)) {
    #find number of neighbors in the db space
    dbn <- (sum(abs(d[which(d$FA_total_no_DB == d[i,"FA_total_no_DB"]),"FA_total_no_C"]-d[i,"FA_total_no_C"]) <= 1)-1)
    d[i,'nieghbours'] <- if(is.na(dbn)){0}else{dbn}
    cn <- (sum(abs(d[which(d$FA_total_no_C == d[i,"FA_total_no_C"]),"FA_total_no_DB"]-d[i,"FA_total_no_DB"]) <= 1)-1)
    d[i,'nieghbours'] <- if(is.na(cn)){d[i,'nieghbours']}else{d[i,'nieghbours']+cn}
    if(d[i,"nieghbours"]<2){
      data[which(data$match_ID == d[i,"match_ID"]),"lpSolve"] <- "Trimmed"
    }
  }
  return(data)
}

# LOB_trim = function(data,species){
#   data <- data[which(data$species == species & data$degree_oxidation == 0),]
#   for (i in 1:length(unique(data$FA_total_no_DB))) {
#     db <- unique(data$FA_total_no_DB)[i]
#     d <- data[which(data$FA_total_no_DB == db & data$lpSolve == "Yes"),]
#     if(nrow(d) > 2){
#       data[data$match_ID == d[which.max(d[,"peakgroup_rt"]),"match_ID"],"lpSolve"] <- "Trimmed"
#       data[data$match_ID == d[which.min(d[,"peakgroup_rt"]),"match_ID"],"lpSolve"] <- "Trimmed"
#     }
#   }
#   return(data)
# }
#

LOB_lpsolveAPI <- function(peakdata, choose_class = NULL, save.files = FALSE, use_ms2 = TRUE, plot_data = FALSE,
                           use_weight = TRUE, hijacking = FALSE, max_solutions = 10, presolve = "mergerows") {

  ### Check Inputs ###

  if (!class(peakdata) %in% c("data.frame","LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }

  #Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    LOBpeaklist <- peakdata
  }else{
    LOBpeaklist <- LOBSTAHS::peakdata(peakdata)
  }

  if (is.null(LOBpeaklist$match_ID)) {
    stop(
      "Peakdata does not contain a 'match_ID' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$compound_name)) {
    stop(
      "Peakdata does not contain a 'compound_name' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$LOBdbase_mz)) {
    stop(
      "Peakdata does not contain a 'LOBdbase_mz' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$peakgroup_rt)) {
    stop(
      "Peakdata does not contain a 'peakgroup_rt' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$FA_total_no_C)) {
    stop(
      "Peakdata does not contain a 'FA_total_no_C' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$FA_total_no_DB)) {
    stop(
      "Peakdata does not contain a 'FA_total_no_DB' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$Flag) & isTRUE(use_ms2)) {
    stop(
      "Peakdata does not contain a 'Flag' column despite 'use_ms2' being set too TRUE. ",
      "Please run LOB_RTFsort on your data.frame or LOBSet to produce Flag column with retention time info."
    )
  }

  if (isFALSE(use_weight)) {
    cat("use_weight is set to FALSE. All lipids will be weighted equally.")
    LOBpeaklist$weights <- rep(1,nrow(LOBpeaklist))
  }

  ### Format our input in a 'run' dataframe

  done <- LOBpeaklist
  done$lpSolve <- rep(NA, length(done$match_ID))
  done$Optim <- rep(NA, length(done$match_ID))
  LOBpeaklist <- LOBpeaklist[which(LOBpeaklist$degree_oxidation == 0), ]

  if (is.null(choose_class) == FALSE) {
    if (any(!(choose_class %in% unique(LOBpeaklist$species)))) {
      stop("Chosen 'choose_class' does not appear in the 'species' column of data.frame.")
    } else {
      LOBpeaklist <- LOBpeaklist[which(LOBpeaklist$species %in% choose_class), ]
    }
  } else {
    LOBpeaklist <- LOBpeaklist[which(is.na(LOBpeaklist$FA_total_no_DB) == FALSE), ]
  }

  for (k in 1:length(unique(LOBpeaklist$species))) {
    LOBrun <- LOBpeaklist[which(LOBpeaklist$species == unique(LOBpeaklist$species)[k]), ]

    LOBrun <- score_grid(LOBrun, species = unique(LOBpeaklist$species)[k], print_out = use_weight, hijacking = hijacking)

    # Put what we need in a dataframe, only include RtF data if it exists
    if (is.null(LOBrun$Flag) == FALSE) {
      PRErun <- data.frame(
        LOBrun$match_ID,
        LOBrun$compound_name,
        LOBrun$LOBdbase_mz,
        LOBrun$peakgroup_rt,
        LOBrun$FA_total_no_C,
        LOBrun$FA_total_no_DB,
        LOBrun$species,
        LOBrun$DNPPE_Factor,
        LOBrun$Flag,
        LOBrun$Flag_RF,
        LOBrun$weights
      )
      # Re-name our column names
      colnames(PRErun) <- c(
        "match_ID",
        "compound_name",
        "LOBdbase_mz",
        "peakgroup_rt",
        "FA_total_no_C",
        "FA_total_no_DB",
        "species",
        "DNPPE_Factor",
        "Flag",
        "Flag_RF",
        "weights"
      )
    } else {
      PRErun <- data.frame(
        LOBrun$match_ID,
        LOBrun$compound_name,
        LOBrun$LOBdbase_mz,
        LOBrun$peakgroup_rt,
        LOBrun$FA_total_no_C,
        LOBrun$FA_total_no_DB,
        LOBrun$species,
        LOBrun$weights
      )
      # Re-name our column names
      colnames(PRErun) <- c(
        "match_ID",
        "compound_name",
        "LOBdbase_mz",
        "peakgroup_rt",
        "FA_total_no_C",
        "FA_total_no_DB",
        "species",
        'weights'
      )
    }

    # If we have elected to use MS2 data, limit our options based on that
    Final_Switch <- FALSE
    if (use_ms2 == TRUE) {
      if (any(PRErun$Flag %in% "5%_rtv" | PRErun$Flag %in% "ms2v")) {
        Final_Switch <- TRUE
        cat("\nPerforming Prescreen of MS2 data for use later for", as.character(unique(LOBpeaklist$species)[k]))

        # find ms2v and %5v points from the data set
        ms2v <- PRErun[which(PRErun$Flag == "ms2v"), ]
        rtv <- PRErun[which(PRErun$Flag == "5%_rtv"), ]
        ms2v <- rbind(ms2v, rtv)

        # Screen for points within that that dont fit and plot the results
        ms2v_screened <- screen_lpsolveAPI(ms2v,done = done,presolve = presolve, max_solutions = max_solutions)
        ms2v_screened_noNA <- ms2v_screened[which(is.na(ms2v_screened$lpSolve) == FALSE & ms2v_screened$species == unique(LOBpeaklist$species)[k]), ]

        if (plot_data == TRUE) {
          showplots(ms2v_screened_noNA, extra = "RtF Confirmed Only")
        }

        # Take only the yes points
        use2screen <- ms2v_screened_noNA[ms2v_screened_noNA$lpSolve == "Yes", ]

        cat("\nWill use", nrow(use2screen), "verified points for screening.")
        cat("\n")

        # For each class eliminate everything that doesn't fit with the yes points

        run <- use2screen
        final_cut <- integer()
        for (i in 1:nrow(run)) {

          # Get our row
          subject <- run[i, ]
          Exclusion_Table <- PRErun

          # Make a table to store our exclusion info
          Exclusion_Table$Exclude <- rep(FALSE, nrow(Exclusion_Table))

          # lets add more info about the peaks to our exclusion table based on info in it
          Exclusion_Table$C_num_diff <- Exclusion_Table$FA_total_no_C - subject$FA_total_no_C
          Exclusion_Table$DB_num_diff <- Exclusion_Table$FA_total_no_DB - subject$FA_total_no_DB
          Exclusion_Table$diff_factor <- Exclusion_Table$C_num_diff / Exclusion_Table$DB_num_diff

          # Lets sort the compounds run above and below our point in terms of rt
          lower_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt < subject$peakgroup_rt), ]
          higher_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt > subject$peakgroup_rt), ]

          # Now find ones that break the rules for lower and higher and set Exclude to TRUE in the Exclusion_Table

          # Exclude the lower
          lower_match_ID <- lower_rt[lower_rt$FA_total_no_C >= subject$FA_total_no_C & lower_rt$FA_total_no_DB <= subject$FA_total_no_DB, "match_ID"]

          lower_diff <- lower_rt[which(lower_rt$diff_factor > 0 & lower_rt$diff_factor <= 1 & lower_rt$C_num_diff < 0), ]
          lower_diff_exclude <- lower_diff[which(lower_diff$C_num_diff == -1 & lower_diff$DB_num_diff == -1), ]
          if (nrow(lower_diff_exclude) > 0) {
            lower_diff <- lower_diff[!(lower_diff$match_ID %in% lower_diff_exclude$match_ID), ]
          }

          final_cut <- c(lower_diff$match_ID, lower_match_ID, final_cut)


          # Exclude the higher
          higher_match_ID <- higher_rt[higher_rt$FA_total_no_C <= subject$FA_total_no_C & higher_rt$FA_total_no_DB >= subject$FA_total_no_DB, "match_ID"]

          higher_diff <- higher_rt[which(higher_rt$diff_factor > 0 & higher_rt$diff_factor <= 1 & higher_rt$C_num_diff > 0), ]
          higher_diff_exclude <- higher_diff[which(higher_diff$C_num_diff == 1 & higher_diff$DB_num_diff == 1), ]
          if (nrow(higher_diff_exclude) > 0) {
            higher_diff <- higher_diff[!(higher_diff$match_ID %in% higher_diff_exclude$match_ID), ]
          }

          final_cut <- c(higher_diff$match_ID, higher_match_ID, final_cut)

          cat("\r")
          flush.console()
          cat("Eliminating points that do not fit. Point number", i, "of", nrow(run), ". Number of points excluded:", length(unique(final_cut)), "...")
        }
        cat("Done!")
        cat("\n")
        ms2_excluded_matchids <- PRErun[PRErun$match_ID %in% final_cut, "match_ID"]
        PRErun <- PRErun[!(PRErun$match_ID %in% final_cut), ]
      } else {
        cat("\nNo points within 5% of a retention time window. Moving to next class")
        next()
      }
    }


    # Screen everything
    cat("\nPerforming screening of complete dataset for", as.character(unique(LOBpeaklist$species)[k]))
    screened <- screen_lpsolveAPI(PRErun,done = done,presolve = presolve, max_solutions = max_solutions)
    if (Final_Switch == TRUE) {
      screened[screened$match_ID %in% ms2_excluded_matchids, "lpSolve"] <- "No"
    }
    screened_plot <- screened[which(is.na(screened$lpSolve) == FALSE & screened$species == unique(LOBpeaklist$species)[k]), ]

    if (plot_data == TRUE) {
      showplots(screened_plot, extra = "Final")
    }

    # if (save.files==TRUE){
    #   ggsave(filename = paste0(as.character(run$species), "_LP_solve.tiff"),
    #          plot = last_plot(),
    #          device = "tiff",
    #          width = 22, height = 17)
    # }
    cat("\n")
    a <- NULL
    if (use_ms2 == TRUE & Final_Switch == TRUE) {
      test <- unique(screened[which(screened$match_ID %in% use2screen$match_ID), "lpSolve"])
      if ("No" %in% test) {
        cat("\n")
        warning("Not all MS2v peaks included in solution.")
        for (a in 1:nrow(screened_plot)) {
          done[which(done$match_ID == screened_plot[a,"match_ID"]), "lpSolve"] <- screened_plot[a,"lpSolve"]
        }
      }
      if ("Maybe" %in% test) {
        cat("\n")
        warning("Not all MS2v peaks included in solution.")
        for (a in 1:nrow(screened_plot)) {
          done[which(done$match_ID == screened_plot[a,"match_ID"]), "lpSolve"] <- screened_plot[a,"lpSolve"]
        }
      }
      cat("\n")
      cat("All MS2v peaks included in solution!")
      for (a in 1:nrow(screened_plot)) {
        done[which(done$match_ID == screened_plot[a,"match_ID"]), "lpSolve"] <- screened_plot[a,"lpSolve"]
      }
    } else {
      for (a in 1:nrow(screened_plot)) {
        done[which(done$match_ID == screened_plot[a,"match_ID"]), "lpSolve"] <- screened_plot[a,"lpSolve"]
      }
      cat("Class Screened without MS2 data.")
    }
    cat("\n")
  }

  if (is.data.frame(peakdata)) {
    return(done)
  }else{
    peakdata@peakdata <- done
    return(peakdata)
  }

}


LOB_rtpattern <- function(peakdata, algorithm = "lpsolve", choose_class = NULL, save.files = FALSE,
                          use_ms2 = TRUE, use_AH = TRUE, plot_data = FALSE, use_weight = TRUE,
                          hijacking = FALSE, iters = 300, mutationChance = 0.01, elitism = T, popSize = 300,
                          max_solutions = 10, presolve = "mergerows", standard_devs = 1) {

  ### Check Inputs ###

  if (!class(peakdata) %in% c("data.frame","LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }

  #Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    LOBpeaklist <- peakdata
  }else{
    LOBpeaklist <- LOBSTAHS::peakdata(peakdata)
  }

  if (is.null(LOBpeaklist$match_ID)) {
    stop(
      "Peakdata does not contain a 'match_ID' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$compound_name)) {
    stop(
      "Peakdata does not contain a 'compound_name' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$LOBdbase_mz)) {
    stop(
      "Peakdata does not contain a 'LOBdbase_mz' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$peakgroup_rt)) {
    stop(
      "Peakdata does not contain a 'peakgroup_rt' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$FA_total_no_C)) {
    stop(
      "Peakdata does not contain a 'FA_total_no_C' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$FA_total_no_DB)) {
    stop(
      "Peakdata does not contain a 'FA_total_no_DB' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$Flag) & isTRUE(use_ms2)) {
    stop(
      "Peakdata does not contain a 'Flag' column despite 'use_ms2' being set too TRUE. ",
      "Please run LOB_RTFsort on your data.frame or LOBSet to produce Flag column with retention time info."
    )
  }

  if (isFALSE(use_weight)) {
    cat("use_weight is set to FALSE. All lipids will be weighted equally.")
    LOBpeaklist$weights <- rep(1,nrow(LOBpeaklist))
    LOBpeaklist$predict <- rep(NA,nrow(LOBpeaklist))
  }

  if (!algorithm %in% c("genalg","lpsolve")) {
    stop(
      "Input 'algorithm' not recognized. Algorithm must be 'genalg' or 'lpsolve'."
    )
  }

  ### Format our input in a 'run' dataframe

  done <- LOBpeaklist
  done$lpSolve <- rep(NA, length(done$match_ID))
  done$Optim <- rep(NA, length(done$match_ID))
  LOBpeaklist <- LOBpeaklist[which(LOBpeaklist$degree_oxidation == 0), ]

  ### Establish the "TestSet" column based on whether we will use AH or MS2 or both.
  # This will be the set of 'good' points all the functions will use or compare with.
  LOBpeaklist$TestSet <- FALSE

  if (use_ms2) {
    LOBpeaklist[which(LOBpeaklist$Flag %in% c('5%_rtv','ms2v')),"TestSet"] <- TRUE
  }

  if (use_AH) {
    LOBpeaklist[which(LOBpeaklist$C2a == 1 | LOBpeaklist$C2b == 1),"TestSet"] <- TRUE
  }

  if (is.null(choose_class) == FALSE) {
    if (any(!(choose_class %in% unique(LOBpeaklist$species)))) {
      stop("Chosen 'choose_class' does not appear in the 'species' column of data.frame.")
    } else {
      LOBpeaklist <- LOBpeaklist[which(LOBpeaklist$species %in% choose_class), ]
    }
  } else {
    LOBpeaklist <- LOBpeaklist[which(is.na(LOBpeaklist$FA_total_no_DB) == FALSE), ]
  }

  for (k in 1:length(unique(LOBpeaklist$species))) {
    LOBrun <- LOBpeaklist[which(LOBpeaklist$species == unique(LOBpeaklist$species)[k]), ]

    if (use_weight) {
      LOBrun <- score_grid(LOBrun, species = unique(LOBpeaklist$species)[k], print_out = plot_data, hijacking = hijacking,standard_devs = standard_devs)
      for (a in 1:nrow(LOBrun)) {
        done[which(done$match_ID == LOBrun[a, "match_ID"]), "MaxOptim_Keep"] <- LOBrun[a,"predict"]
      }
    }

    # Put what we need in a dataframe
      PRErun <- data.frame(
        match_ID = LOBrun$match_ID,
        compound_name = LOBrun$compound_name,
        LOBdbase_mz = LOBrun$LOBdbase_mz,
        peakgroup_rt = LOBrun$peakgroup_rt,
        FA_total_no_C = LOBrun$FA_total_no_C,
        FA_total_no_DB = LOBrun$FA_total_no_DB,
        species = LOBrun$species,
        TestSet = LOBrun$TestSet,
        weights = LOBrun$weights,
        predict = LOBrun$predict
      )

    # If we have elected to use MS2 data, or AH data limit our options based on that
    Final_Switch <- FALSE
    if (use_ms2|use_AH) {
      if (any(PRErun$TestSet)) {
        Final_Switch <- TRUE
        cat("\nPerforming Prescreen of MS2 data for use later for", as.character(unique(LOBpeaklist$species)[k]))

        # find TestSet from the data set
        v <- PRErun[which(PRErun$TestSet), ]

        # Screen for points within that that dont fit and plot the results
        v_screened <- screen_lpsolveAPI(v, done = done, presolve = presolve, max_solutions = max_solutions)
        v_screened_noNA <- v_screened[which(is.na(v_screened$lpSolve) == FALSE & v_screened$species == unique(LOBpeaklist$species)[k]), ]

        if (plot_data == TRUE) {
          showplots(v_screened_noNA, extra = "RtF Confirmed Only")
        }

        # Take only the yes points
        use2screen <- v_screened_noNA[v_screened_noNA$lpSolve == "Yes", ]

        cat("\nWill use", nrow(use2screen), "verified points for screening.")
        cat("\n")

        # For each class eliminate everything that doesn't fit with the yes points

        run <- use2screen
        final_cut <- integer()
        for (i in 1:nrow(run)) {

          # Get our row
          subject <- run[i, ]
          Exclusion_Table <- PRErun

          # Make a table to store our exclusion info
          Exclusion_Table$Exclude <- rep(FALSE, nrow(Exclusion_Table))

          # lets add more info about the peaks to our exclusion table based on info in it
          Exclusion_Table$C_num_diff <- Exclusion_Table$FA_total_no_C - subject$FA_total_no_C
          Exclusion_Table$DB_num_diff <- Exclusion_Table$FA_total_no_DB - subject$FA_total_no_DB
          Exclusion_Table$diff_factor <- Exclusion_Table$C_num_diff / Exclusion_Table$DB_num_diff

          # Lets sort the compounds run above and below our point in terms of rt
          lower_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt < subject$peakgroup_rt), ]
          higher_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt > subject$peakgroup_rt), ]

          # Now find ones that break the rules for lower and higher and set Exclude to TRUE in the Exclusion_Table

          # Exclude the lower
          lower_match_ID <- lower_rt[lower_rt$FA_total_no_C >= subject$FA_total_no_C & lower_rt$FA_total_no_DB <= subject$FA_total_no_DB, "match_ID"]

          lower_diff <- lower_rt[which(lower_rt$diff_factor > 0 & lower_rt$diff_factor <= 1 & lower_rt$C_num_diff < 0), ]
          lower_diff_exclude <- lower_diff[which(lower_diff$C_num_diff == -1 & lower_diff$DB_num_diff == -1), ]
          if (nrow(lower_diff_exclude) > 0) {
            lower_diff <- lower_diff[!(lower_diff$match_ID %in% lower_diff_exclude$match_ID), ]
          }

          final_cut <- c(lower_diff$match_ID, lower_match_ID, final_cut)


          # Exclude the higher
          higher_match_ID <- higher_rt[higher_rt$FA_total_no_C <= subject$FA_total_no_C & higher_rt$FA_total_no_DB >= subject$FA_total_no_DB, "match_ID"]

          higher_diff <- higher_rt[which(higher_rt$diff_factor > 0 & higher_rt$diff_factor <= 1 & higher_rt$C_num_diff > 0), ]
          higher_diff_exclude <- higher_diff[which(higher_diff$C_num_diff == 1 & higher_diff$DB_num_diff == 1), ]
          if (nrow(higher_diff_exclude) > 0) {
            higher_diff <- higher_diff[!(higher_diff$match_ID %in% higher_diff_exclude$match_ID), ]
          }

          final_cut <- c(higher_diff$match_ID, higher_match_ID, final_cut)

          cat("\r")
          flush.console()
          cat("Eliminating points that do not fit. Point number", i, "of", nrow(run), ". Number of points excluded:", length(unique(final_cut)), "...")
        }
        cat("Done!")
        cat("\n")
        ms2_excluded_matchids <- PRErun[PRErun$match_ID %in% final_cut, "match_ID"]
        PRErun <- PRErun[!(PRErun$match_ID %in% final_cut), ]
      } else {
        cat("\nNo points within test set. Moving to next class")
        next()
      }
    }

    # Screen everything
    cat("\nPerforming screening of complete dataset for", as.character(unique(LOBpeaklist$species)[k]))
    if (algorithm == "genalg") {
      screened <- screen_ga(PRErun, done = done, iters = iters, mutationChance = mutationChance,
                            elitism = elitism, popSize = popSize, plot_data = plot_data)
    }

    if (algorithm == "lpsolve") {
    screened <- screen_lpsolveAPI(PRErun,done = done,presolve = presolve, max_solutions = max_solutions)
    }

    if (Final_Switch == TRUE) {
      screened[screened$match_ID %in% ms2_excluded_matchids, "lpSolve"] <- "No"
    }

    screened_plot <- screened[which(is.na(screened$lpSolve) == FALSE & screened$species == unique(LOBpeaklist$species)[k]), ]

    if (plot_data == TRUE) {
      showplots(screened_plot, extra = "Final")
    }

    cat("\n")
    a <- NULL
    if (use_ms2 == TRUE & Final_Switch == TRUE) {
      test <- unique(screened[which(screened$match_ID %in% use2screen$match_ID), "lpSolve"])
      if (any(c("No", "Maybe") %in% test)) {
        cat("\n")
        warning("Not all MS2v peaks included in solution.")
      } else {
        cat("\n")
        cat("All MS2v peaks included in solution!")
      }
    } else {
      cat("Class Screened without MS2 data.")
    }
    for (a in 1:nrow(screened_plot)) {
      done[which(done$match_ID == screened_plot[a, "match_ID"]), "lpSolve"] <- screened_plot[a, "lpSolve"]
    }
    cat("\n")
  }

  if (is.data.frame(peakdata)) {
    return(done)
  }else{
    peakdata@peakdata <- done
    return(peakdata)
  }

}


LOB_rtpattern2 <- function(peakdata, algorithm = "lpsolve", choose_class = NULL, save.files = FALSE,
                          use_ms2 = TRUE, use_AH = TRUE, plot_data = FALSE,
                          hijacking = FALSE, iters = 300, mutationChance = 0.01, elitism = T, popSize = 300,
                          max_solutions = 10, presolve = "mergerows", standard_devs = 1) {

  ### Check Inputs ###

  if (!class(peakdata) %in% c("data.frame","LOBSet")) {
    stop(
      "Input 'peakdata' is not an 'data.frame' or 'LOBset' object.\n",
      "Please use one of these formats for your peakdata."
    )
  }

  #Rename our peak list so we can modify it and keep the complete one
  if (is.data.frame(peakdata)) {
    LOBpeaklist <- peakdata
  }else{
    LOBpeaklist <- LOBSTAHS::peakdata(peakdata)
  }

  if (is.null(LOBpeaklist$match_ID)) {
    stop(
      "Peakdata does not contain a 'match_ID' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$compound_name)) {
    stop(
      "Peakdata does not contain a 'compound_name' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$LOBdbase_mz)) {
    stop(
      "Peakdata does not contain a 'LOBdbase_mz' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$peakgroup_rt)) {
    stop(
      "Peakdata does not contain a 'peakgroup_rt' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$FA_total_no_C)) {
    stop(
      "Peakdata does not contain a 'FA_total_no_C' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$FA_total_no_DB)) {
    stop(
      "Peakdata does not contain a 'FA_total_no_DB' column.",
      "Please use a data.frame generated by 'getLOBpeaklist' or a LOBSet."
    )
  }

  if (is.null(LOBpeaklist$Flag) & isTRUE(use_ms2)) {
    stop(
      "Peakdata does not contain a 'Flag' column despite 'use_ms2' being set too TRUE. ",
      "Please run LOB_RTFsort on your data.frame or LOBSet to produce Flag column with retention time info."
    )
  }

  if (!algorithm %in% c("genalg","lpsolve")) {
    stop(
      "Input 'algorithm' not recognized. Algorithm must be 'genalg' or 'lpsolve'."
    )
  }

  ### Format our input in a 'run' dataframe

  done <- LOBpeaklist
  done$lpSolve <- rep(NA, length(done$match_ID))
  done$Optim <- rep(NA, length(done$match_ID))
  LOBpeaklist <- LOBpeaklist[which(LOBpeaklist$degree_oxidation == 0), ]
  LOBpeaklist$weights <- rep(1,nrow(LOBpeaklist)) # this is a hold over from how I changed the order of functions, see if you can clean up

  ### Establish the "TestSet" column based on whether we will use AH or MS2 or both.
  # This will be the set of 'good' points all the functions will use or compare with.
  LOBpeaklist$TestSet <- FALSE

  if (use_ms2) {
    LOBpeaklist[which(LOBpeaklist$Flag %in% c('5%_rtv','ms2v')),"TestSet"] <- TRUE
  }

  if (use_AH) {
    LOBpeaklist[which(LOBpeaklist$C2a == 1),"TestSet"] <- TRUE
  }

  if (is.null(choose_class) == FALSE) {
    if (any(!(choose_class %in% unique(LOBpeaklist$species)))) {
      stop("Chosen 'choose_class' does not appear in the 'species' column of data.frame.")
    } else {
      LOBpeaklist <- LOBpeaklist[which(LOBpeaklist$species %in% choose_class), ]
    }
  } else {
    LOBpeaklist <- LOBpeaklist[which(is.na(LOBpeaklist$FA_total_no_DB) == FALSE), ]
  }

  for (k in 1:length(unique(LOBpeaklist$species))) {
    LOBrun <- LOBpeaklist[which(LOBpeaklist$species == unique(LOBpeaklist$species)[k]), ]

    # Put what we need in a dataframe
    PRErun <- data.frame(
      match_ID = LOBrun$match_ID,
      compound_name = LOBrun$compound_name,
      LOBdbase_mz = LOBrun$LOBdbase_mz,
      peakgroup_rt = LOBrun$peakgroup_rt,
      FA_total_no_C = LOBrun$FA_total_no_C,
      FA_total_no_DB = LOBrun$FA_total_no_DB,
      species = LOBrun$species,
      TestSet = LOBrun$TestSet,
      weights = LOBrun$weights
    )

    # If we have elected to use MS2 data, or AH data limit our options based on that
    Final_Switch <- FALSE
    if (use_ms2|use_AH) {
      if (any(PRErun$TestSet)) {
        Final_Switch <- TRUE
        cat("\nPerforming Prescreen of MS2 data for use later for", as.character(unique(LOBpeaklist$species)[k]))

        # find TestSet from the data set
        v <- PRErun[which(PRErun$TestSet), ]

        # Screen for points within that that dont fit and plot the results
        v_screened <- screen_lpsolveAPI(v, done = done, presolve = presolve, max_solutions = max_solutions)
        v_screened_noNA <- v_screened[which(is.na(v_screened$lpSolve) == FALSE & v_screened$species == unique(LOBpeaklist$species)[k]), ]

        if (plot_data == TRUE) {
          showplots(v_screened_noNA, extra = "RtF Confirmed Only")
        }

        # Take only the yes points
        use2screen <- v_screened_noNA[v_screened_noNA$lpSolve == "Yes", ]

        cat("\nWill use", nrow(use2screen), "verified points for screening.")
        cat("\n")

        # For each class eliminate everything that doesn't fit with the yes points

        run <- use2screen
        final_cut <- integer()
        for (i in 1:nrow(run)) {

          # Get our row
          subject <- run[i, ]
          Exclusion_Table <- PRErun

          # Make a table to store our exclusion info
          Exclusion_Table$Exclude <- rep(FALSE, nrow(Exclusion_Table))

          # lets add more info about the peaks to our exclusion table based on info in it
          Exclusion_Table$C_num_diff <- Exclusion_Table$FA_total_no_C - subject$FA_total_no_C
          Exclusion_Table$DB_num_diff <- Exclusion_Table$FA_total_no_DB - subject$FA_total_no_DB
          Exclusion_Table$diff_factor <- Exclusion_Table$C_num_diff / Exclusion_Table$DB_num_diff

          # Lets sort the compounds run above and below our point in terms of rt
          lower_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt < subject$peakgroup_rt), ]
          higher_rt <- Exclusion_Table[which(Exclusion_Table$peakgroup_rt > subject$peakgroup_rt), ]

          # Now find ones that break the rules for lower and higher and set Exclude to TRUE in the Exclusion_Table

          # Exclude the lower
          lower_match_ID <- lower_rt[lower_rt$FA_total_no_C >= subject$FA_total_no_C & lower_rt$FA_total_no_DB <= subject$FA_total_no_DB, "match_ID"]

          lower_diff <- lower_rt[which(lower_rt$diff_factor > 0 & lower_rt$diff_factor <= 1 & lower_rt$C_num_diff < 0), ]
          lower_diff_exclude <- lower_diff[which(lower_diff$C_num_diff == -1 & lower_diff$DB_num_diff == -1), ]
          if (nrow(lower_diff_exclude) > 0) {
            lower_diff <- lower_diff[!(lower_diff$match_ID %in% lower_diff_exclude$match_ID), ]
          }

          final_cut <- c(lower_diff$match_ID, lower_match_ID, final_cut)


          # Exclude the higher
          higher_match_ID <- higher_rt[higher_rt$FA_total_no_C <= subject$FA_total_no_C & higher_rt$FA_total_no_DB >= subject$FA_total_no_DB, "match_ID"]

          higher_diff <- higher_rt[which(higher_rt$diff_factor > 0 & higher_rt$diff_factor <= 1 & higher_rt$C_num_diff > 0), ]
          higher_diff_exclude <- higher_diff[which(higher_diff$C_num_diff == 1 & higher_diff$DB_num_diff == 1), ]
          if (nrow(higher_diff_exclude) > 0) {
            higher_diff <- higher_diff[!(higher_diff$match_ID %in% higher_diff_exclude$match_ID), ]
          }

          final_cut <- c(higher_diff$match_ID, higher_match_ID, final_cut)

          cat("\r")
          flush.console()
          cat("Eliminating points that do not fit. Point number", i, "of", nrow(run), ". Number of points excluded:", length(unique(final_cut)), "...")
        }
        cat("Done!")
        cat("\n")
        ms2_excluded_matchids <- PRErun[PRErun$match_ID %in% final_cut, "match_ID"]
        PRErun <- PRErun[!(PRErun$match_ID %in% final_cut), ]
      } else {
        cat("\nNo points within test set. Moving to next class")
        next()
      }
    }

    # Screen everything
    cat("\nPerforming screening of complete dataset for", as.character(unique(LOBpeaklist$species)[k]))
    if (algorithm == "genalg") {
      screened <- screen_ga(PRErun, done = done, iters = iters, mutationChance = mutationChance,
                            elitism = elitism, popSize = popSize, plot_data = plot_data)
    }

    if (algorithm == "lpsolve") {
      screened <- screen_lpsolveAPI(PRErun,done = done,presolve = presolve, max_solutions = max_solutions)
    }

    if (Final_Switch == TRUE) {
      screened[screened$match_ID %in% ms2_excluded_matchids, "lpSolve"] <- "No"
    }

    #screened <- LOB_trim(screened,species = unique(LOBpeaklist$species)[k])
    screened$TestSet <- FALSE
    screened[which(screened$lpSolve == "Yes"),"TestSet"] <- TRUE
    screened <- score_grid(screened, species = unique(LOBpeaklist$species)[k], print_out = TRUE, hijacking = hijacking,standard_devs = standard_devs)

    for (a in 1:nrow(screened)) {
      done[which(done$match_ID == LOBrun[a, "match_ID"]), "MaxOptim_Keep"] <- screened[a,"predict"]
      done[which(done$match_ID == LOBrun[a, "match_ID"]), "MaxOptim_weights"] <- screened[a,"weights"]
    }

    screened_plot <- screened[which(is.na(screened$lpSolve) == FALSE & screened$species == unique(LOBpeaklist$species)[k]), ]

    if (plot_data == TRUE) {
      showplots(screened_plot, extra = "Final")
    }

    cat("\n")
    a <- NULL
    if (use_ms2 == TRUE & Final_Switch == TRUE) {
      test <- unique(screened[which(screened$match_ID %in% use2screen$match_ID), "lpSolve"])
      if (any(c("No", "Maybe") %in% test)) {
        cat("\n")
        warning("Not all MS2v peaks included in solution.")
      } else {
        cat("\n")
        cat("All MS2v peaks included in solution!")
      }
    } else {
      cat("Class Screened without MS2 data.")
    }
    for (a in 1:nrow(screened_plot)) {
      done[which(done$match_ID == screened_plot[a, "match_ID"]), "lpSolve"] <- screened_plot[a, "lpSolve"]
    }
    cat("\n")
  }

  if (is.data.frame(peakdata)) {
    return(done)
  }else{
    peakdata@peakdata <- done
    return(peakdata)
  }

}
