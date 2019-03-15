# Add final column to LP-solved and flagged document

Combine_Filters <- function(solved_flagged_lobset){
  
  
  if (is.null(solved_flagged_lobset$Flag)) {
    
    stop("Looks like you're not using a DNPPE retention time factor flagged dataset.\n",
         "Make sure you're using a set that's been filtered by the RF database.")
  }
  
  if (is.null(solved_flagged_lobset$lpSolve)) {
    
    stop("Looks like your dataset has not been through the lpSolve algorithm.",
         "Make sure you're using a set that has an lpSolve column.")
  }
  
  solved_flagged_lobset$Points <- rep(0, length(solved_flagged_lobset$match_ID))
  
  for (i in 1:length(solved_flagged_lobset$Points)){
   if (solved_flagged_lobset$Flag[i] == "ms2v"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 4
    }
  
    if (solved_flagged_lobset$Flag[i] == "5%_rtv"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 3
    }
  
    if (solved_flagged_lobset$Flag[i] == "10%_rtv"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 2
    }
    
    if (solved_flagged_lobset$Flag[i] == "Unknown"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 1
    }
    
    if (solved_flagged_lobset$lpSolve[i] == "Good_Fit"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 2
    }
    
    if (is.na(solved_flagged_lobset$lpSolve[i]) == TRUE){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 0
    }else if (solved_flagged_lobset$lpSolve[i] == "Good_Fit"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 2
    }else if (solved_flagged_lobset$lpSolve[i] == "Maybe"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 1
    }else if (solved_flagged_lobset$lpSolve[i] == "No_Fit"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 0
    }
    
    
    if (solved_flagged_lobset$Flag[i] == "Red"){
      solved_flagged_lobset$Points[i] <-  solved_flagged_lobset$Points[i] + 0
    }
   
    
    
    
  }
  return(solved_flagged_lobset)
}

# 
# flags
# 
#   4    ms2v
#   3    5% - no ms2
#   2    10% - no ms2
#   1    unknown rtf
#   0    red
# 
#   2    lpSolve- fit
#   1    lpSolve maybe
#   0    lpSolve red
