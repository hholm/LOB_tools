# Merge RT Factor Dataframes
# 2/19/19

# Match rows from the new ok'ed data to the main Hummel database
# by LOBdbase_mz, and add the new values into new columns



#setwd("C:/Users/TSQ/Desktop/Daniel Lowenstein/Older_Projects/RT_Factors/")



########################
# Time to just write a function to do it 12/5/19
########################

# How it will work:
# Load the current Hummel database and your new verified LOBSTAHS output into R
# The inputs to this function will be those two sheets, plus a short identifying string to go in the header
# that lets us know what dataset these rt's came from
# The function will add two new columns: rt_"name" and DNPPE_rtf_"name", 
# and it will edit the "mean_DNPPE_Factor" by pulling from every column that starts with "DNPPE_rtf"

Update_Hummel_database <- function(NewLobset, OldDatabase, NewName){
  
  # create two new column names with the new name string
  new_column_1 <- paste0("rt_", NewName)
  new_column_2 <- paste0("DNPPE_rtf_", NewName)
  
  # create a new database, adding a couple of columns at the end, then rename those columns with the new column names you just made
  NewDatabase <- OldDatabase
  NewDatabase$new_column_a <- rep(NA, length(NewDatabase[1]))
  NewDatabase$new_column_b <- rep(NA, length(NewDatabase[1]))
  names(NewDatabase)[names(NewDatabase) == "new_column_a"] <- new_column_1
  names(NewDatabase)[names(NewDatabase) == "new_column_b"] <- new_column_2
  
  
  # set up a new df to get average retention times of each compound
  retention_times <- data.frame(compound_name = unique(NewLobset$compound_name),
                                times = rep(NA, length(unique(NewLobset$compound_name))))
  retention_times$times <- sapply(retention_times$compound_name, function(x){mean(NewLobset$peakgroup_rt[which(NewLobset$compound_name == x)])})
  
  # get DNPPE rt from your new lobset to calc rtf in a minute
  DNPPE_rt <- as.numeric(retention_times$times[which(retention_times$compound_name == "DNPPE")])
  
  # get a list of columns in the database with rtfs in it to calculate the average
  factor_columns <- which(grepl(x = colnames(NewDatabase), pattern = "^DNPPE_rtf_"))
  
  
  # in order
  # enter the new average retention time of each compound in your lobset into the second to last column in the new database 
  # divide that by your DNPPE rt to calculate your rtf
  for (i in 1:length(retention_times$compound_name)){
    NewDatabase[which(as.character(NewDatabase$compound_name) == as.character(retention_times$compound_name[i])), (length(colnames(NewDatabase))-1)] <- retention_times[i, 2]
    NewDatabase[which(as.character(NewDatabase$compound_name) == as.character(retention_times$compound_name[i])), (length(colnames(NewDatabase)))] <- (retention_times[i, 2]/DNPPE_rt)
    
  }
  
  # calculate average rtfs for each row
  # and add sd, min, and max values to make plotting and quality control easier
  # may still add the same values in seconds instead of rtfs to make plots easier to read.
  for(j in 1:length(NewDatabase$compound_name)){
    if(!all(is.na(NewDatabase[j, factor_columns]))){
      NewDatabase$Mean_DNPPE_Factor[j] <- mean(as.numeric(NewDatabase[j, factor_columns]), na.rm = TRUE)
      NewDatabase$StDev_DNPPE_Factor[j] <- sd(as.numeric(NewDatabase[j, factor_columns]), na.rm = TRUE)
      NewDatabase$Min_DNPPE_Factor[j] <- min(as.numeric(NewDatabase[j, factor_columns]), na.rm = TRUE)
      NewDatabase$Max_DNPPE_Factor[j] <- max(as.numeric(NewDatabase[j, factor_columns]), na.rm = TRUE)
    }else{
      next
    }
  }
  
  return(NewDatabase)
}

####################################
## Testing
####################################
# 
# OldDatabase <- read.csv("C:/Users/TSQ/Desktop/Daniel Lowenstein/Older_Projects/RT_Factors/Hummel RtF Master Database - rtf_data.csv")
# NewLobset <- read.csv("C:/Users/TSQ/Desktop/Daniel Lowenstein/Henry_Antarctica/Henry_Antarctica_Coded.csv")
# new_name <- "NAAMES_3_4"
# 
#   
# #########################
# 
# # create two new column names with the new name string
# new_column_1 <- paste0("rt_", new_name)
# new_column_2 <- paste0("DNPPE_rtf_", new_name)
# 
# # create a new database, adding a couple of columns at the end, then rename those columns with the new column names you just made
# NewDatabase <- OldDatabase
# NewDatabase$new_column_a <- rep(NA, length(NewDatabase[1]))
# NewDatabase$new_column_b <- rep(NA, length(NewDatabase[1]))
# names(NewDatabase)[names(NewDatabase) == "new_column_a"] <- new_column_1
# names(NewDatabase)[names(NewDatabase) == "new_column_b"] <- new_column_2
# 
# 
# # set up a new df to get average retention times of each compound
# retention_times <- data.frame(compound_name = unique(NewLobset$compound_name),
#                               times = rep(NA, length(unique(NewLobset$compound_name))))
# retention_times$times <- sapply(retention_times$compound_name, function(x){mean(NewLobset$peakgroup_rt[which(NewLobset$compound_name == x)])})
# 
# # get DNPPE rt from your new lobset to calc rtf in a minute
# DNPPE_rt <- as.numeric(retention_times$times[which(retention_times$compound_name == "DNPPE")])
# 
# # get a list of columns in the database with rtfs in it to calculate the average
# factor_columns <- which(grepl(x = colnames(NewDatabase), pattern = "^DNPPE_rtf_"))
# 
# 
# # in order
# # enter the new average retention time of each compound in your lobset into the second to last column in the new database 
# # divide that by your DNPPE rt to calculate your rtf
# # calculate a mean DNPPE_rtf[actor] for each compound
# for (i in 1:length(retention_times$compound_name)){
#   NewDatabase[which(as.character(NewDatabase$compound_name) == as.character(retention_times$compound_name[i])), (length(colnames(NewDatabase))-1)] <- retention_times[i, 2]
#   NewDatabase[which(as.character(NewDatabase$compound_name) == as.character(retention_times$compound_name[i])), (length(colnames(NewDatabase)))] <- (retention_times[i, 2]/DNPPE_rt)
# 
# }
# 
# for(j in 1:length(NewDatabase$compound_name)){
#   if(!all(is.na(NewDatabase[j, factor_columns]))){
#     NewDatabase$Mean_DNPPE_Factor[j] <- mean(as.numeric(NewDatabase[j, factor_columns]), na.rm = TRUE)
#     NewDatabase$StDev_DNPPE_Factor[j] <- sd(as.numeric(NewDatabase[j, factor_columns]), na.rm = TRUE)
#     NewDatabase$Min_DNPPE_Factor[j] <- min(as.numeric(NewDatabase[j, factor_columns]), na.rm = TRUE)
#     NewDatabase$Max_DNPPE_Factor[j] <- max(as.numeric(NewDatabase[j, factor_columns]), na.rm = TRUE)
#   }else{
#     next
#   }
# }
# 
# 
# ###################
# # plot testing
# ###################
# 
# lipidclass <- new_test %>%
#   mutate(C_DB = paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"))) %>%
#   filter(is.na(FA_total_no_C) == FALSE)
# 
# 
# lipid_classes <- unique(lipidclass$species)
# 
# 
# for (i in 1:length(lipid_classes)){
#   Lipid <- lipidclass %>%
#     filter(species == paste(lipid_classes[i]))
#   
#   print(ggplot(Lipid, aes(x = Mean_DNPPE_Factor, y = LOBdbase_mz, color =  FA_total_no_DB))+
#           geom_point()+
#           geom_text(aes(label = C_DB, hjust = 1, vjust = 2))+
#           ggtitle(paste0("M/Z vs. RT in Hummel Database Update", lipid_classes[i])))+
#     geom_errorbarh(aes(xmax = Max_DNPPE_Factor, xmin = Min_DNPPE_Factor))
#   
#   ggsave(filename = paste0(lipid_classes[i], "_Hummel_Dbase_Updates.tiff"),
#          plot = last_plot(),
#          device = "tiff",
#          width = 22, height = 17)
#   
#   
#   
# }
