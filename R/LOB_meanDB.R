
### Double Bond Average Computation ###

LOB_meanDB <- function(LOBpeaklist, samples, class = NULL){
  
#Cut down to certian classes if desired 
  if(!is.null(class)){
      run <- LOBpeaklist[which(LOBpeaklist$species %in% class),]
  }else{
      run <- LOBpeaklist
  }

done <- data.frame()

for (j in 1:length(unique(run$species))){
 
  # Get rows of certain class, Get sample columns that are samples 
  run_samps <- LOBpeaklist[which(LOBpeaklist$species==unique(run$species)[j]),samples]
  
  # Make a matrix of percents
  run_samps <- t(t(run_samps)/rowSums(t(run_samps)))
  run_samps <- as.data.frame(run_samps)
  
  ### Add the db values back in.
  run_samps$FA_total_no_DB <- LOBpeaklist[which(LOBpeaklist$species==unique(run$species)[j]),"FA_total_no_DB"]
  
  ### Now a for loop to sum our percents with same DBs
  ### First make a storage frame
  storage<- data.frame(matrix(nrow=(ncol(run_samps)), ncol=16))
  
  for (i in 0:15) {
    summed<- as.data.frame(colSums(run_samps[which(run_samps[,"FA_total_no_DB"]==i),]))
    summed["FA_total_no_DB",]<-i
    colnames(summed)<- i
    storage[,i+1] <-summed
    colnames(storage)[i+1] <- i
    row.names(storage) <- row.names(summed)
  }
  
  ###add a column of 0s for the averages
  storage$dbaverage <- rep(0)
  
  ###calculate the weighted average for each sample (applyed by row)
  for(y in 1:nrow(storage)-1){
    av<-(storage[y,"0"]*0)+(storage[y,"1"]*1)+(storage[y,"2"]*2)+(storage[y,"3"]*3)+(storage[y,"4"]*4)+(storage[y,"5"]*5)+(storage[y,"6"]*6)+(storage[y,"7"]*7)+(storage[y,"8"]*8)+(storage[y,"9"]*9)+(storage[y,"10"]*10)+(storage[y,"11"]*11)+(storage[y,"12"]*12)+(storage[y,"13"]*13)+(storage[y,"14"]*14)+(storage[y,"15"]*15)
    
    ###store the average
    storage[["dbaverage"]][y] <- av
  }

  if(j==1){
    done <- data.frame(storage$dbaverage)
    colnames(done) <- paste(unique(run$species)[j],"_db_mean",sep = "")
  }else{
    out <- data.frame(storage$dbaverage)
    colnames(out) <- paste(unique(run$species)[j],"_db_mean",sep = "")
    done <- cbind(done,out)
  }
}

rownames(done) <- rownames(storage)
return(done)

}



