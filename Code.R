#######################
#### Preparation ######
#######################

# if you don't have bain installed, install bain , and load library
# install.packages("bain")


##################################
### Analysis individual BFs ######
##################################

library("bain")

#read in data
data <- read.table("data.txt", header = TRUE)
# data <- read.table(file = "https://raw.githubusercontent.com/fayetteklaassen/gpbf/master/data.txt", header = TRUE)
# Determine the number of unique ppnrs = the number of cases
names(data)
# [1] "ppnr"           "TimePerception" "Valence"        "Arousal"       
# [5] "Condition" 
N <- length(unique(data$ppnr)) 

# create an empty list to store results
results <- vector("list", length = N) 
names(data)

set.seed(7561) # seed to create replicable results

for(i in 1:N) {  # loop over N individuals
  data_i <- data[data$ppnr == i,]   # subset data for ppnr == i
  
  fit_i <- lm(formula = TimePerception ~ Condition + Valence + Arousal, 
              data = data_i)   # execute linear model 
  # save the results of Bain analysis. 
  results[[i]] <-  bain(fit_i, "Condition=0 & Valence>0 & Arousal>0; 
                            Condition>0 & Valence>0 & Arousal>0") 
}

###########################
### viewing the output ####
###########################

names(results[[1]]) # view the names of the bain output for first person ([[1]]).
# view the output of fit and BFmatrix
results[[1]]$fit 
results[[1]]$BFmatrix 

output <- matrix(0, nrow = N, ncol = 2) # create output table with N rows and 4 columns
colnames(output) <- c("BF1c", "BF12") # name the columns of the output

for(i in 1:N){ # loop over persons
  BFtab <- results[[i]]$fit # obtain the fit table of person i
  # compute relevant BFs
  BF1c <- results[[i]]$fit[1,7]
  BF12 <- results[[i]]$BFmatrix[1,2]
  # save the 4 BFs in the i-th row of the output matrix
  output[i,] <- c(BF1c,BF12)
}

output # view the final output

###########################
##### gpbf function #######
###########################
gPBF <- function(BFs){
  N  <- ifelse(is.null(nrow(BFs)), length(BFs), nrow(BFs))
  
  res <- apply(BFs, 2, function(x){
    GP <- prod(x) ^ (1 / N)
    ER <- abs((GP < 1) - sum(x > 1)/N)
    SR <- ifelse(GP < 1, 
                 sum(x < GP) / N, 
                 sum(x > GP) / N)
    c(GP, ER, SR)
  })
  
  rownames(res) <- c("Geometric Product", "Evidence Rate", "Stability Rate")
  out <- list("GPBF" = res, "BFs" = BFs, "N" = N)
  class(out) <- "gPBF"
  return(out)
}
##############################
### obtaining the gpbfs ######
##############################
gpout <- gPBF(output)
gpout

save(output, gpout, file ="results.RData")

