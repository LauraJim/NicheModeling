# Laura Jim?nez
# Firs version: September, 2018
# Last modification: October, 2020

### R code to create a matrix with the following columns: long, lat, bio1, bio12 (or more variables)

### Input data:
### (1) a csv file with the species name in the first column and
### the geographic coordinates of the occurrence points (in columns 2 and 3),
### (2) as many raster files as environmental variables to be included in the study

# Load packages
library(dismo)
library(tools)

# Read files with environmental layers
bio1 <- raster("./ClimateData10min/bio1WH.asc")
bio12 <- raster("./ClimateData10min/bio12WH.asc")

# Plot the environmental layers to check that we upload the correct ones
x11()
plot(bio12)

# Create a single raster with as many layers as environmental variables
bios <- stack(bio1,bio12)

# Read the names of all the csv files, each one containing occurrences of a single species
filelst <- list.files("./02_MaxentModels/Occ_Eval_indep")
spnames <- file_path_sans_ext(filelst)

# For each of the species: read the geographic coordinates, extract its environmental
# values, create a matrix with all the information, and save it into a new csv file

# save the number of occurrences of each species
nocc <- vector(mode = "numeric")

for(i in 1:length(spnames)){
  # Read file with occurrence points
  occ <- read.csv(paste0("./02_MaxentModels/Occ_Eval_indep/",filelst[i]),header=T)# [,c(5,3,4)] for Eval_Indep
  # Extract environmental values for each occurrence point
  occ1 <- extract(bios,occ[,2:3])
  occ2 <- na.omit(cbind(occ,occ1))
  nocc[i] <- nrow(occ2)
  # change column names following same order than raster stack
  colnames(occ2) <- c(colnames(occ),"bio1","bio12")
  # Save the resulting matrix in a new csv file
  write.csv(occ2,file=paste0("./02_MaxentModels/Occ_Eval_indep_Bios/",spnames[i],"_bios.csv"),row.names = F)
}

# Save a file with species names and number of occurrences
write.csv(cbind(spnames,nocc),file="./02_MaxentModels/Occ_Eval_indep_Bios/Spnames_occ_numbers.csv",row.names = F)
