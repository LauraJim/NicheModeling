# Function: "priorpar"
# Laura Jimenez & Carola Franzen
# July 2021
#

# Description: ---------------------------------


# Main: How the functions work (examples) --------------------------------------

# Working directory and libraries -----------------------------

# Read functions that create accumulation curves of occurrences
source(".\\Functions\\apriori_pars_bayesianNF.R")


# Example: Threnetes ruckeri (two environmental layers) ------

# read matrix with tolerance range for two environmental conditions of a species
#  (these are made up example values) 
limits <- read.csv("./T_ruckeri_tolerances.csv")

# apply the function
priorpar(limits[,2:3],6,2)


## plot

# Read the environmental variables and occurrence data
# Geo-referenced data (longitud,latitude) + Environmental observations (Comp1,Comp2)

# read background data
envall <- read.csv("Threnetes_ruckeri_M_GE.csv",header=T) 
head(envall)

# read occurrence data
data <- read.csv("Threnetes_ruckeri_occ_GE.csv",header=T)
head(data)

# specify species name
rotule <- "Threnetes_ruckeri"
# select species color
spcol <- "#FF00AAFF"

