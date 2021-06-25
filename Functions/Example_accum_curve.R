# Laura Jim?nez
# April, 2021

# Working directory and libraries ----------------------------------------------
library(dismo)
library(maptools)
# NOTE: package 'maptools' is also needed (used inside function)

# Read functions that create accumulation curves of occurrences
source(".\\Rcode-hyperTest\\Occurrence_accum_curve.R")
source(".\\Rcode-hyperTest\\Occurrence_accum_curve1.R")
source(".\\Rcode-hyperTest\\compare_models.R")

# Species name
spname <- "Hamadryas_glauconome"
spna <- "Hamadryas_g"


# Maxent model -----------------------------------------------------------------
# Occurrence data with model output
spp.occSm <- read.csv(".\\Rcode-hyperTest\\Example\\Hamadryas_maxent_medium_occ.csv",header=T)
spp.occBm <- read.csv(".\\Rcode-hyperTest\\Example\\Hamadryas_maxent_big_occ.csv",header=T)

# Output for Ellipse Model for two study areas
spp.sm <- read.csv(".\\Rcode-hyperTest\\Example\\Hamadryas_maxent_medium_M.csv",header=T)
spp.bm <- read.csv(".\\Rcode-hyperTest\\Example\\Hamadryas_maxent_big_M.csv",header=T)

# Apply evaluation method
mod1.smallm <- accum.occ(spname,output.mod=spp.sm,occ.pnts=spp.occSm,null.mod="hypergeom",conlev=0.95)
mod2.bigm <- accum.occ(spname,output.mod=spp.bm,occ.pnts=spp.occBm,null.mod="hypergeom",conlev=0.95)

# Bioclim model ----------------------------------------------------------------
# Occurrence data with model output
spp.occSb <- read.csv(".\\Rcode-hyperTest\\Example\\Hamadryas_bioclim_medium_occ.csv",header=T)
spp.occBb <- read.csv(".\\Rcode-hyperTest\\Example\\Hamadryas_bioclim_big_occ.csv",header=T)

# Output for Ellipse Model for two study areas
spp.sb <- read.csv(".\\Rcode-hyperTest\\Example\\Hamadryas_bioclim_medium_M.csv",header=T)
spp.bb <- read.csv(".\\Rcode-hyperTest\\Example\\Hamadryas_bioclim_big_M.csv",header=T)

# Apply evaluation method
mod1.smallb <- accum.occ(spname,output.mod=spp.sb,occ.pnts=spp.occSb,null.mod="hypergeom",conlev=0.95)
mod2.bigb <- accum.occ(spname,output.mod=spp.bb,occ.pnts=spp.occBb,null.mod="hypergeom",conlev=0.95)

# Comparing the models ---------------------------------------------------------
# create a list with all the matrices using only the first two columns 
modelS <- list(mod1.smallm[,1:2], mod1.smallb[,1:2])
modelB <- list(mod2.bigm[,1:2], mod2.bigb[,1:2])

# Apply the function that draws all the models in a single plot
comp.accplot(mods=modelS,nocc=mod1.small[nrow(mod1.small),1],
             ncells=mod1.small[nrow(mod1.small),2],
             sp.name=spname,mods.names=c("Maxent","BioClim"),alpha=0.95)

# Apply same function but adjust the x-axis when using big M
comp.accplot(mods=modelB,nocc=mod2.big[nrow(mod2.big),1],
             ncells=mod2.big[nrow(mod2.big),2],sp.name=spname,
             mods.names=c("Maxent","BioClim"),alpha=0.95)

#
#
# This function is needed to create a table that is an occurrence table with an 
# additional column
#
# Create matrices "output.mod" for MaxEnt and BioClim --------------------------
get.table <- function(occ,columns,cnames,fname){
  # Convert to points so we have long,lat values in the first two columns
  mat1 <- rasterToPoints(columns)
  colnames(mat1) <- cnames
  head(mat1)
  
  # Save the resulting matrix
  write.csv(mat1,file=paste0(".\\Rcode-hyperTest\\Example\\",fname,"_M.csv"),row.names = F)
  
  # Now repeat the previous steps with the occurrence points
  occ1 <- extract(columns,occ[,2:3])
  occ2 <- cbind(occ[,2:3],occ1)
  colnames(occ2) <- c("long","lat",sdm.name,"bio1","bio12")
  head(occ2)
  
  # Save the resulting matrix
  write.csv(occ2,file=paste0(".\\Rcode-hyperTest\\Example\\",fname,"_occ.csv"),row.names = F)

  print("Done!")
}


# END