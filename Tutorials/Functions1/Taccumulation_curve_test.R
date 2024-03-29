# Laura Jimenez & Jorge Soberon
# First version: May, 2018

# Function 1: "get.table"--------------------------------------------------------
# Description:
# The function "get.table" will extract data from a raster that has information
# on suitable niches for a species and create a table with the extracted 
# information.
#
# Parameters:
# G.occ = a matrix with two columns with the coordinates of a species occurrence
# suit.Estck = a rasterstack that contains calculated suitability areas for a 
#           species and environmental layers, such as temperature or precipitation
#
# Output:
# The output is a matrix with coordinates, values from applied models such as the
# Mahalanobis distance, and environmental values. The matrix will have three 
# columns for the coordinates and model values and as many columns as there are
# environmental layers.

get.table <- function(G.occ,suit.Estck){
  # Convert to points so we have long,lat values in the first two columns
  mat1 <- rasterToPoints(suit.Estck)
  # order of index, ask order of rows order (range of indexes)
  iord <- order(mat1[,3], decreasing = T)
  # create new matrix with new order by suitability rows (high to low)
  mat2 <- mat1[iord,]
  
  # Now repeat the previous steps with the occurrence points
  occ1 <- extract(suit.Estck,G.occ)
  occ2 <- na.omit(cbind(G.occ,occ1))
  # sort the values of vector
  # order of index, ask order of rows order (range of indexes)
  iord2 <- order(occ2[,3], decreasing = T) 
  occ3 <- occ2[iord2,]
  colnames(occ3) <- colnames(mat2)
  mat3 <- rbind(occ3, mat2)
  mat4 <- cbind(mat3, Type=c(rep(1, nrow(occ3)), rep(0, nrow(mat2))))
  
  return(mat4)
}

# Function 2: "get.curve"--------------------------------------------------------
# Description:
# The function "get.curve" classifies and accumulates cells for the study area 
# and occurrences of a species based on the suitability that was previously 
# calculated by using a specific model, such as Mahalanobis distance. The 
# accumulated cells are then counted for the intervals of a curve.
#
# Parameters:
# occ.suit = a matrix with longitude and latitude of a species' occurrence and 
#             correlating suitability index, environmental conditions, and type
#             that has been classified and ordered 
# mod.suit = a matrix with longitude and latitude of the study area/subregions  
#             and correlating suitability index, environmental conditions, and
#             type that has been classified and ordered 
#
# Output:
# The first output gives the number of cells in the study are, the number of 
# occurrence points, and the probability of selecting an occurrence point with 
# the used model. The second output gives the number of accumulated number of 
# cells in the subregions and the accumulated number of occurrences in the 
# subregions. The second output provides the intervals between steps of a curve.


get.curve <- function(occ.suit,mod.suit){
  # Number of cells in study area
  nmod <- length(mod.suit)
  # Number of occurrences
  nocc <- length(occ.suit)
  # Prevalence = number of occurrences / number of cells in G
  preval <- nocc / (nmod + nocc)
  # Use vector of suitability values for the occurrences to classify all the cells
  # according to the partition of the interval (0,1) defined by that vector
  brks <- unique(occ.suit)
  if(brks[length(brks)] != 1){ # 1 is not included
    if(brks[1] != 0){ # 0 is not included
      mod.cls <- tabulate(cut(mod.suit,breaks=c(0,brks,1),labels=F,right=T))
    } else{ # 0 is included
      mod.cls <- tabulate(cut(mod.suit,breaks=c(brks,1),labels=F,right=T))
    }
  } else{ # 1 is included
    if(brks[1] != 0){ # 0 is not included
      mod.cls <- tabulate(cut(mod.suit,breaks=c(0,brks),labels=F,right=T))
    } else{ # 0 is included
      mod.cls <- tabulate(cut(mod.suit,breaks=brks,labels=F,right=T))
    }
  }
  # calculate the accumulated number of cells in the subregions
  mod.acc <- c(0,cumsum(rev(mod.cls)))
  # Count the number of occurrence points in each subregion
  counts <- vector(mode="numeric")
  dupl <- duplicated(occ.suit)*1 # ==1 if value is duplicated, ==0 otherwise
  for(i in 1:nocc){
    if(dupl[i]==0){ # a value can be replicated two or more times
      counts <- c(counts,1)
    } else {
      nn <- length(counts)
      counts[nn] <- counts[nn] + 1
    }
  }
  # calculate the accumulated number of occurrences in the subregions
  occ.acc <- c(0,cumsum(counts),rep(nocc,length(mod.acc)-length(counts)-1))
  # select values that contain important information
  last <- sum(occ.acc < nocc)
  if(mod.acc[last] == nmod){
    ntt <- 1:(last+1)
  } else{
    ntt <- 1:length(mod.acc)
  }
  
  return(list(out1=c(nmod,nocc,preval),out2=cbind(mod.acc[ntt],occ.acc[ntt])))
}

# Function 3: plot.aco
# Plots the accumulation curve of occurrences given a confidence level and the
# output of the function get.curve

plot.aco <- function(species,aco.curve,conlev,model){
  # Print the important values
  print(paste("Number of cells in study area:",aco.curve$out1[1]),quote=F)
  print(paste("Number of occurrence points:",aco.curve$out1[2]),quote=F)
  print(paste("Probability of selecting an occurrence point:",
              round(aco.curve$out1[3],4)),quote=F)
  
  # Calculate values of the confidence intervals around the random line using the null model
  if(conlev > 0){
    conlev1 <- (1 - conlev) / 2
    if(model == "binomial"){
      infs <- qbinom(conlev1,aco.curve$out2[,1],aco.curve$out1[3])
      sups <- qbinom(conlev1,aco.curve$out2[,1],aco.curve$out1[3],lower.tail = F)
    }
    if(model == "hypergeom") {
      infs <- qhyper(conlev1,m=aco.curve$out1[2],n=aco.curve$out1[1],
                     k=aco.curve$out2[,1])
      sups <- qhyper(conlev1,m=aco.curve$out1[2],n=aco.curve$out1[1],
                     k=aco.curve$out2[,1],lower.tail = F)
    }
  }
  nsub0 <- length(aco.curve$out2[,1])
  # Plot
  #x11()
  plot(aco.curve$out2[,1],aco.curve$out2[,1]*aco.curve$out1[3],type="b",col="red",
       xlab="Number of cells (ordered from most to least suitable)",
       ylab="Number of ccurrences",
       main="Accumulation of curve",xlim=c(0,aco.curve$out2[,1][nsub0]),
       ylim=c(0,aco.curve$out1[2]),lwd=2)
  # confidence intervals from hypergeom/binomial distribution
  if(conlev > 0){
    segments(aco.curve$out2[,1],infs,aco.curve$out2[,1],sups,col = "gray")
    points(aco.curve$out2[,1],infs,pch=19,col="grey25")
    points(aco.curve$out2[,1],sups,pch=19,col="grey25")
    if(model == "binomial") legmod <- paste("Binomial CI, p =",conlev)
    if(model == "hypergeom") legmod <- paste("Hypergeometric CI, p =",conlev)
  }
  # under non-random selection hypothesis
  lines(aco.curve$out2[,1],aco.curve$out2[,2],type="o",col="blue",lwd=2)
  if(nsub0<=50){
    text(aco.curve$out2[,1],aco.curve$out2[,2],labels=aco.curve$out2[,2],pos=2)
  } else {
    rind <- seq(1,nsub0,by=20) #%#
    text(aco.curve$out2[,1][rind],aco.curve$out2[,2][rind],
         labels=aco.curve$out2[,2][rind],pos=2)
  }
  legend("bottomright",legend=c(species,"Random counts",legmod,"SDM counts"),
         lwd=2,col=c("white","red","gray","blue"),bty="n")
}


### FUNCTION 4: accum.occ (v.3) ---------------------------------------------------

# ARGUMENTS:
# 'sp.name' -> character chain with species name, used for plot legends
# 'mod.Ecoords' -> raster with raw output from SDM
# 'G.occ' -> csv file with 3 columns and as many rows as presences of the species.
#               First column contains the name of the species and columns 2:3 contain the lon,lat coordinates.
# 'null.mod' -> indicate the distribution used as null model, possible values are "binomial" and "hypergeom".
# 'conlev' -> probability that indicates the confidence level to be drawn around the null model.
#             When 'flag=F', 'conlev' must be equal to 0.
# 'flag' -> if FALSE, values are caluculated and returned but not plotted

## Difference between 'accum.occ' and 'accum.occ1': 'mod.Ecoords' and 'occ.pnts' format is different, and,
## the second one don't produce plots in environmental space.

accum.occ <- function(names,G.occ,suit.Estck,null.mod="hypergeom",clev=0,flag=T){
  #
  table <- get.table(G.occ,suit.Estck)
  mod.ord <- table[table$Type==0,]
  occ.ord <- table[table$Type==1,]
  #
  curve <- get.curve(occ.suit=occ.ord[,3],mod.suit=mod.ord[,3])
  # Number of cells in study area
  nmod0 <- curve$out1[1]
  # Number of occurrences
  nocc0 <- curve$out1[2]
  # Prevalence = number of occurrences / number of cells in G
  preval0 <- curve$out1[3]
  # Accumulated number of cells in the subregions
  mod.acc0 <- curve$out2[,1]
  # Accumulated number of occurrences in the subregions
  occ.acc0 <- curve$out2[,2]

  # Now make all the plots
  if(flag==T){
    # before making the plots, we will use shades of gray to identify the different subareas
    nsub <- length(mod.acc0)
    cols <- gray((0:nsub/nsub)) #zero indicates black, and one indicates white
    ci <- vector("numeric")
    for (i in 1:nsub){
      # black indicates highly suitable, and white indicates highly unsuitable
      if(i==nsub){
        c <- rep(cols[nsub],nmod0-length(ci))
        ci <-c(ci,c)
      } else{
        c <- rep(cols[i],mod.acc0[i+1]-mod.acc0[i])
        ci <-c(ci,c)
      }
    }
    # we will use the world map from 'maptools'
    data("wrld_simpl", package="maptools")
    # but, before plotting we need to crop the world map using mod.Ecoords create the clipping polygon
    mod.ext <- bbox(SpatialPoints(mod.ord[,1:2]))
    ###
    # Plot 1: subregions in geographic space
    x11()
    plot(wrld_simpl,xlim=mod.ext[1,],ylim=mod.ext[2,],col="wheat1",axes=T,bg="azure2",
         main="Subregions in Geographical Space",xlab="Longitude",ylab="Latitude")
    # add points with corresponding gray shades
    points(mod.ord[,1:2],pch=15,col=ci,cex=0.5)
    # add occurrences
    points(occ.ord[,1:2],pch=19,col="red")
    ###
    # Plot 2: subregions in environmental space
    # NOTE: we plot first the points with low suitability (white/light grey) so they do not 
    # hide the points with high suitability (dark grey/black)
    cis <- c(rev(ci),rep(2,nocc0)) # grey shades for subregions + red for occurrence points
    if((ncol(table)-4)>2){
      nc <- 4:(ncol(table) - 1)
      # build a matrix with background and presence points, sorted for visualization
      plotm <- rbind(as.matrix(mod.ord[nmod0:1,nc]),as.matrix(occ.ord[,nc]))
      pch.occ <- ifelse(nocc0<=50,19,20)
      pcht <- c(rep(18,length(ci)),rep(pch.occ,nocc0))
      mypanel <- function(x,y,col=cis,pch=pcht,bgcolor="steelblue4"){
        # function used to color the background in the panels from function pairs
        ll <- par("usr")
        rect(ll[1], ll[3], ll[2], ll[4], col=bgcolor)
        points(x,y,col=col,pch=pch,bg=bgcolor)
      }
      x11()
      pairs(plotm,panel=mypanel,main="Subregions in Environmental Space")
    } else{
      if(length(names)==3){
        x11()
        plot(table[,4], table[,5], type="n", main="Subregions in Environmental Space", 
            xlab= names[2], ylab=names[3])
        u <- par("usr")
        rect(u[1], u[3], u[2], u[4], col = "steelblue4", border = "red")
        points(mod.ord[nmod0:1,4:5],pch=15,col=rev(ci),cex=0.5)
        # add occurrences
        points(occ.ord[,4:5],pch=19,col="red")
        legend("topleft",legend=c(names[1],"Occurence points","Suitability of points in M"),
               text.col="white",pch=c(19,19,15),col=c("steelblue4","red","grey"),bty="n")
      } else{
        print("The vector 'names' should have the species names and the names of the two environmental variables")
      }
    }
    ###
    # Plot 3: comparison among counts under random selection hypothesis
    x11()
    plot.aco(names[1],curve,clev,null.mod)
  }
  
  # Return coordinates of accumulation curve and corresponding percentages
  resul <- cbind(occ.acc0,mod.acc0,round((occ.acc0/nocc0)*100,2),round((mod.acc0/nmod0)*100,2))
  colnames(resul) <- c("No.occurrences","No.cells","%Gained Occ","%Area")
  return(resul)
}


### FUNCTION 5: COMPARING MODELS ------------------------------------------------

# Description: compare the occurrence-accumulation curves of different SDM/ENMs 
# for a single species

# ARGUMENTS:
# 'mods' -> list with as much elements as models to be compared, each element 
#           must be the resulting matrix of values from the function 'accum.occ' 
#           applied to the same occurrence data (first two columns only)
# 'nocc' -> number of occurrence points
# 'ncells' -> number of cells in M
# 'sp.name' -> character chain with the name of the species under study
# 'mods.names' -> character vector with the names of the SDMs to be compared, 
#                 must have the same length as 'mods'
# 'alpha' -> values between 0 and 1 representing the probability of the CI for the null model 

compare.aco <- function(mods,nocc,ncells,xrange=c(0,ncells),sp.name,mods.names,alpha){
  # number of models to compare
  nmods <- length(mods)
  if(length(mods.names)==nmods){
    # calculate the prevalence of the species in M
    preval <- nocc / ncells
    # Plot the curves of the models to be compared
    x11()
    plot(0:nocc,0:nocc,type="n",xlim=xrange,ylim=c(0,nocc),
         main=substitute(expr = italic(tl),env = list(tl=sp.name)),
         xlab="Number of cells (ordered from most to least suitable)",
         ylab="Number of ccurrences")
    # add the curves for each model and calculate the max value in x-axis
    colmod <- colorRampPalette(c("orange","royalblue"))(nmods)
    pt <- 15:(15+nmods)
    for (i in 1:nmods) {
      lines(mods[[i]][,2],mods[[i]][,1],type="o",pch=pt[i],lwd=2,col=colmod[i])
    }
    # add the line of random classification and confidence bands
    pnts <- floor(seq(0,xrange[2],length=nocc)) #ifelse(xrange[2]==ncells,ncells,max(xmax[2:nmods]))
    a1 <- (1-alpha)/2
    infs <- qhyper(a1,m=nocc,n=ncells-nocc,k=pnts)
    sups <- qhyper(1-a1,m=nocc,n=ncells-nocc,k=pnts)
    lines(0:ncells,0:ncells*preval,col="red",lwd=2)
    for (j in 1:nocc) {
      segments(pnts[j],infs[j],pnts[j],sups[j],col = "gray")
      points(pnts[j],infs[j],pch=20,col="grey25")
      points(pnts[j],sups[j],pch=20,col="grey25")
    }
    # add a legend to identify the different lines
    legend("bottomright",legend=c(mods.names,"Random counts","Hypergeometric-CI"),
           pch=c(pt,NA,NA),col=c(colmod,"red","gray"),lwd=3)
  } else{
    print("Warning! 'mods' and 'mods.names' should have the same length")
  }
}


## read libraries ---------------
library(dismo)
library(maptools)
library(rgdal)
library(raster)
# NOTE: package 'maptools' is also needed (used inside function)

# END