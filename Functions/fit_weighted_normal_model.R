# Functions: rs.inE, negloglike
# Laura Jimenez
# First version: February 2020
# Last version: June 2021

# Description: ----------------------------------------------------------
# Maximum likelihood approach of the fundamental niche estimation problem using 
# a weighted distribution where the weights represent the availability of 
# environmental combinations inside M. This approach contains three functions.

## Parameters rs.inE:
# region = a shapefile of the study area (polygon)
# N = the sample size
# Estck = a rasterstack that contains at least two layers with environmental data

## Parameters negloglike
# guess = a vector of length 5 when d=2, it contains the mu and A values as elements
# sam1 = matrix containing the original sample of environmental combinations that 
#       correspond to presences
# sam2 = matrix containing a second random sample of environmental combinations 
#       which come from the area of study (M)

# rasterstack, shapefile and occurrence points

# Calling packages


# FUNCTIONS -------------------------------------------------------------

# Get a random sample of points inside the polygon that delimits M and extract their environmental values
## region should be a shapefile containing the polygon that defines the area of study
## N sample size
## Estck must be a stack of environmental raster layers
rs.inE <- function(region,N,Estck){
  # crop and mask the environmental layers with the polygon
  clip <- mask(crop(Estck,region),region)
  # get rid of cells with NA values = indices
  ind <- which(!is.na(clip[[1]][]))
  # get a random sample of indices
  sam <- sample(ind,N,replace = T)
  # choose the points corresponding to the selected indices
  pnts <- clip[][sam,]
  return(pnts)
}

# Negative log-likelihood function for theta=(mu,A)
## guess -- is a vector of length 5 when d=2, it contains the mu and A values as elements
## sam1 -- matrix containing the original sample of environmental combinations that correspond to presences
## sam2 -- matrix containing a second random sample of environmental combinations which come from the area of study (M)
negloglike <- function(guess,sam1,sam2){
  # define the parameters of interest from the guess parameter
  mu <- guess[1:2]
  A <- matrix(c(guess[3],guess[4],guess[4],guess[5]),nrow=2,ncol=2)
  # original sample size: number of presence points in the sample
  n <- nrow(sam1)
  # function that calculates quadratic terms, use inverse of matrix
  quad <- function(xi) { ax<-as.matrix(xi - mu); t(ax) %*% A %*% ax }
  q1 <- apply(sam1, 1, quad) # quadratic terms of presence points
  q2 <- apply(sam2, 1, quad) # quadratic terms of M points
  # negative log-likelihood value
  S <- 0.5*sum(q1) + n*log(sum(exp(-0.5*q2)))
  return(S)
}

# maximum likelihood
fitNiche <- function(E.occ, E.samM) {
  # calculate mu
  mu.ini <- colMeans(E.occ)
  # calculate A (covariance)
  Sig.ini <- cov(E.occ)
  # invert matrix sig.ini
  A.ini <- chol2inv(chol(Sig.ini))
  # whole vector of inicial values
  vals.ini <- c(mu.ini, A.ini[1,1], A.ini[1,2], A.ini[2,2])#c(mu.ini,A.ini[1,1],A.ini[1,2],A.ini[2,2])
  # fix the values of the samples used to evaluate the neg-log-likelihood
  like.fn <- function(theta){ negloglike(theta, sam1= E.occ, sam2= E.samM) } 
  find.mle <- optim(par=vals.ini, fn=like.fn, method="Nelder-Mead")
  mle <- find.mle$par
  mle.mu <- mle[1:2]
  mle.A <- matrix(c(mle[3:4],mle[4:5]),nrow=2,ncol=2)
  mle.Sig <- tryCatch(expr={chol2inv(chol(mle.A))}, error= function(e){NULL})

  # change column names for mle.Sig
  if(!is.null(mle.Sig)){
  colnames(mle.Sig) <- colnames(Sig.ini)
  rownames(mle.Sig) <- rownames(Sig.ini)
  }
    
  # wn = weighted normal distribution
  return(list(wn.mu = mle.mu, wn.sigma = mle.Sig, maha.mu = mu.ini, maha.sigma = Sig.ini))
}

# MAIN ------------------------------------------------------------------------

## Read libraries
# library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)
# library(maptools)
# data("wrld_simpl")

## Read datasets -------


# Read environmental layers
bio1 <- raster("./ClimateData10min/bio1WH.asc")
bio12 <- raster("./ClimateData10min/bio12WH.asc")
bios <- raster::stack(bio1,bio12)




## Set parameters -------

# Set the sample size (>= 10,000) that will be used to approximate expected value
N <- 10000



# Read presence points and M polygon
sp.occ <- read.csv("./Catasticta_nimbice_occ_GE.csv",header=T)[,-(1:2)]

M.shp <- readOGR("./Shapefiles","nimbice3")

# get a random sample of points in M and extract its corresponding environmental values
sam.Mpnts <- rs.inE(region = M.shp, N = N, Estck = bios)


# Looking for the MLE of mu and A --------------------------
# in tutorial add cache=TRUE to avoid running the function every time it knits
ml <- fitNiche(E.occ = sp.occ, E.samM = sam.Mpnts)

ml.table <- do.call("cbind",ml)
colnames(ml.table) <- c("wn.mu", "wn.sigma1", "wn.sigma2", "maha.mu", "maha.sigma1", "maha.sigma2")

write.csv(ml.table,"./Results/Catasticta_nimbice_Estimated_parameters.csv")

# get the ellipse defined by the ml estimators
el <- ellipse::ellipse(x=ml[[2]], centre=ml[[1]], level=0.99)
# get the ellipse from a multivarite normal model / mahalanobis distance method
el.ml <- ellipse::ellipse(x=ml[[4]], centre=ml[[3]], level=0.99)


# colorpalette
colpal <- c("grey70", "chartreuse4", "coral2", "cadetblue3")

# plot will be saved as .png
png(paste0("./Results/Catasticta_nimbice","_mle.png"),width = 2300, height = 2300, 
    res = 600, pointsize = 6)
x11()
plot(sam.Mpnts,col=colpal[1],pch=1, xlab="Annual mean temperature (°C*10)", 
     ylab="Annual precipitation (mm)")
# add presence points to the plot
points(sp.occ,col=colpal[3],pch=20,cex=1.5) # presences used in model
# ellipse maha
lines(el.ml,col=colpal[4],lwd=2)
# ellipse wn
lines(el,col=colpal[2],lwd=2)
sp.leg <- paste("Catasticta nimbice","(",nrow(sp.occ),")")
legend("topleft",legend = c(sp.leg,"Points inside M","Presences",
                            "Ellipse from Mahalanobis method",
                            "Ellipse from weighted-normal model"),
       pch=c(NA,1,19,NA,NA),col = c("white", colpal[1], colpal[3], colpal[4], colpal[2]),
       lwd=c(NA,NA,NA,2,2),bty = "n")
# finish saving png
dev.off()


## Example 2: Threnetes ruckeri in ggplot

# read input data
sp.occ2 <- read.csv("./Threnetes_ruckeri_occ_GE.csv",header=T)[,-(1:2)]

M.shp2 <- readOGR("./Shapefiles","Threnetes_ruckeri")

# get a random sample of points in M and extract its corresponding environmental values
sam.Mpnts2 <- rs.inE(region = M.shp2, N = 5000, Estck = bios)

# use function
ml2 <- fitNiche(E.occ = sp.occ2, E.samM = sam.Mpnts2)

# change function into proper table and rename column names
ml.table2 <- do.call("cbind",ml2)
colnames(ml.table2) <- c("wn.mu", "wn.sigma1", "wn.sigma2", "maha.mu", "maha.sigma1", "maha.sigma2")

# df.ml2 <- as.data.frame(ml.table2)

# write table as a csv
write.csv(ml.table2,"./Results/Threnetes_ruckeri_Estimated_parameters.csv")

# get the ellipse defined by the ml estimators
el2 <- ellipse::ellipse(x=ml2[[2]], centre=ml2[[1]], level=0.99, npoints = 500)

df.el2 <- as.data.frame(el2)
colnames(df.el2) <- c("Temperature", "Precipitation")



# get the ellipse from a multivarite normal model / mahalanobis distance method
el.ml2 <- ellipse::ellipse(x=ml2[[4]], centre=ml2[[3]], level=0.99)

df.elml <- as.data.frame(el.ml2)
colnames(df.elml) <- c("Temperature", "Precipitation")



## plot in ggplot
# prepare data as a dataframe for ggplot
bckgrnd <- data.frame(Temperature = sam.Mpnts2[,1], Precipitation = sam.Mpnts2[,2])
species <- data.frame(Temperature = sp.occ2[,1], Precipitation = sp.occ2[,2])
data <- cbind(rbind(bckgrnd[,1:2], species[,1:2]), c(rep(1,nrow(bckgrnd)),rep(2,nrow(species))))
data2 <- data.frame(Temperature = data[, 1], Precipitation = data[, 2], 
                    Type = data[,3]) 

# plot
x11()
##### doesn't work :(
ggplot(data2, aes(x = Temperature, y = Precipitation)) +
  geom_point(aes(color = factor(Type), shape = factor(Type))) +
  scale_shape_manual(values=c(1, 19), guide = FALSE) +
  scale_color_manual(name= "Data",
                     labels= c("Background", "Presence"),
                     values= c("1"=colpal[1], "2"= colpal[3])) +
  theme(legend.position = c(.05, .95), # for x, value of 0 puts it  to the 
        # left side, value of 1 to the right, for y, value of 0 puts it to 
        # the bottom, # value of 1 puts it to the top
        legend.justification = c("left", "top")) +
  scale_x_continuous("Annual mean temperature (°C*10)") +
  scale_y_continuous("Annual precipitation (mm)") +
  # geom_path(data = df_el2, lineend="butt", linejoin="round", linemitre=1, show.legend = "Ellipse from weighted-normal model")
  geom_path(data = df.el2, color = colpal[2], size = 1.2) +
  geom_path(data = df.elml, color = colpal[4], size = 1.2)



# END