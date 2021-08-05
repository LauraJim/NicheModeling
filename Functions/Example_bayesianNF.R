# Function: "Bayesian model for the fundamental Niche construction"
# Laura Jimenez & Carola Franzen
# July 2021
#

# Description: ---------------------------------


# Main: How the functions work (examples) --------------------------------------

# Working directory and libraries --------------------------------------

# Read functions and libraries
source(".\\Functions\\apriori_pars_bayesianNF.R")
source(".\\Functions\\Nf_Model_functions.R")

library(ggplot2)
library(ggpubr)

# Example: Threnetes ruckeri (two environmental layers) ----------------

## Preparing the parameters

# Read the environmental variables and occurrence data
# Geo-referenced data (longitud,latitude) + Environmental observations (Comp1,Comp2)
# read background data
envall <- read.csv("Threnetes_ruckeri_M_GE.csv",header=T) 
# read occurrence data
data <- read.csv("Threnetes_ruckeri_occ_GE.csv",header=T)

# Fixing the species and environmental variables to work with
# Indicate which columns from the occurrence and background files contain the environmental variables
DefineSp( env = envall, data.sp = data, Comp1=c(3,3), Comp2=c(4,4))
# Now we have fixed: Comp1, Comp2, env.d, env.sp, n, m, N, Et

# specify species name
rotule <- "Threnetes_ruckeri"
# select species color for plots
spcol <- "royalblue3"

# read matrix with tolerance range for two environmental conditions of a species
#  (these are made up example values) 
limits <- read.csv("./T_ruckeri_tolerances.csv")

# apply the function priorpar to receive boundaries for a species
boundary <- priorpar(limits[,2:3],6,2)
mu0 <<- boundary$mu0
A0 <<- boundary$A0
CholA0 <- chol(A0)
Sigma0 <<- chol2inv(CholA0)    # precision matrix
W0 <<- boundary$W0

el<-ellipse::ellipse(x=boundary$Sigma0,centre = boundary$mu0,level=0.95)



# Define a valid interval for mu, depending on the rage of the environmental variables
b1 <- 10 
b2 <- 100
mu.lim <<- c(min(Et[,1],limits[1,2])-b1,max(Et[,1],limits[1,3])+b1,min(Et[,2],limits[2,2])-b2,max(Et[,2],limits[2,3])+b2)


# First Plot: GSpace and ESpace ----------------------
### Plot the geographical and environmental spaces with the occurrence of the species on top
x11()
par(mfrow=c(1,2))
## Geographical Space:
# Plot the geographical locations of reported presences of the species, on top of the global environmental variables
plot(envall[,1],envall[,2],pch=".",col=1,xlab="Longitude",ylab="Latitude",main="Geographical space")
points(data[,1],data[,2],pch=20,col=spcol)

## Environmental Space:
# Plot environmental variables of species data and the location of reported presences of the species on top
plot(env.d, pch=".", col=1, xlim=mu.lim[1:2], ylim=mu.lim[3:4], xlab="Mean Annual Temperature (°C *10)",
     ylab="Total Annual Precipitation (mm)",main="Environmental space")
points(env.sp, pch=20, col=spcol)
rect(xleft = limits[1,2], xright= limits[1,3], ybottom = limits[2,2], ytop = limits[2,3], border="gold",lwd=2)
legend("topleft",legend=c("Species presences:",rotule,"Tolerance ranges"),pch=c(20,NA,NA),
       bty="n",lty=c(0,0,1),col=c(spcol,"white","gold"),lwd=2)
# Plot the contour line in the environmental space:
lines(el,col="gold",lwd=2)

# For the precision matrix, it is a Wishart distribution with
# alpha <- 2 #shape parameter, default, do not move for now
# W <- alpha*A0  #scale matrix W.
W <- boundary$W0
CholW <- chol(W)
Winv <- chol2inv(CholW)
W
# The prior expected value for the precision matrix is E(A) = W


# in ggplot

# Geographic Space:
# create new data-frame with combined coordinates of random background points
# and occurrence. 
bckgrnd1 <- data.frame(longitude = envall[, 1], latitude = envall[, 2])
occ1 <- data.frame(longitude = data[, 1], latitude = data[, 2])
data1 <- cbind(rbind(bckgrnd1[,1:2],occ1[,1:2]),
              # an extra column is added to differentiate bckgrnd and E.occ 
              # (1 for bckgrnd, 2 for E.occ)
              c(rep(1,nrow(bckgrnd1)),rep(2,nrow(occ1))))
# rename columns
data2 <- data.frame(Longitude = data1[, 1], Latitude = data1[, 2], 
                    Type = data1[,3])

# create plot for G-Space
p1 <- ggplot(data2, aes(x = Longitude, y = Latitude, color = factor(Type), 
                        shape = factor(Type))) +
  geom_point(alpha=0.5) +
  scale_shape_manual(values=c(20,19), guide = FALSE) +
  scale_color_manual(values= c("grey45", spcol), guide = FALSE)


# Environmental Space:
bckgrnd2 <- data.frame(Temperature = envall[,3], Precipitation = envall[,4])
species <- data.frame(Temperature = data[,3], Precipitation = data[,4])
data3 <- cbind(rbind(bckgrnd2[,1:2], species[,1:2]), c(rep(1,nrow(bckgrnd2)),rep(2,nrow(species))))
data4 <- data.frame(Temperature = data3[, 1], Precipitation = data3[, 2], 
                    Type = data3[,3]) 

df.el <- as.data.frame(el)
colnames(df.el) <- c("Temperature", "Precipitation")

p2 <- ggplot(data4, aes(x = Temperature, y = Precipitation, xlab = "Mean Annual Temperature", ylab = "Total Annual Precipitation")) +
  geom_point(aes(color = factor(Type), shape = factor(Type)), alpha = 0.6) +
  scale_shape_manual(values=c(20, 19), guide = FALSE) +
  scale_color_manual(name= "Species: Threnetes ruckeri",
                     labels= c("Background", "Presence"),
                     values= c("grey45", spcol)) +
  theme(legend.position = c(.05, .95), # for x, value of 0 puts it  to the 
        # left side, value of 1 to the right, for y, value of 0 puts it to 
        # the bottom, # value of 1 puts it to the top
        legend.justification = c("left", "top")) +
  scale_x_continuous("Annual mean temperature (°C*10)") +
  scale_y_continuous("Annual precipitation (mm)") +
  geom_path(data = df.el, color = "goldenrod1", size = 1.2) +
  geom_rect(xmin = limits[1,2], xmax = limits[1,3], ymin = limits[2,2], 
            ymax = limits[2,3], color = "goldenrod1", alpha = 0, size = 1.2)

x11()
ggarrange(p1, p2, ncol = 2, nrow = 1)


## -------------
# Number of environmental variables in the study
 dd <- 2
 
 