# Function: "priorpar"
# Laura Jimenez & Carola Franzen
# July 2021
#

# Description ----------------------------------

# The function *priorpar* uses data of environmental tolerance limits (maximum
# and mininimum) of a species to determine the values of the parameters that 
# define the a priori distribution. 

## Parameters:
# tolran = a matrix that contains the minimum and maximum values of the tolerance 
#         limits of a species (two columns for min and max and as many rows as  
#         environmental variables are used)
# nsd = number of standard deviations covered by a tolerance range
# alpha = number of degrees of freedom in the a priori Wishart distribution

## Output:
# The output is a list containing a vector of length equals to the number of
# dimensions in the environmental space, and two squared matrices with that
# same number of rows/columns. The first two elements represent the a priori
# parameters of *mu* and the third element is the scale matrix that
# defines the a priori distribution of W0.

# Function's code: priorpar ------------------

priorpar <- function(tolran,nsd,alpha=2)
{
  if (ncol(tolran)==2){ 
    mu.prior <- vector(mode="numeric",length=nrow(tolran))
    s.prior <- mu.prior
    for (j in 1:nrow(tolran)) {
      if(tolran[j,2] > tolran[j,1]){
        # Estimate mean and variance for every environmental variable
        mu.prior[j] <- tolran[j,1] + (tolran[j,2]-tolran[j,1])/2
        s.prior[j] <- ((tolran[j,2]-tolran[j,1])/nsd)^2
      } else { 
        warning("the maximum values should be greater than the minimum values")
      }
    }
    # Fix value of mu0
    mu0 <<- mu.prior
    
    # Define covariance matrix for the a priori distribution of the mean
    sigma0 <- diag(s.prior)
    
    # Fix value of Choleski decomposition of covariance matrix
    CholSigma0 <<- chol(sigma0)   
    
    # Define the scale parameter of the Wishart distribution
    W <- alpha*chol2inv(CholSigma0)
    
    # Fix value of Choleski decomposition of the scale parameter
    CholW <<- chol(W)
  }
}

# Main: How to use "priorpar --------------

# Example: Threnetes ruckeri

# read matrix with tolerance ranges for two environmental variables
limits <- read.csv("./T_ruckeri_tolerances.csv")

# apply function using only the columns with the tolerance limits
priorpar(limits[,2:3],nsd=6)
# Note that nsd=6 means that the tolerance ranges will cover
# 6 standard deviations around the mean value, mu0

# after applying this function, the following objects became global variables
mu0
CholSigma0
CholW

# Change:
# write csv for accumulation curve comparison (needed for nicheG)
#write.csv(boundary,file=paste0("./Results/Threnetes_ruckeri","_bayesian.csv"),row.names = F)

# End
