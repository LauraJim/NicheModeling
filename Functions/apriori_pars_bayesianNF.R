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
        # Estimate mean and variance of every environmental variable
        mu.prior[j] <- tolran[j,1] + (tolran[j,2]-tolran[j,1])/2
        s.prior[j] <- ((tolran[j,2]-tolran[j,1])/nsd)^2
      } else { 
        warning("the maximum values should be greater than the minimum values")
      }
    }
    # Define covariance matrix for the a priori distribution of the mean
    sigma0 <- diag(s.prior)
    A <- chol2inv(chol(sigma0))    # precision matrix
    W <- alpha*A
    CholW <- chol(W)
    Winv <- chol2inv(CholW)
    
  }
  return(list(mu0 = mu.prior, A0 = A, W0 = Winv))
}

# Main: How to use "priorpar --------------

# Example: Threnetes ruckeri

# read matrix with tolerance range for two environmental conditions 
#  (these are made up example values) of a species
limits <- read.csv("./T_ruckeri_tolerances.csv")

priorpar(limits[,2:3],6,2)

