#
#
#

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
    A0 <- chol2inv(chol(sigma0))    # precision matrix
  }
  return(list(mu0 = mu.prior, Sigma0 = sigma0, W0 = alpha*A0))
}

# Example

limits <- read.csv("./T_ruckeri_tolerances.csv")

priorpar(limits[,2:3],6,2)

