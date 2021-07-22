# Laura Jimenez & Jorge Soberon
# First version: May, 2018

# Function 1: "get.table"--------------------------------------------------------
# Description:
# The function "get.table" will extract data from a raster that has information
# on suitable niches for a species and create a table with the extracted information.
#
# Parameters:
# G.occ = a matrix with three columns, where the second and third column contain
#         the coordinates in longitude and latitude
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

# END