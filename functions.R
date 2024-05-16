# load R libraries necessary to use filterSg function
library(signal)
library(plyr)

#---------------------------------------------------------------------
#------------To build a function to remove these noisy portions-------
#---------------------------------------------------------------------

# The function requires two inputs:
#spectra (the dataset upon which we want to apply the spectral trimming) and
#wavlimits (the region of interest that we want to extract). The trimming function
#is called trimSpec.

# function for trimming spectra or targeting a specific spectral region of interest
trimSpec <- function(spectra, wavlimits){
  datawavs <- as.numeric(colnames(spectra))
  # set the limits
  limits <- which(datawavs %in% wavlimits)
  # mention the index that we keep from the matrix
  keptIndex <- seq(limits[1], limits[2], 1)
  # keep the index selected previously
  trimmedSpectra <- spectra[, keptIndex]
  # return the trimmed spectra
  keptNames <- datawavs[keptIndex]
  colnames(trimmedSpectra) <- keptNames
  return(trimmedSpectra)
}

#----------------------------------------------
#--------------Savitzky-Golay Filtering--------
#----------------------------------------------

# function for applying Savitzky-Golay smoothing filter
filterSg <- function(spectra, w, k, m) {
  spectra <- as.matrix(spectra)
  
  ## run filter, the window size in the sgolayfilt function is called n and
  ## the polynomial order is called p
  sg <- aaply(spectra, 1, sgolayfilt, n = w, p = k, m = m) # p = polynomial order, 
                                                          # n = filter length (odd)	
                                                          # m = m-th derivative

  ## arrange appropriately if a single sample
  if (nrow(spectra) == 1) {
    sg <- matrix(sg, dim(spectra))
  }
  ## return data frame
  sg <- as.data.frame(sg)
  colnames(sg) <- colnames(spectra)
  return(sg)
}

#----------------------------------------------------
#--------------------Scatter Correction--------------
#----------------------------------------------------

# Standard Normal Variate

# function for applying standard normal variate transformation
snvBLC <- function(spectra) {
  spectra <- as.matrix(spectra)
  snvMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
  
  # apply the standardization to each row
  for (i in 1:nrow(spectra)) {
    snvMat[i, ] <- (spectra[i, ] - mean(spectra[i, ]))/sd(spectra[i,])
  }
  snvMat<- as.data.frame(snvMat)
  colnames(snvMat) <- colnames(spectra)
  return(snvMat)
}

#-------------------------------------------------------
#------------Multiple scatter plot correction-----------
#-------------------------------------------------------

# Fearn (2008) recommends applying SNV scatter correction after
# filtering the noise from the spectra

# function for applying multiplicative scatter correction
mscBLC <- function(spectra) {
  
  # first calculate a mean spectrum.
  meanSpec <- as.matrix(colMeans(spectra))
  mscMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
  spectra <- as.matrix(spectra)
  
  # make a loop over each row
  for (i in 1:nrow(spectra)) {
    
    # determine the slope and intercept coefficients
    specLM <- lm(spectra[i, ] ~ meanSpec)
    specCE <- t(as.matrix(specLM$coefficients))
    
    # adjust the spectra
    mscMat[i, ] <- t(as.matrix((spectra[i, ] - specCE[1, 1])/specCE[1,2]))
  }
  mscMat <- as.data.frame(mscMat)
  colnames(mscMat) <- colnames(spectra)
  return(mscMat)
}

#------------------------------------------------------
#--------------------Detrending------------------------
#------------------------------------------------------

#A practical alternative to SNV and MSC is to remove the mean value or a linear
#trend from the spectra. This method is called detrending.

# function for detrending a matrix of spectra

detrendSpc <- function(spectra) {
  
  # load the required package
  require(pracma)
  
  detrendMat <- matrix(NA, ncol = ncol(spectra), nrow = nrow(spectra))
  spectra <- as.matrix(spectra)
  
  # make a loop over each row
  for (i in 1:nrow(spectra)) {
    
    # detrend each spectra, specify the linear model
    specLM <- pracma::detrend(spectra[i, ], tt = "linear")
    
    # take the values and store in the matrix
    detrendMat[i, ] <- as.numeric(specLM[,1])
  }
  detrendMat<- as.data.frame(detrendMat)
  colnames(detrendMat) <- colnames(spectra)
  return(detrendMat)
}

# reduce the column dimension of the spectra by averaging
compSpec <- function(spectra, w) {
  # ensure that the number of columns can be divided by the window size
  if(ncol(spectra)%%w != 0)
  {stop("Error: choose a compatible window size")
  }else{
    compMat <- matrix(NA, ncol = (ncol(spectra))/w, nrow = nrow(spectra))
    cc <- 1# loop over each column
    for (i in 1:ncol(compMat)) {
      compMat[, i] <- rowMeans(spectra[, cc:(cc + (w - 1))])
      cc <- cc + w
    }# provide the new column name (wavelength)
    colab = seq(as.numeric(colnames(spectra)[1]),
                as.numeric(colnames(spectra)[ncol(spectra)]), by = w)
    compMat <- as.data.frame(compMat)
    colnames(compMat) <- colab
  }
  return(compMat)
}

# Root mean square error
RMSE <- function(obs, pred){
  sqrt(mean((pred - obs)^2, na.rm = TRUE))
}


# Mean error (ME) or bias
ME <- function(obs, pred){
  mean(pred - obs, na.rm = TRUE)
}

# Squared correlation coefficient
r2 <- function(obs, pred){
  cor(pred, obs, method = "spearman", use = "pairwise.complete.obs")^2
}

# Coefficient of determination (R2)
R2 <- function(obs, pred){# sum of the squared error
  SSE <- sum((pred - obs) ^ 2, na.rm = T)# total sum of squares
  SST <- sum((obs - mean(obs, na.rm = T)) ^ 2, na.rm = T)
  R2 <- 1 - SSE/SST
  return(R2)
}

# Ratio of performance to deviation (RPD) 
RPD <- function(obs, pred){
  sdObs <- sd(obs)
  RMSE <- sqrt(mean((pred - obs)^2))
  rpd <- sdObs/RMSE
  return(rpd)
}

# Ratio of performance to inter-quartile distance (RPIQ) 
RPIQ <- function(obs, pred){
  q25 <- as.numeric(quantile(obs)[2])
  q75 <- as.numeric(quantile(obs)[4])
  iqDist <- q75 - q25
  RMSE <- sqrt(mean((pred - obs)^2))
  rpiq <- iqDist/RMSE
  return(rpiq)
}

modelPerformance <- function(obs, pred) {
  
  # Set options to avoid scientific notation temporarily
  old_options <- options(scipen = 999)
  
  # Calculate Root Mean Square Error
  rmse <- sqrt(mean((pred - obs)^2, na.rm = TRUE))
  
  # Calculate Mean Error
  me <- mean(pred - obs, na.rm = TRUE)
  
  # Calculate Squared Correlation Coefficient
  r2 <- cor(pred, obs, method = "spearman", use = "pairwise.complete.obs")^2
  
  # Calculate Coefficient of Determination
  sse <- sum((pred - obs) ^ 2, na.rm = TRUE)  # sum of squared errors
  sst <- sum((obs - mean(obs, na.rm = TRUE)) ^ 2, na.rm = TRUE)  # total sum of squares
  R2 <- 1 - sse / sst
  
  # Calculate Ratio of Performance to Deviation
  sd_obs <- sd(obs, na.rm = TRUE)
  rpd <- sd_obs / rmse
  
  # Calculate Ratio of Performance to Inter-Quartile Distance
  q25 <- quantile(obs, probs = 0.25, na.rm = TRUE)
  q75 <- quantile(obs, probs = 0.75, na.rm = TRUE)
  iq_dist <- q75 - q25
  rpiq <- iq_dist / rmse
  
  # Return all metrics as a named vector
  metrics <- c(RMSE = rmse, ME = me, r2 = r2, R2 = R2, RPD = rpd, RPIQ = rpiq)
  metrics <- round(metrics, 4)
  return(metrics)
}

