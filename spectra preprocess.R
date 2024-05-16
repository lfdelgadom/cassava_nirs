# specify all the packages used
myPackages <- c("tidyverse", "readxl", "clhs","prospectr")

# define which packages are not installed in the current computer
notInstalled <- myPackages[!(myPackages%in%installed.packages()[ , "Package"])]
# install the missing packages
if(length(notInstalled)>0) install.packages(notInstalled)

# Loading packages
library(tidyverse)
library(readxl)
library(clhs)
library(prospectr)

# The following dataset is about phisycal data and wavelengths with NIRs 
CASSAVA_BASE <- read.csv("Data/Data_raw.csv", header = T, sep = ";")

names(CASSAVA_BASE)[1:14]

# put the spectra into a single dataframe
spec_Cassava <- CASSAVA_BASE[grep("X", names(CASSAVA_BASE), value = TRUE)] 

# remove the spectra from the current dataframe
CASSAVA_BASE <- CASSAVA_BASE[ , -which(names(CASSAVA_BASE) %in% grep("X", names(CASSAVA_BASE),
                                                                     value = TRUE))]

# add the spectra to the dataframe
CASSAVA_BASE$spc <- spec_Cassava 

#The dataset is now much easier to read.
names(CASSAVA_BASE)

CASSAVA_BASE$spc

# We would like to remove the X in front of the column names of the spectra
# wavelengths. This will make easier plotting in a later step.

# take each column name from the spectra dataset
oldNames <- grep("X", names(CASSAVA_BASE$spc), value = TRUE) 

wavelength <- as.numeric(grep("X", names(CASSAVA_BASE$spc), value = TRUE))

# remove the "X" and make a numeric vector
wavelength <- as.numeric(substring(grep("X", names(CASSAVA_BASE$spc), value = T), 2, 20))

# change the name of the columns of the spectra
colnames(CASSAVA_BASE$spc) <- wavelength

# display the first ten column names
colnames(CASSAVA_BASE$spc)[1:10]

#-----------------------------
#----Plotting the spectra-----
#-----------------------------

# spectra base on Dry matter content

library(RColorBrewer)
# take three colours from the "RdBu" colour scale
cols = brewer.pal(3, "RdBu")
# make a function which convert the three colours to a gradient of colours
pal = colorRampPalette(cols)
# order the values of clay from the dataframe to find their increasing order
order = findInterval(CASSAVA_BASE$DM, sort(CASSAVA_BASE$DM))

png(paste0("images/","raw_spectra.jpg"), width = 1080, height = 800, 
    units = "px", pointsize = 12)

# plot spectra
matplot(
  x = colnames(CASSAVA_BASE$spc), y = t(CASSAVA_BASE$spc),
  xlab = "Wavelength (nm)",
  ylab = "log 1/R",
  main = "NIRs in cassava",
  type = "l",
  lty = 1,
  col = pal(nrow(CASSAVA_BASE))[order]
)
legend('topleft',
       border = "black",
       bty='n',
       lty = c(1, 1),
       # cex = 0.9,
       col = pal(2),
       # box.lty = 0,
       legend = c("Low dry matter", "High dry matter")
)
dev.off()


#---------------------------------------
#-----------Noise Removal---------------
#---------------------------------------

#We plot a single spectrum to see whether we need to remove the noisy spectra at
#wavelengths from 400 to 408 nm and 2492 to 2500 nm.

plot(colnames(CASSAVA_BASE$spc), CASSAVA_BASE$spc[1, ],
     type = "l",
     ylab = "Absorbance", xlab = "Wavelength /nm",
     main = "NIR spectral cassava (noisy sections)",
     col = rgb(
       red = 0.5,
       green = 0.5,
       blue = 0.5, alpha = 1
     )
)
rect(
  xleft = 360,
  xright = 450,
  ybottom = 0.20,
  ytop = 0.35,
  border = "red",
  lwd = 1,
  lty = "dashed"
)
rect(
  xleft = 2400, xright = 2530,
  ybottom = 1.75, ytop = 1.9,
  border = "red", lwd = 1, lty = "dashed"
)


#------------------------------------------------------
#-------Trimming the absorbance spectra----------------
#------------------------------------------------------


CASSAVA_BASE$spcT <- trimSpec(spectra = CASSAVA_BASE$spc, 
                              wavlimits = range(408, 2492))

# Plotting the spectrum as before is the same, except that 
# now we need to specify the new wavelength sequence of the trimmed spectra

wavs <- seq(408, 2492, by = 2)

#---------------------------------------------
#-------ploting the trimmed spectra-----------
#---------------------------------------------

plot(wavs, CASSAVA_BASE$spcT[1, ], type = "l",
     ylab = "Absorbance", xlab = "Wavelength /nm",
     main = "NIR spectral cassava (trimmed spectra)",
     col = rgb(red = 0.5, green = 0.5, blue = 0.5, alpha = 1))

#-------------------------------------------------------------------
#------------------Moving Window Average----------------------------
#-------------------------------------------------------------------

# Each wavelength value is taken as the average of the neighboring wavelengths. 
# The original spectra are smoothed, which reduces the information content but also 
# the noise that it contains too.

# The number of wavelengths lost is obtained by (w−1)/2
# where w is the window size.

# add some random noise to simulate an example of noisy spectra
CASSAVA_BASE$spctNoisy <- CASSAVA_BASE$spcT + rnorm(ncol(CASSAVA_BASE$spcT), 0, 0.003)

# specify the window size
windowMa <- 5

# apply the moving average
CASSAVA_BASE$spctMa <- movav(X = CASSAVA_BASE$spctNoisy, w = windowMa)

# plot the noisy spectrum
plot(names(CASSAVA_BASE$spctNoisy[1,]), CASSAVA_BASE$spctNoisy[1,],
     main = "Noisy (red line) and smoothed (black line)",
     type = "l",
     ylab = "Absorbance", xlab = "Wavelength /nm",
     sub = "Moving window average of size 5",
     col = rgb(red = 1, green = 0, blue = 0, alpha = 1))

# add the smoothed spectrum (black line)
lines(names(CASSAVA_BASE$spctMa[1,]), CASSAVA_BASE$spctMa[1, ],
      col = rgb(red = 0.3, green = 0.3, blue = 0.3, alpha = 1))

# add a legend
legend(400, 2, border = "white",
       cex=0.9,
       box.lty = 0,
       legend = c("Raw", "Moving average"),
       lty = c(1, 1),
       col = 2:1)

#----------------------------------------------
#--------------Savitzky-Golay Filtering--------
#----------------------------------------------

# Essentially this function uses the sgolayfilt function from the signal
# package. It uses the aaply function from the plyr package to apply the filter
# over the entire spectra collection


#-----------------------------------
#---------------Derivate 0----------

# filter spectra dataset
# comparar con la el spectro original adicionando un eje adicional a la derecha

# k = polynomial order, 
# w = filter length (odd)	
# m = m-th derivative

# filter + 1st derivative
CASSAVA_BASE$spctSg <- filterSg(CASSAVA_BASE$spcT,
                                w = 11,
                                k = 2,
                                m = 1)

# filter + 2nd derivative
CASSAVA_BASE$spctSg2 <- filterSg(CASSAVA_BASE$spcT,
                                 w = 11,
                                 k = 2,
                                 m = 2)


matplot(names(CASSAVA_BASE$spctSg), t(CASSAVA_BASE$spctSg),
        main = "Savitzky-Golay Filtering",
        type = "l",
        ylab = "Absorbance", xlab = "Wavelength /nm",
        sub = "Savitzky-Golay Filtering",
        col = rgb(red = 1, green = 0, blue = 0, alpha = 1))


# add the second derivative spectrum
matlines(names(CASSAVA_BASE$spctSg2), t(CASSAVA_BASE$spctSg2),
         col = rgb(red = 0.3, green = 0.3, blue = 0.3, alpha = 1))

# add a legend
legend('topleft',
       border = "black",
       bty='n',
       box.lty = 0,
       legend = c("1st derivative", "2nd derivative"),
       lty = c(1, 1),
       col = 2:1)

#----------------------------------------------------
#--------------Scatter Correction SNV ---------------
#----------------------------------------------------

# Standard Normal Variate

# apply the standard normal variate transformation
CASSAVA_BASE$specSnvC <- snvBLC(CASSAVA_BASE$spcT)


# plot the scatter-corrected spectra (contains negative values)
matplot(names(CASSAVA_BASE$specSnvC), t(CASSAVA_BASE$specSnvC),
        type = "l",
        ylab = "Absorbance", xlab = "Wavelength /nm",
        main = "Standard Normal Variate (SNV)",
        col = rgb(red = 1, green = 0, blue = 0, alpha = 1))

# plot the original spectra
matlines(names(CASSAVA_BASE$spcT), t(CASSAVA_BASE$spcT),
         type = "l",
         col = rgb(red = 0.3, green = 0.3, blue = 0.3, alpha = 1))

# add a legend
legend('topleft',
       border = "black",
       bty='n',
       box.lty = 0,
       legend = c("SNV", "Raw"),
       lty = c(1, 1),
       col = 2:1)

#-------------------------------------------------------
#------------Multiple scatter plot correction-----------
#-------------------------------------------------------

# function for applying multiplicative scatter correction
# apply the multiplicative scatter correction
CASSAVA_BASE$specMscC <- mscBLC(CASSAVA_BASE$spcT)

# plot the multiplicative scatter-corrected spectra
plot(names(CASSAVA_BASE$specMscC[1,]), CASSAVA_BASE$specMscC[1,],
     type = "l",
     ylab = "Absorbance", xlab = "Wavelength /nm",
     main = "Multiplicative scatter-corrected spectra",
     col = rgb(red = 1, green = 0, blue = 0, alpha = 1))

# add the raw spectra to the plot
lines(names(CASSAVA_BASE$spcT[1,]), CASSAVA_BASE$spcT[1,],
      col = rgb(red = 0.3, green = 0.3, blue = 0.3, alpha = 1))

# add a legend
legend('topleft',
       border = "black",
       bty='n',
       box.lty = 0,
       legend = c("MSC", "Raw"),
       lty = c(1, 1),
       col = 2:1)


#-----------------------------------------------------------------
#--------------------MSC, SNV, Raw----------------------------------
#-----------------------------------------------------------------

# plot the multiplicative scatter-corrected spectra (red)
matplot(names(CASSAVA_BASE$specMscC), t(CASSAVA_BASE$specMscC),
        type = "l",
        ylab = "Absorbance", xlab = "Wavelength /nm",
        main = "MSC, SNV and raw spectra",
        col = rgb(red = 1, green = 0, blue = 0, alpha = 1),
        ylim = c(-2, 2))

# add the SNV spectra to the plot (black)
matlines(names(CASSAVA_BASE$specSnvC), t(CASSAVA_BASE$specSnvC),
         col = rgb(red = 0.3, green = 0.3, blue = 0.3, alpha = 1),
         lty = 1)


# add a legend
legend('topleft',
       border = "black",
       bty='n',
       box.lty = 0,
       legend = c("MSC", "SNV", "Raw"),
       lty = c(1, 1, 1),
       col = c("red", "black"))

# El MSC es un suavisado de la señal y el SNV es una correción del 
# ruido.

# El MSC es necesario cuando no tengo tanto ruido como liquido.
#------------------------------------------------------
#--------------------Detrending------------------------
#------------------------------------------------------

#A practical alternative to SNV and MSC is to remove the mean value or a linear
#trend from the spectra. This method is called detrending.

# detrend the spectra
CASSAVA_BASE$specDT <- detrendSpc(CASSAVA_BASE$spcT)

# plot the detrended spectra
matplot(names(CASSAVA_BASE$specDT), t(CASSAVA_BASE$specDT),
        type = "l",
        ylab = " ", xlab = "Wavelength /nm",
        col = rgb(red = 1, green = 0, blue = 0, alpha = 1),
        ylim=c(-0.3, 2))

# add the original spectra
matlines(names(CASSAVA_BASE$spcT), t(CASSAVA_BASE$spcT),
         col = rgb(red = 0.3, green = 0.3, blue = 0.3, alpha = 1))

# add the remove linear trend
matlines(names(CASSAVA_BASE$spcT[1, ]),
         t(CASSAVA_BASE$spcT[1, ] - CASSAVA_BASE$specDT[1, ]),
         lty = 5)
# add a legend
legend("topleft",
       legend = c("raw", "detrended", "linear trend"),
       lty = c(1, 1, 5), col = c(1, 2, 1))


# SNV + detrend

CASSAVA_BASE$spc_SNV_DT <- snvBLC(CASSAVA_BASE$spcT) %>% detrendSpc()

matplot(names(CASSAVA_BASE$spc_SNV_DT), t(CASSAVA_BASE$spc_SNV_DT),
        type = "l",
        ylab = " ", xlab = "Wavelength /nm",
        col = rgb(red = 1, green = 0, blue = 0, alpha = 1),
        main = "SNV + Detrend") 

#------------------------------------------------------
#--------------------Derivatives-----------------------
#------------------------------------------------------

# Converting the spectra to first- or second-order derivatives aims at accentuating the
#absorbance features contained in the spectra. Derivatives also remove both additive
#and multiplicative effects on the spectra

# The first- and second-order derivatives of the smoothed spectra are computed using
# the filterSg function, by specifying this time m = 1 or m = 2 to indicate
# whether we want to take the first- or second-order derivative of the smoothed
# spectra

library(plyr)
# take first derivative of spectra
CASSAVA_BASE$specDeriv1 <- filterSg(CASSAVA_BASE$spcT,
                                    w = 3,
                                    k = 2,
                                    m = 1)

# take second derivative of spectra
CASSAVA_BASE$specDeriv2 <- filterSg(CASSAVA_BASE$spcT,
                                    w = 3,
                                    k = 2,
                                    m = 2)


# plot the first order derivative spectra
matplot(names(CASSAVA_BASE$specDeriv1), t(CASSAVA_BASE$specDeriv1),
        type = "l",
        ylab = " ",
        xlab = "Wavelength /nm",
        main = "First order derivative | Second order derivative",
        col = rgb(red = 1, green = 0, blue = 0, alpha = 1))

# plot the second order derivative spectra
matlines(names(CASSAVA_BASE$specDeriv2), t(CASSAVA_BASE$specDeriv2),
         col = rgb(red = 0.3, green = 0.3, blue = 0.3, alpha = 1))

# add a legend

legend("topleft",
       cex=0.7,
       legend = c("First order derivative", "Second order derivative"),
       lty = c(1, 1),
       col = 2:1)

# SNV + detrend

CASSAVA_BASE$spc_SNV_DT_1 <- snvBLC(CASSAVA_BASE$spcT) %>% 
  detrendSpc() %>% filterSg(w = 3, k = 2, m = 1) %>% movav(w = 4)

#----------------------------------------------------------------------
#---------------------Centring and Standardizing-----------------------
#----------------------------------------------------------------------

# centre the spectra wavelengths
CASSAVA_BASE$specNorm <- scale(CASSAVA_BASE$spcT, center = TRUE, scale = FALSE)

# standardize the spectra wavelengths
CASSAVA_BASE$specSdt <- scale(CASSAVA_BASE$spcT, center = TRUE, scale = TRUE)


# We can now plot with their original spectra (Fig. 5.12).
# plot the original spectra
plot(colnames(CASSAVA_BASE$spcT[1,]), CASSAVA_BASE$spcT[1,],
     type = "l",
     ylab = " ",
     xlab = "Wavelength /nm",
     col = rgb(red = 1, green = 0, blue = 0, alpha = 1),
     ylim= c(-0.7, 2))

# add to the plot the centred spectra

matlines(names(CASSAVA_BASE$specNorm[1,]), CASSAVA_BASE$specNorm[1,],
         col = rgb(red = 0.3, green = 1, blue = 0.3, alpha = 1))

# add to the plot the standardized spectra
matlines(names(CASSAVA_BASE$specSdt[1,]), CASSAVA_BASE$specSdt[1,],
         col = rgb(red = 0.3, green = 0.3, blue = 1, alpha = 1))

# add a legend
legend("topleft",
       cex=0.7,
       legend = c("Raw", "Centred", "Standardized"),
       lty = c(1, 1, 1),
       col = 2:4)
