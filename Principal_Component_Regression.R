# Load necessary functions and packages
source("https://raw.githubusercontent.com/lfdelgadom/cassava_nirs/main/functions.R")
library(here)

# Load transformed spectra data
sel_file_use <- list.files(here("output"))
load(here("output", sel_file_use))

# Apply SNV + detrend + 1st derivative
CASSAVA_BASE$spc_SNV_DT_1derived <- filterSg(
  CASSAVA_BASE$spc_SNV_DT,
  w = 3,
  k = 2,
  m = 1
)

# Split data into calibration (75%) and validation (25%) sets
set.seed(19101991)
calId <- sample(1:nrow(CASSAVA_BASE), size = round(0.75 * nrow(CASSAVA_BASE)))
datC <- CASSAVA_BASE[calId, ]
datV <- CASSAVA_BASE[-calId, ]

# Plot histograms for Dry Matter (DM) content in calibration and validation sets
png(paste0("images/","hist.jpg"), width = 1080, height = 800, 
    units = "px", pointsize = 25)
par(mfrow = c(1, 2))
hist(datC$DM, main = "Calibration Set", xlab = "Dry Matter")
hist(datV$DM, main = "Validation Set", xlab = "Dry Matter")
dev.off()

# Principal Component Regression
pcspectra <- prcomp(datC$spc_SNV_DT_1derived, center = TRUE, scale. = TRUE)
v <- pcspectra$sdev^2
cumv <- 100 * cumsum(v) / sum(v)

# Plot cumulative percentage of variance explained by PCs
png(paste0("images/","scree.jpg"), width = 1080, height = 800, 
    units = "px", pointsize = 25)
plot(
  cumv[1:20],
  type = "b",
  xlab = "Principal Component",
  ylab = "Cumulative Variance Explained (%)"
)
dev.off()

# Specify number of components to use
ncp <- 12

# Fit linear model using PCs
sdata <- as.data.frame(pcspectra$x[, 1:ncp])
DM_pc_CModel <- lm(datC$DM ~ ., data = sdata)
summary(DM_pc_CModel)

# Predictions for calibration dataset
DM_pc_CPred <- predict(DM_pc_CModel, sdata)

# Predictions for validation dataset
pcspectraV <- predict(pcspectra, datV$spc_SNV_DT_1derived)
sdataNew <- as.data.frame(pcspectraV[, 1:ncp])
DM_pc_VPred <- predict(DM_pc_CModel, sdataNew)

# Plot observed vs. predicted DM for calibration and validation sets
png(paste0("images/","obs_vs_pred.jpg"), width = 1080, height = 800, 
    units = "px", pointsize = 25)
par(mfrow = c(1, 2))
plot(
  datC$DM, DM_pc_CPred,
  xlab = "Observed DM",
  ylab = "Predicted DM",
  xlim = c(30, 45),
  ylim = c(30, 45),
  pch = 16,
  main = "Calibration"
)
abline(0, 1)

plot(
  datV$DM, DM_pc_VPred,
  xlab = "Observed DM",
  ylab = "Predicted DM",
  xlim = c(30, 45),
  ylim = c(30, 45),
  pch = 16,
  main = "Validation"
)
abline(0, 1)
dev.off()

# Evaluate model performance
modelPerformance(datC$DM, DM_pc_CPred)
modelPerformance(datV$DM, DM_pc_VPred)

# do it again with water absorption trait


