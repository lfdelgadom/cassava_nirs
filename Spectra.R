#------NIR SPECTRA STRUCTURE-------

# Load necessary libraries
library(readr)         # For reading CSV files
library(tidyverse)     # For data manipulation and visualization

# Read the CSV file containing NIR spectra data
ss <- read_csv("https://raw.githubusercontent.com/sdjbrown/NIR-spectra-structure/master/spectra-structure.csv")

# Define colors for the plot, with darker colors indicating greater molar absorptivity
priorityCols <- rev(c("#00626D", "#008B93", "#7BB1B6", "#C0D1D2", "#E2E2E2"))

# Define the wavelengths for the vertical lines in the plot
waveLine <- seq(1000, 3100, by = 100)

# Reverse the order of the spectra for plotting
ord <- max(ss$ord) - ss$ord

# Categorize the molar absorption values for color coding in the plot
absorb <- as.numeric(cut(log10(ss$molarAbsorption), 4)) + 1
absorb[is.na(absorb)] <- 1

# Set up the output image for the plot
png(file="NIR-spectra-structure-chart.png", height=11.69, width=16.54, units = "in", res = 300)

# Set up the plot margins
par(mar = c(5,0.5,4,0.5))

# Initialize an empty plot with specified x and y limits and labels
plot(c(450, 3200), range(as.numeric(ord)), xlab = "Wavelength (nm)", type = "n", yaxt = "n", xaxt = "n", bty = "n", ylab = NA, main = "Spectra-structure correlation chart for the near-infrared region.", sub = "Data from Goddu and Delker 1960. Darker colours indicate greater molar absorbtivity.")

# Add vertical lines at specified wavelengths
abline(v = waveLine, col = priorityCols[1])

# Customize tick marks on the x-axis
par(tcl = -0.2)
axis(1, at = seq(1000, 3100, by = 50), labels = FALSE, lwd = 1, lwd.ticks = 1)
par(tcl = -0.5)
axis(1, at = waveLine, labels = waveLine, lwd = 0, lwd.ticks = 2)

# Draw segments for each spectrum, color-coded by molar absorptivity
segments(ss$wavelengthStart, ord, ss$wavelengthEnd, col = priorityCols[absorb], lwd = 10)

# Add text labels to the plot
text(rep(975, length(unique(ord))), unique(ord), labels = ss$description1[match(unique(ord), ord)], pos = 2)

# Display the output image
file.show("NIR-spectra-structure-chart.png")
