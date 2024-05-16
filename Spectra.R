#------NIR SPECTRA STRUCTURE-------
library(readr)
library(tidyverse)
ss <- read_csv("https://raw.githubusercontent.com/sdjbrown/NIR-spectra-structure/master/spectra-structure.csv")

#Data from Goddu and Delker 1960. Darker colours indicate greater molar absorbtivity.  
#Create variables to assist in plotting
priorityCols <- rev(c("#00626D", "#008B93", "#7BB1B6", "#C0D1D2", "#E2E2E2"))
waveLine <- seq(1000, 3100, by = 100)
ord <- max(ss$ord) - ss$ord
absorb <- as.numeric(cut(log10(ss$molarAbsorption), 4)) + 1
absorb[is.na(absorb)] <- 1

png(file="NIR-spectra-structure-chart.png", height=11.69, width=16.54, units = "in", res = 300)
par(mar = c(5,0.5,4,0.5))
plot(c(450,3200), range(as.numeric(ord)), xlab = "Wavelength (nm)", type = "n", yaxt = "n", xaxt = "n", bty = "n", ylab = NA, main = "Spectra-structure correlation chart for the near-infrared region.", sub = "Data from Goddu and Delker 1960. Darker colours indicate greater molar absorbtivity.")
abline(v = waveLine, col = priorityCols[1])
par(tcl= -0.2)
axis(1, at = seq(1000, 3100, by = 50), labels = FALSE, lwd = 1, lwd.ticks = 1)
par(tcl= -0.5)
axis(1, at = waveLine, labels = waveLine, lwd = 0, lwd.ticks = 2)
segments(ss$wavelengthStart, ord, ss$wavelengthEnd, col = priorityCols[absorb], lwd = 10)
text(rep(975, length(unique(ord))), unique(ord), labels = ss$description1[match(unique(ord), ord)], pos = 2)
file.show("NIR-spectra-structure-chart.png")


# Prepare the data
ss <- ss %>% 
  mutate(
    ord = max(ord) - ord,
    absorb = cut(log10(molarAbsorption), breaks = 4, labels = FALSE),
    color = rev(c("#00626D", "#008B93", "#7BB1B6", "#C0D1D2", "#E2E2E2"))[absorb]
  )

# Create the plot
p <- ggplot(ss, aes(x = wavelengthStart, xend = wavelengthEnd, y = ord, yend = ord, color = color)) +
  geom_segment(size = 10) +
  scale_color_identity() +
  geom_text(data = unique(ss[, c("ord", "description1")]), aes(x = 975, y = ord, label = description1), hjust = 1) +
  theme_minimal() +
  labs(title = "Spectra-structure correlation chart for the near-infrared region",
       subtitle = "Data from Goddu and Delker 1960. Darker colours indicate greater molar absorbtivity.",
       x = "Wavelength (nm)", y = NULL) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(breaks = seq(1000, 3100, by = 100)) +
  coord_cartesian(xlim = c(1000, 3100))
