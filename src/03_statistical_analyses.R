## ---------------------------------------------------------------------------##
##                                                                            ##
##  Title: STATISTICAL ANALYSES                                               ##
##                                                                            ##
##  Description: Statistical analyses to test the relationship between        ##
##  morphological variability, environmental suitability and distance from    ##
##  the range edge throughout the geographic range of the Gray slender        ##
##  opossum, Marmosops incanus. Models describing the relationship between    ##
##  dependent-independent variables were selected using the Akaike criteria.  ##
##                                                                            ##
##  Created by: Alan Braz (@brazagm)                                          ##
##                                                                            ##
## ---------------------------------------------------------------------------##
setwd("./")

library(adehabitatHR)
library(geosphere)
library(psych)
library(quantreg)
library(raster)
library(rgdal)
library(rgeos)
library(scales)
library(spatstat)


# 1. FUNCTIONAL VARIABILITY  -----------------------------------------------####
## Import data
# Obs: functional traits are calculated as the linear measures without correction
# because the values were extracted from the columns with values calculated before correction (i.e. nas.gls, can.bam,...)
data <- read.csv("./data/data_size-corrected.csv") # size-corrected morphological data
str(data)

## Set the name of morphological traits
trait_names <- c("gls", "cb", "nas", "bbc", "zb", "poc", "ic", "lbasi", "bps", "bam", "can", "pl", "max", "sm", "mad", "alpcor", "larcon") # name of linear measures

# Remove measures which measurement error > 20% and have sexual dimorphism in shape
trait_names <- trait_names[!trait_names %in% c("poc","bps","bam","larcon")]

## Set the name of functional traits
fun_names <- c("size", "shape", "nas.gls", "can.bam", "pl.cb", "bam.pl", "sm.max", "alpcor.mad") # name of functional traits - please, consult Material & Methods for further details on the calculations of functional traits

## Data filtering
# Only adult specimens of Marmosops incanus, with sex identified and populations defined
y <- subset(data, subset = sp == "Marmosops incanus" & idade == "adulto" & pop != "" & medido == TRUE)
y$pop <- as.factor(y$pop) # 'pop' as factor
populations <- levels(y$pop)  # population names


## 1.1. CALCULATE MORPHOLOGICAL VARIABILITY  ####
# Morphological variability was calculated based on the multidimensional approach described in Boucher et al. (2013)

## Standardize morphological measures
# Standardize the species' traits to mean = 0 and sd = 1, using the entire sample
y_std <- subset(y, select = c("pop", "lon", "lat", "sexo", trait_names, fun_names)) # Select only essential columns
y_std[, c(trait_names, fun_names)] <- psych::rescale(y_std[, c(trait_names, fun_names)], mean = 0, sd = 1) # standardize

## Calculate the variability extent (sensu Boucher et al. 2013)
# Create a dataframe for data input
# The variability extent is calculated based on two diferent measures
# 'raw' = linear measures of the skull corrected by size
# 'fun' = functional traits based on ratios of linear measures (e.g. can:bam, nas:gls)
df <- data.frame(population = populations, extent.raw = NA, integr.raw = NA, extent.fun = NA)

## Extract the variability of each population, for each type of measures (linear and functional trait)
for(i in populations){
  
  ## 1) Calculate the variability extension of linear measures
  yi <- subset(y_std, subset = pop == i, select = c(trait_names))
  yi <- na.omit(yi) # omit NAs
  df[which(df$population == i), "n.raw"] <- nrow(yi) # sample size
  
  # var-covar matrix for linear measures
  varcor <- cov(yi)
  df[which(df$population == i), "extent.raw"] <- sum(diag(varcor)) # extract the extent value
  
  
  ## 2) Calculate the variability integration of linear measures
  if(nrow(yi) > 1){
    df[which(df$population == i), "integr.raw"] <- var(eigen(varcor)$values)
  }
  
  ## 2) Calculate the variability extension of functional traits
  yi <- subset(y_std, subset = pop == i, select = c(fun_names))
  yi <- na.omit(yi) # omit NAs
  df[which(df$population == i), "n.fun"] <- nrow(yi) # sample size
  
  
  # var-covar matrix for functional traits
  varcor <- cov(yi)
  df[which(df$population == i), "extent.fun"] <- sum(diag(varcor)) # extract the extent value
}


# 2. PREDICTOR VARIABLES  -------------------------------------------------####
## 2.1. IMPORT ENVIRONMENTAL SUITABILITY  ####
## Import niche models
# Ecological niche models were calibrated according to Braz et al. (2020) "Interspecific competition constrains local abundance in highly suitable areas". Ecography 43: 1560–1570
list <- list.files("./results/niche_models", pattern = "avg.asc", full.names = TRUE)
sdm <- stack(list) # import rasters
names(sdm) <- gsub("_avg", "", names(sdm)) # remove '_avg' from raster names
projection(sdm) <- crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

## Rescale suitability scores ranging from 0 to 1 for each model
for(i in names(sdm)){
  if(i == "mahal"){
    x <- values(sdm[[i]])
    x <- -(x - 1) # Distance value
    x <- sqrt(x)
    x <- scales::rescale(x, to = c(0, 1), from = range(x, na.rm = TRUE, finite = TRUE))
    x <- -(x -1)
    values(sdm[[i]]) <- x
  }else{
    x <- values(sdm[[i]])
    x <- scales::rescale(x, to = c(0, 1), from = range(x, na.rm = TRUE))
    values(sdm[[i]]) <- x
  }
}


## 2.2. EXTRACT SUITABILITY VALUES  ####
# Because some populations include two or more nearby localities, the same population can have more than one suitability value
# In this case, the value is the simple average between all these different values

# Calculate the climatic/suitability for each locality
for(i in populations){

  # Extract values of environmental suitability
  for(x in names(sdm)){
    
    coord <- y[which(y$pop == i), c("lon", "lat")]
    suit <- extract(sdm[[x]], coord)
    suit <- unique(suit)
    suit <- mean(suit)
    
    df[which(df$population == i), x] <- suit
    df[which(df$population == i), "lon"] <- coord[1,"lon"]
    df[which(df$population == i), "lat"] <- coord[1,"lat"]
  }
}


## 2.3. DISTANCE FROM THE RANGE EDGE  ####
## Import occurrence records
# Records are filtered by a minimum distance of 10 km between them to avoid spatial bias in sample
points <- read.csv("./data/data_records_thin1.csv") # import csv
points <- SpatialPoints(points[, c("lon","lat")], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) # as spatialpoints

# Create a Minimum Convex Polygon (MCP)
limits <- readOGR("./data/shapes/borders/", layer = "TM_WORLD_BORDERS-0.3")
mcp <- mcp(points, percent = 100)
intersect <- intersect(limits, mcp)
df_points <- SpatialPointsDataFrame(df[,c("lon","lat")], df)
dist_line <- dist2Line(df_points, intersect) # calculate distance from points to MCP
df_points@data$distance <- dist_line[,"distance"]
df <- df_points@data


# 3. STATISTICAL ANALYSES  -------------------------------------------------####
## 3.1. TEST THE VARIABILITY-SUITABILITY RELATIONSHIP  ####
# Select predictors
pred_names <- c("domain","mahal","maxent", "distance")

# Data filtering
xy <- subset(df, subset = n.raw >= 5) # Only populations with sample size >= 5 individuals

# Something gone wrong: Catas Altas and Reserva Vale have NA values for suitabilities. Add values manually
xy[which(xy$population == "Catas Altas"), c("domain","mahal","maxent")] <- c(0.8806858, 0.9085349, 0.9675414)
xy[which(xy$population == "Reserva Vale"), c("domain","mahal","maxent")] <- c(0.9235663, 0.9540392, 0.5223794)

# Remove outliers (source from https://www.r-bloggers.com/identify-describe-plot-and-remove-the-outliers-from-the-dataset/)
# Explaratory analyses on residual distribution
source("http://goo.gl/UUyEzD")
outlierKD(xy, extent.raw) # YES: Remove 1 outlier from 26 observations


### 3.1.1 Test for mean differences and variance equality in variability between low/high populations ####
# Classify each locality as low/high suitability for each ENM
for(i in pred_names){
  
  # Low/high suitability was based on the median
  suit <- which(xy[,i] > median(xy[,i]))
  unsuit <- which(xy[,i] <= median(xy[,i]))
  
  # Insert a column in 'incanus' for low/high suitability for each ENM
  j <- paste(i, "binary", sep = "_")
  
  xy[suit, j] <- "high"
  xy[unsuit, j] <- "low"
}

# Create a dataframe for results
results <- data.frame(row.names = c("levene.F", "levene.p", "welch.t", "welch.p", "anova.F", "anova.p"))

# Analyses for each predictor
for(i in pred_names){
  library(car) # required for Levene's test
  
  # Characterize each sample as low/high environmental suitability
  y_low <- xy[which(xy[, paste(i, "binary", sep = "_")] == "low"),"extent.raw"]
  y_high <- xy[which(xy[, paste(i, "binary", sep = "_")] == "high"),"extent.raw"]
  x <- xy[, paste(i, "binary", sep = "_")]
  x <- as.factor(x)
  y <- xy$extent.raw
  
  # Calculate Levene's test, Welch's test (t-test with unequal variances), and ANOVA
  lev <- leveneTest(y, group = x, center= mean)
  welch <- t.test(y_low, y_high, alternative = "two.sided", var.equal = FALSE)
  anova <- aov(y ~ x); anova <- summary(anova)
  
  # Compute the results
  results[c("levene.F","levene.p"),i] <- c(lev$`F value`[1], lev$`Pr(>F)`[1])
  results[c("welch.t","welch.p"),i] <- c(welch$statistic, welch$p.value)
  results[c("anova.F","anova.p"),i] <- c(anova[[1]][["F value"]][1], anova[[1]][["Pr(>F)"]][1])
}

# Export results
write.csv(results, "./results/04_median_groups.csv")


### 3.1.2. Test for linear relationship between variability and predictor variables ####
# Test if morphological variability data follows a gamma distribution using Kolmogorov-Smirnov test
gamma <- fitdistr(na.omit(xy$extent.raw), "gamma") 
ks.test(na.omit(xy$extent.raw), "pgamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

## Fit linear regressions 
for(i in pred_names){
  
  if(i == "distance"){
    x <- log(xy[,i]) # log-transformation for distance from the range edge
  }else{
    x <- xy[,i]
  }
  
  # Model fitting
  lm <- lm(xy$extent.raw ~ x) # fit linear models
  glm <- glm(xy$extent.raw ~ x, family = Gamma(link = "log")) # fit GLM gamma
  
  if(i != "distance"){
    lm_distance <- lm(xy$extent.raw ~ x + log(xy$distance))
    assign(paste('lm_dist', i, sep = "_"), lm_distance)
    glm_distance <- glm(xy$extent.raw ~ x + xy$distance, family = Gamma(link = "log"))
    assign(paste('glm_dist', i, sep = "_"), glm_distance)
  }
  
  # Assign models to an object
  assign(paste('lm', i, sep = "_"), lm)
  assign(paste('glm', i, sep = "_"), glm)
}

# Model selection results
lm_sel <- MuMIn::model.sel(lm_domain, lm_dist_domain, lm_mahal, lm_dist_mahal, lm_maxent, lm_dist_maxent, lm_distance, lm(xy$extent.raw ~ 1))

glm_sel <- MuMIn::model.sel(glm_domain, glm_dist_domain, glm_mahal, glm_dist_mahal, glm_maxent, glm_dist_maxent, glm_distance, glm(xy$extent.raw ~ 1))

## Test if the dependency of lm residuals to the explanatory variable (i.e. heteroskedasticity) using Breusch-Pagan test
lmtest::bptest(lm_domain)
lmtest::bptest(lm_mahal)
lmtest::bptest(lm_maxent)
lmtest::bptest(lm_distance)

## Although a visual inspection of the xy plot suggests a certain degree of heteroskedasticity, variability does not seems to follow a gamma distribution and Breusch-Pagan test does not reject the homoskedasticity
## For this reason, results were based on simple linear regressions
write.csv(lm_sel, "./results/05_model_selection.csv")


# 4. FIGURES  --------------------------------------------------------------####
## 4.1. Graphic parameters  ####
# Create a function to automate the procedures
myPlot <- function(x, x.name, y.name){
  par(mfrow = c(1, 1),
  oma = c(0, 0, 0, 0), mar = c(4, 1, 1, 4), mgp = c(1.8, 0.5, 0),
  cex.lab = 1.4, family = "sans", font.lab = 2)
  
  plot(xy$extent.raw ~ x, pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = x.name)
  axis(4)
  mtext(y.name, 4, line = 2, cex = 1.4, family = "sans", font = 2)
}

## 4.2. Figure 3 - Relationship between morphological variability and explanatory variables ####
# Figure 3A - Domain
pdf("./results/figures/Figure 3A.pdf",  width = 5, height = 4)

myPlot(x = xy$domain, x.name = "Environmental suitability", y.name = "Morphological variability")
abline(lm_domain, col = "black", lty = 1, lwd = 2)
legend("topleft", legend = "R² = 0.16
P = 0.04
n = 25", bty = "n", text.font = 1, adj = c(0,0), cex = 1, inset = c(0, 0.1))

dev.off()

# Figure 3B - Mahalanobis distance
pdf("./results/figures/Figure 3B.pdf",  width = 5, height = 4)

myPlot(x = xy$mahal, x.name = "Environmental suitability", y.name = "Morphological variability")

dev.off()

# Figure 3C - Maxent
pdf("./results/figures/Figure 3C.pdf",  width = 5, height = 4)

myPlot(x = xy$maxent, x.name = "Environmental suitability", y.name = "Morphological variability")

dev.off()

# Figure 4 - Distance from the range edge
pdf("./results/figures/Figure 4.pdf",  width = 5, height = 4)

myPlot(x = log(xy$distance), x.name = expression(bold("Distance log"[10]*"(m)")), y.name = "Morphological variability")
abline(lm_distance, col = "black", lty = 1, lwd = 2)
legend("topleft", legend = "R² = 0.23
P = 0.02
n = 25", bty = "n", text.font = 1, adj = c(0,0), cex = 1, inset = c(0, 0.1))

dev.off()


## 4.3. Figure for Supplementary Material - Relationships without outliers  ####
# Select predictors
pred_names <- c("domain","mahal","maxent", "distance")

# Data filtering
xy_outlier <- subset(df, subset = n.raw >= 5) # Only populations with sample size >= 5 individuals

# Something gone wrong: Catas Altas and Reserva Vale have NA values for suitabilities. Add values manually
xy_outlier[which(xy_outlier$population == "Catas Altas"), c("domain","mahal","maxent")] <- c(0.8806858, 0.9085349, 0.9675414)
xy_outlier[which(xy_outlier$population == "Reserva Vale"), c("domain","mahal","maxent")] <- c(0.9235663, 0.9540392, 0.5223794)

# Graphic parameters
pdf("./results/figures/Figure S2.pdf",  width = 6, height = 5)

par(mfrow = c(2, 2),
    oma = c(2, 2, 0, 0), mar = c(3, 1, 1, 1), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.1, family = "sans", font.lab = 2)

# Create a function to automate the procedures
myPlot <- function(x, x.name, y.name){
  plot(xy_outlier$extent.raw ~ x, pch = 1, cex = 1.2,
       ylab = NA, xlab = x.name)
  #mtext(y.name, 4, line = 2, cex = 1.4, family = "sans", font = 2)
}

# Domain
myPlot(x = xy_outlier$domain, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Domain", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = xy_outlier$mahal, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Mahalanobis", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))
abline(lm_mahal, col = "black", lty = 1, lwd = 2)
legend("topleft", legend = "R² = 0.14
P = 0.04
n = 26", bty = "n", text.font = 1, adj = c(0,0), cex = 1, inset = c(0, 0.1))

myPlot(x = xy_outlier$maxent, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Maxent", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = log(xy_outlier$distance), x.name = expression(bold("Distance log"[10]*"(m)")), y.name = "Morphological variability")

title(ylab = "Morphological variability", cex.lab = 1.1, line = 0.5, outer = TRUE)

dev.off()


## 4.4. Correlation between environmental suitability and distance from the range edge  ####
# Set graphic parameters
pdf("./results/figures/Figure S4.pdf",  width = 10, height = 4)

par(mfrow = c(1,3), oma = c(3, 3, 0, 0), mar = c(1, 1, 1, 1), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.4, family = "sans", font.lab = 2)

# Create a function to automate the procedures
myPlot <- function(x, x.name, y.name){
  plot(x ~ log(xy$distance), pch = 1, cex = 1.2,
       ylab = "", xlab = x.name)
}

# Domain
myPlot(x = xy$domain, x.name = "", y.name = "")
legend("bottomright", legend = "Domain", bty = "n", text.font = 2, adj = c(0,0), cex = 1.5, inset = c(0.05, 0))

# Mahalanobis distance
myPlot(x = xy$mahal, x.name = "", y.name = "")
legend("bottomright", legend = "Mahalanobis", bty = "n", text.font = 2, adj = c(0,0), cex = 1.5, inset = c(0.05, 0))


# Maxent
myPlot(x = xy$maxent, x.name = "", y.name = "")
legend("bottomright", legend = "Maxent", bty = "n", text.font = 2, adj = c(0,0), cex = 1.5, inset = c(0.05, 0))

title(ylab = "Environmental suitability", cex.lab = 1.4, line = 1, outer = TRUE)
title(xlab = expression(bold("Distance log"[10]*"(m)")), cex.lab = 1.4, line = 1, outer = TRUE)

dev.off()


# 5. SPATIAL AUTOCORRELATION   ---------------------------------------------####
library(ape)
library(MuMIn)
library(nlme)

# Get lon~lat coordinates
df <- data.frame(lon = xy$lon, lat = xy$lat)
df <- SpatialPointsDataFrame(xy[, c("lon","lat")], xy, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")) # spatial points

# Remove Mocambinho population (= outlier sample)
df@data <- df@data[!df@data$pop == "Mocambinho",]

# Set linear model names
lm_names <- c("lm_domain", "lm_mahal", "lm_maxent", "lm_distance", "lm_dist_domain", "lm_dist_mahal", "lm_dist_maxent")

# Get linear model residuals
for(i in lm_names){
  model <- get(i)
  
  df@data[,paste(i, "resid", sep = "_")] <- model$residuals
}

# Calculate Moran's I
sample_dist <- as.matrix(dist(df@data[,c("lon","lat")]))
sample_dist <- 1/sample_dist
diag(sample_dist) <- 0

table <- data.frame(model = lm_names, observed = NA, expected = NA, sd = NA, p.value = NA)

for(i in lm_names){
  data <- df@data[, paste(i, "resid", sep = "_")]
  
  moran <- Moran.I(data, sample_dist, alternative = "greater")
  
  table[which(table$model == i), "observed"] <- moran$observed
  table[which(table$model == i), "expected"] <- moran$expected
  table[which(table$model == i), "sd"] <- moran$sd
  table[which(table$model == i), "p.value"] <- moran$p.value
}

write.csv(table, "./results/06_spatial_autocorr.csv")
