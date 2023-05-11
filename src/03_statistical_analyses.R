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
library(tidyverse)


# 0. DATA INPUT AND OUTPUT  ------------------------------------------------####
## 0.1. Inputs  ####
in_size_corrected <- "./data/data_size-corrected.csv"

## 0.2. Outputs  ####
out_population_data <- "./processed_data/03.1_populational_data.csv"
out_median <- "./results/03.2_median_groups.csv"
out_mean <- "./results/03.3_trait_mean_results.csv"
out_extent <- "./results/03.4_model_selection_extent.csv"
out_integr <- "./results/03.5_model_selection_integr.csv"

## 0.3. Set variable names  ####
## Set the name of morphological traits
trait_names <- c("gls", "cb", "nas", "bbc", "zb", "poc", "ic", "lbasi", "bps", "bam", "can", "pl", "max", "sm", "mad", "alpcor", "larcon") # name of linear measures

# Remove measures which measurement error > 20% and have sexual dimorphism in shape
trait_names <- trait_names[!trait_names %in% c("poc","bps", "bam","larcon")]


# 1. FUNCTIONAL VARIABILITY  -----------------------------------------------####
## Import size-corrected morphological data
data <- read.csv(in_size_corrected)

## Data filtering
# Only adult specimens of Marmosops incanus, with sex identification and populations location
y <- subset(data, subset = sp == "Marmosops incanus" & idade == "adulto" & pop != "" & medido == TRUE) %>%
  droplevels(.)

# Get population names
populations <- levels(droplevels(y$pop))


## 1.1. Calculate morphological variability  ####
# Morphological variability was calculated based on the multidimensional approach described in Boucher et al. (2013)

## Standardize morphological measures
y_std <- y %>%
  select(pop, lon, lat, sexo, ano, trait_names) # only essential columns
# Standardize the species' traits to mean = 0 and sd = 1, using the entire sample
y_std[, trait_names] <- psych::rescale(y_std[, trait_names], mean = 0, sd = 1)

## Create a dataframe for data input
df <- data.frame(population = populations, n = NA, extent = NA, integr = NA, year.min = NA, year.max = NA, temp.interval = NA, lon = NA, lat = NA)

## Extract the variability extent of each population
for(i in populations){
  
  # Get the individuals from population 'i'
  yi <- subset(y_std, subset = pop == i, select = c(trait_names)) %>%
    na.omit(.) # omit NAs
  
  # Get row index of population 'i'
  row_index <- which(df$population == i)
  
  # Get sample size of population 'i'
  df$n <- replace(df$n, row_index, nrow(yi))
  
  # Create var-covariation matrix
  varcor <- cov(yi)
  
  # For populations with more than one specimen...
  if(nrow(yi) > 1){
    ### 1.1.1. Extract the extent of morphological variability  ####
    df$extent <- replace(df$extent, row_index, sum(diag(varcor)))
    
    ### 1.1.2. Extract the morphological integration  ####
    df$integr <- replace(df$integr, row_index, var(eigen(varcor)$values))
  }
  
  # Record min/max year for the population
  yi_min <- min(subset(y_std, pop == i, ano))
  df$year.min <- replace(df$year.min, row_index, yi_min)
  yi_max <- max(subset(y_std, pop == i, ano))
  df$year.max <- replace(df$year.max, row_index, yi_max)
  
  # Record time interval for the population
  df$temp.interval <- replace(df$temp.interval, row_index, yi_max - yi_min)
  
  # Record lon/lat for the population
  yi_lon <- as.numeric(names(which.max(table(subset(y_std, pop == i, lon)))))
  df$lon <- replace(df$lon, row_index, yi_lon)
  yi_lat <- as.numeric(names(which.max(table(subset(y_std, pop == i, lat)))))
  df$lat <- replace(df$lat, row_index, yi_lat)
}


## 1.2. Calculate trait mean and variance  ####
for(t in trait_names){
  
  df[paste0(t, ".mean")] <- NA
  df[paste0(t, ".var")] <- NA
  
  for(i in populations){
    
    ### 1.2.1. Extract trait mean of each linear measures  ####
    y_mean <- subset(y, subset = pop == i, select = t)
    y_mean <- mean(na.omit(y_mean[,1])) # omit NAs
    df[which(df$population == i), paste0(t, ".mean")] <- y_mean
  
    ### 1.2.2. Extract trait variance of each linear measure  ####
    y_var <- subset(y, subset = pop == i, select = t)
    y_var <- var(na.omit(y_var[,1])) # omit NAs
    df[which(df$population == i), paste0(t, ".var")] <- y_var
    
    }
}


# 2. PREDICTOR VARIABLES  --------------------------------------------------####
## 2.1. Import environmental suitability rasters  ####
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


## 2.2. Extract suitability values  ####
# Because some populations include two or more nearby localities, the same population can have more than one suitability value
# In this case, the value is the simple average between all these different values

# Calculate the climatic/suitability for each locality
for(i in populations){

  # Extract values of environmental suitability
  for(x in names(sdm)){
    
    coord <- y[which(y$pop == i), c("lon", "lat")]
    suit <- raster::extract(sdm[[x]], coord)
    suit <- unique(suit)
    suit <- mean(suit)
    
    df[which(df$population == i), x] <- suit
    df[which(df$population == i), "lon"] <- coord[1,"lon"]
    df[which(df$population == i), "lat"] <- coord[1,"lat"]
  }
}


## 2.3. Extract distance from the range edge  ####
## Import occurrence records
# Records are filtered by a minimum distance of 10 km between them to avoid spatial bias in sample
points <- read.csv("./data/data_records_thin1.csv") # import csv
points <- SpatialPoints(points[, c("lon","lat")], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) # as spatialpoints

# Create a Minimum Convex Polygon (MCP)
limits <- readOGR("./data/shapes/", layer = "TM_WORLD_BORDERS-0.3")
mcp <- mcp(points, percent = 100)
intersect <- raster::intersect(limits, mcp)
df_points <- SpatialPointsDataFrame(df[,c("lon","lat")], df)
dist_line <- dist2Line(df_points, intersect) # calculate distance from points to MCP
df_points@data$distance <- dist_line[,"distance"]
df <- df_points@data

# export a dataframe with variability indices and predictor variables per population
write_csv(df, out_population_data)


# 3. STATISTICAL ANALYSES  -------------------------------------------------####
## 3.1. Data preparation and exploratory analysis  ####
# Set predictor variable names
pred_names <- c("domain","mahal","maxent", "distance")

# Data filtering
# Select only populations with sample size >= 5 individuals
xy <- subset(df, subset = n >= 5) 

# Something gone wrong: Catas Altas and Reserva Vale have NA values for suitabilities. Add values manually
xy[which(xy$population == "Catas Altas"), c("domain","mahal","maxent")] <- c(0.8806858, 0.9085349, 0.9675414)
xy[which(xy$population == "Reserva Vale"), c("domain","mahal","maxent")] <- c(0.9235663, 0.9540392, 0.5223794)

# Remove outliers (source from https://www.r-bloggers.com/identify-describe-plot-and-remove-the-outliers-from-the-dataset/)
# Explaratory analyses on residual distribution
source("http://goo.gl/UUyEzD")
outlierKD(xy, extent) # YES: Remove 1 outlier from 26 observations

# remove Mocambinho
xy <- subset(xy, population != "Mocambinho")

#### 3.1.1. Exploratory analyses  #####
## Variability extent data distribution
# seems to follow the Gaussian distribution
var <- xy$extent
hist(var, 15, freq = FALSE); lines(density(var))
shapiro.test(var)

## Morphological integration distribution
var <- xy$integr
hist(var, 15, freq = FALSE); lines(density(var))
shapiro.test(var)

# Test if morphological variability data follows a gamma distribution using Kolmogorov-Smirnov test
var <- xy$integr

gamma <- fitdistr(na.omit(var), "gamma") 
ks.test(na.omit(var), "pgamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


## 3.2. Are variability correlated with temporal interval?  ####
### 3.2.1. Correlation test for variability extent  ####
cor.test(xy$extent, xy$temp.interval, method = "pearson")
plot(xy$temp.interval, xy$extent)

### 3.2.2. Correlation test for morphological integration  ####
cor.test(xy$integr, xy$temp.interval, method = "pearson")
plot(xy$temp.interval, xy$integr)


### 3.3. Test for mean differences and variance equality in variability between low/high populations ####
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
  y_low <- xy[which(xy[, paste(i, "binary", sep = "_")] == "low"),"extent"]
  y_high <- xy[which(xy[, paste(i, "binary", sep = "_")] == "high"),"extent"]
  x <- xy[, paste(i, "binary", sep = "_")]
  x <- as.factor(x)
  y <- xy$extent
  
  # Calculate Levene's test, Welch's test (t-test with unequal variances), and ANOVA
  lev <- leveneTest(y, group = x, center = "mean")
  welch <- t.test(y_low, y_high, alternative = "two.sided", var.equal = FALSE)
  anova <- aov(y ~ x); anova <- summary(anova)
  
  # Compute the results
  results[c("levene.F","levene.p"),i] <- c(lev$`F value`[1], lev$`Pr(>F)`[1])
  results[c("welch.t","welch.p"),i] <- c(welch$statistic, welch$p.value)
  results[c("anova.F","anova.p"),i] <- c(anova[[1]][["F value"]][1], anova[[1]][["Pr(>F)"]][1])
}

# Export results
write.csv(results, out_median)


### 3.4. Test for linear relationship between mean and predictor variables ####
results_mean <- data.frame(row.names = c("domain_slope", "domain_r-squared", "domain_p-value", "dist_slope", "dist_r-squared", "dist_p-value", "both_domain_slope", "both_dist_slope", "both_r-squared", "both_p-values"))

for(t in paste0(trait_names, ".mean")){
  
  # Model fitting
  lm <- lm(xy[,t] ~ xy$domain) # fit linear models
  assign(paste('lm', t, sep = "_"), lm)
  results_mean[c("domain_slope","domain_r-squared","domain_p-value"),t] <- c(lm$coefficients[2], summary(lm)$coefficients[2,4], summary(lm)$r.squared)
  
  lm_distance <- lm(xy[,t] ~ log(xy$distance))
  assign(paste('lm_dist', t, sep = "_"), lm_distance)
  results_mean[c("dist_slope","dist_r-squared","dist_p-value"),t] <- c(lm_distance$coefficients[2], summary(lm_distance)$coefficients[2,4], summary(lm_distance)$r.squared)
  
  lm_both <- lm(xy[,t] ~ xy$domain + log(xy$distance))
  assign(paste('lm_both', t, sep = "_"), lm_both)
  results_mean[c("both_domain_slope","both_dist_slope","dist_r-squared","dist_p-value"),t] <- c(lm_both$coefficients[2], lm_both$coefficients[3], summary(lm_both)$coefficients[2,4], summary(lm_both)$r.squared)
  
}

# Export results
write.csv(results_mean, out_mean)


### 3.4. Test for linear relationship between variability and predictor variables ####
## Fit linear regressions and glm with gamma distribution
for(i in pred_names){
  
  if(i == "distance"){
    x <- log(xy[,i]) # log-transformation for distance from the range edge
  }else{
    x <- xy[,i]
  }
  
  for(resp in c("extent", "integr")){
  
    if(resp == "extent"){
      y <- xy$extent
    }else{
      y <- xy$integr
    }
    
    # Model fitting
    lm <- lm(y ~ x) # fit linear models
    glm <- glm(y ~ x, family = Gamma(link = "log")) # fit GLM gamma

    if(i != "distance"){
      lm_x_distance <- lm(y ~ x + log(xy$distance))
      assign(paste(resp, 'lm_dist', i, sep = "_"), lm_x_distance)
    
      glm_x_distance <- glm(y ~ x + xy$distance, family = Gamma(link = "log"))
      assign(paste(resp, 'glm_dist', i, sep = "_"), glm_x_distance)
      }
  
  # Assign models to an object
  assign(paste(resp, 'lm', i, sep = "_"), lm)
  assign(paste(resp, 'glm', i, sep = "_"), glm)
  }
}

## Replace integration ~ suitability models by gls with autocorrelation due to Moran's I test for spatial correlation
integr_lm_domain <- gls(integr~domain, xy, correlation = corExp(form = ~lon+lat, nugget = TRUE))
integr_lm_mahal <- gls(integr~mahal, xy, correlation = corExp(form = ~lon+lat, nugget = TRUE))
integr_lm_maxent <- gls(integr~maxent, xy, correlation = corExp(form = ~lon+lat, nugget = TRUE))

## Test if the dependency of lm residuals to the explanatory variable (i.e. heteroskedasticity) using Breusch-Pagan test
lmtest::bptest(extent_lm_domain)
lmtest::bptest(extent_lm_mahal)
lmtest::bptest(extent_lm_maxent)
lmtest::bptest(extent_lm_distance)

lmtest::bptest(integr_lm_domain)
lmtest::bptest(integr_lm_mahal)
lmtest::bptest(integr_lm_maxent)
lmtest::bptest(integr_lm_distance)

# Model selection results
extent_lm_sel <- MuMIn::model.sel(extent_lm_domain,
                                  extent_lm_dist_domain,
                                  extent_lm_mahal,
                                  extent_lm_dist_mahal,
                                  extent_lm_maxent,
                                  extent_lm_dist_maxent,
                                  extent_lm_distance,
                                  lm(xy$extent ~ 1))

integr_lm_sel <- MuMIn::model.sel(integr_lm_domain,
                                  integr_lm_dist_domain,
                                  integr_lm_mahal,
                                  integr_lm_dist_mahal,
                                  integr_lm_maxent,
                                  integr_lm_dist_maxent,
                                  integr_lm_distance,
                                  lm(xy$integr ~ 1))

extent_glm_sel <- MuMIn::model.sel(extent_glm_domain,
                                  extent_glm_dist_domain,
                                  extent_glm_mahal,
                                  extent_glm_dist_mahal,
                                  extent_glm_maxent,
                                  extent_glm_dist_maxent,
                                  extent_glm_distance,
                                  glm(xy$extent ~ 1, family = Gamma(link = "log")))

integr_glm_sel <- MuMIn::model.sel(integr_glm_domain,
                                  integr_glm_dist_domain,
                                  integr_glm_mahal,
                                  integr_glm_dist_mahal,
                                  integr_glm_maxent,
                                  integr_glm_dist_maxent,
                                  integr_glm_distance,
                                  glm(xy$integr ~ 1, family = Gamma(link = "log")))


## Although a visual inspection of the xy plot suggests a certain degree of heteroskedasticity, variability does not seems to follow a gamma distribution and Breusch-Pagan test does not reject the homoskedasticity
## For this reason, results were based on simple linear regressions
write.csv(extent_lm_sel, out_extent)
write.csv(integr_lm_sel, out_integr)


# 4. FIGURES  --------------------------------------------------------------####
## 4.1. Figure 3 - Relationship between morphological variability and explanatory variables ####
### Figure 3A - Variability ~ Domain
pdf("./results/figures/Figure 3A.pdf",  width = 8, height = 3.5)

par(mfrow = c(1, 2),
    oma = c(0, 0, 0, 0), mar = c(4, 1, 1, 4), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.4, family = "sans", font.lab = 2)

plot(xy$extent ~ xy$domain, pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = "Environmental suitability")
axis(4)
mtext("Morphological variability", 4, line = 2, cex = 1.4, family = "sans", font = 2)

abline(extent_lm_domain, col = "black", lty = 1, lwd = 2)
legend("topleft", legend = "R² = 0.17
P = 0.04
n = 25", bty = "n", text.font = 1, adj = c(0,0), cex = 1, inset = c(0, 0.1))


plot(xy$integr ~ xy$domain, pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = "Environmental suitability")
axis(4)
mtext("Morphological integration", 4, line = 2, cex = 1.4, family = "sans", font = 2)

dev.off()

# Figure 3B - Variability ~ Mahalanobis distance
pdf("./results/figures/Figure 3B.pdf",  width = 8, height = 3.5)

par(mfrow = c(1, 2),
    oma = c(0, 0, 0, 0), mar = c(4, 1, 1, 4), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.4, family = "sans", font.lab = 2)

plot(xy$extent ~ xy$mahal, pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = "Environmental suitability")
axis(4)
mtext("Morphological variability", 4, line = 2, cex = 1.4, family = "sans", font = 2)

plot(xy$integr ~ xy$mahal, pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = "Environmental suitability")
axis(4)
mtext("Morphological integration", 4, line = 2, cex = 1.4, family = "sans", font = 2)

dev.off()

# Figure 3C - Variability ~ Maxent
pdf("./results/figures/Figure 3C.pdf",  width = 8, height = 3.5)

par(mfrow = c(1, 2),
    oma = c(0, 0, 0, 0), mar = c(4, 1, 1, 4), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.4, family = "sans", font.lab = 2)

plot(xy$extent ~ xy$maxent, pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = "Environmental suitability")
axis(4)
mtext("Morphological variability", 4, line = 2, cex = 1.4, family = "sans", font = 2)

plot(xy$integr ~ xy$maxent, pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = "Environmental suitability")
axis(4)
mtext("Morphological integration", 4, line = 2, cex = 1.4, family = "sans", font = 2)

dev.off()

# Figure 4 - Distance from the range edge
pdf("./results/figures/Figure 4.pdf",  width = 8, height = 3.5)

par(mfrow = c(1, 2),
    oma = c(0, 0, 0, 0), mar = c(4, 1, 1, 4), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.4, family = "sans", font.lab = 2)

plot(xy$extent ~ log(xy$distance), pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = expression(bold("Distance log"[10]*"(m)")))
axis(4)
mtext("Morphological variability", 4, line = 2, cex = 1.4, family = "sans", font = 2)

abline(extent_lm_distance, col = "black", lty = 1, lwd = 2)
legend("topleft", legend = "R² = 0.23
P = 0.02
n = 25", bty = "n", text.font = 1, adj = c(0,0), cex = 1, inset = c(0, 0.1))

plot(xy$integr ~ log(xy$distance), pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = expression(bold("Distance log"[10]*"(m)")))
axis(4)
mtext("Morphological integration", 4, line = 2, cex = 1.4, family = "sans", font = 2)

abline(integr_lm_distance, col = "black", lty = 1, lwd = 2)
legend("topleft", legend = "R² = 0.21
P = 0.02
n = 25", bty = "n", text.font = 1, adj = c(0,0), cex = 1, inset = c(0, 0.1))

dev.off()


## 4.3. Figure for Supplementary Material - Relationships with outliers  ####
# Select predictors
pred_names <- c("domain","mahal","maxent", "distance")

# Data filtering
xy_outlier <- subset(df, subset = n >= 5) # Only populations with sample size >= 5 individuals

# Something gone wrong: Catas Altas and Reserva Vale have NA values for suitabilities. Add values manually
xy_outlier[which(xy_outlier$population == "Catas Altas"), c("domain","mahal","maxent")] <- c(0.8806858, 0.9085349, 0.9675414)
xy_outlier[which(xy_outlier$population == "Reserva Vale"), c("domain","mahal","maxent")] <- c(0.9235663, 0.9540392, 0.5223794)

### 4.3.1. Figure S1A - Extent of morphological variability  ####
# Graphic parameters
pdf("./results/figures/Figure S2A.pdf",  width = 6, height = 5)

par(mfrow = c(2, 2),
    oma = c(2, 2, 0, 0), mar = c(3, 1, 1, 1), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.1, family = "sans", font.lab = 2)

# Create a function to automate the procedures
myPlot <- function(x, x.name, y.name){
  plot(xy_outlier$extent ~ x, pch = 1, cex = 1.2,
       ylab = NA, xlab = x.name)
  #mtext(y.name, 4, line = 2, cex = 1.4, family = "sans", font = 2)
}

# Domain
myPlot(x = xy_outlier$domain, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Domain", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = xy_outlier$mahal, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Mahalanobis", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = xy_outlier$maxent, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Maxent", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = log(xy_outlier$distance), x.name = expression(bold("Distance log"[10]*"(m)")), y.name = "Morphological variability")

title(ylab = "Morphological variability", cex.lab = 1.1, line = 0.5, outer = TRUE)

dev.off()


### 4.3.2. Figure S1B - Morphological integration  ####
# Graphic parameters
pdf("./results/figures/Figure S2B.pdf",  width = 6, height = 5)

par(mfrow = c(2, 2),
    oma = c(2, 2, 0, 0), mar = c(3, 1, 1, 1), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.1, family = "sans", font.lab = 2)

# Create a function to automate the procedures
myPlot <- function(x, x.name, y.name){
  plot(xy_outlier$integr ~ x, pch = 1, cex = 1.2,
       ylab = NA, xlab = x.name)
  #mtext(y.name, 4, line = 2, cex = 1.4, family = "sans", font = 2)
}

# Domain
myPlot(x = xy_outlier$domain, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Domain", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = xy_outlier$mahal, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Mahalanobis", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = xy_outlier$maxent, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Maxent", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = log(xy_outlier$distance), x.name = expression(bold("Distance log"[10]*"(m)")), y.name = "Morphological integration")

title(ylab = "Morphological integration", cex.lab = 1.1, line = 0.5, outer = TRUE)

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
df <- SpatialPointsDataFrame(xy[, c("lon","lat")], xy, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs")) # spatial points

# Remove Mocambinho population (= outlier sample)
df@data <- df@data[!df@data$pop == "Mocambinho",]

# Set linear model names
extent_glm_names <- c("extent_glm_domain", "extent_glm_mahal", "extent_glm_maxent", "extent_glm_distance", "extent_glm_dist_domain", "extent_glm_dist_mahal", "extent_glm_dist_maxent")

integr_glm_names <- c("integr_glm_domain", "integr_glm_mahal", "integr_glm_maxent", "integr_glm_distance", "integr_glm_dist_domain", "integr_glm_dist_mahal", "integr_glm_dist_maxent")


## 5.1. Test for spatial autocorrelation using Moran's I index  #####
# Get linear model residuals for the extent models
for(i in extent_glm_names){
  model <- get(i)
  
  df@data[,paste(i, "resid", sep = "_")] <- model$residuals
}

# Calculate Moran's I
sample_dist <- as.matrix(dist(df@data[,c("lon","lat")]))
sample_dist <- 1/sample_dist
diag(sample_dist) <- 0

table <- data.frame(model = extent_glm_names, observed = NA, expected = NA, sd = NA, p.value = NA)

for(i in extent_glm_names){
  data <- df@data[, paste(i, "resid", sep = "_")]
  
  moran <- Moran.I(data, sample_dist, alternative = "greater")
  
  table[which(table$model == i), "observed"] <- moran$observed
  table[which(table$model == i), "expected"] <- moran$expected
  table[which(table$model == i), "sd"] <- moran$sd
  table[which(table$model == i), "p.value"] <- moran$p.value
}

# Get linear model residuals for the integration models
for(i in integr_glm_names){
  model <- get(i)
  
  df@data[,paste(i, "resid", sep = "_")] <- model$residuals
}

# Calculate Moran's I for the integration models
sample_dist <- as.matrix(dist(xy[,c("lon","lat")]))
sample_dist <- 1/sample_dist
diag(sample_dist) <- 0

table <- data.frame(model = integr_glm_names, observed = NA, expected = NA, sd = NA, p.value = NA)

for(i in integr_glm_names){
  data <- df@data[, paste(i, "resid", sep = "_")]
  
  moran <- Moran.I(data, sample_dist, alternative = "greater")
  
  table[which(table$model == i), "observed"] <- moran$observed
  table[which(table$model == i), "expected"] <- moran$expected
  table[which(table$model == i), "sd"] <- moran$sd
  table[which(table$model == i), "p.value"] <- moran$p.value
}

write.csv(table, "./results/06_spatial_autocorr_integr.csv")


## 5.2. Control spatial autocorrelation in models that rejected the Moran's test  ####
# Model selection to identify the best correlation function
# Test with domain, mahalanobis or maxent
model.sel(gls(integr ~ maxent, df, correlation = corExp(form = ~lon + lat, nugget=T)),
          gls(integr ~ maxent, df, correlation = corGaus(form = ~lon + lat, nugget=T)),
          gls(integr ~ maxent, df, correlation = corSpher(form = ~lon + lat, nugget=T)),
#          gls(integr ~ domain, df, correlation = corLin(form = ~lon + lat, nugget=T)),
          gls(integr ~ maxent, df, correlation = corRatio(form = ~lon + lat, nugget=T))
          )

# exponential correlation function is the best model!
# Check variogram of each model
plot(nlme:::Variogram(gls(extent ~ domain, df), form = ~lat + lon, resType = "normalized"))

## 5.3. Check spatial autocorrelation in a single variable  ####
library(pgirmess)

# Import dataframe with 3 columns (lat, lon and variable)
xy <- na.omit(xy)
coords <- xy[, c("lat","lon")] # lat and lon in columns 1 and 2

# Moran's I results
corr_log <- correlog(coords, xy$extent, nbclass = 7)
plot(corr_log) #correlogram


