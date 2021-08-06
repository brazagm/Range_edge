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

##++++++++++++++++++++++++++++++##
## 1. FUNCTIONAL VARIABILITY  ####
##++++++++++++++++++++++++++++++##
## Import data
# Obs: functional traits are calculated as the linear measures without correction
# because the values were extracted from the columns with values calculated before correction (i.e. nas.gls, can.bam,...)
data <- read.csv("./Data_Size-corrected.csv") # size-corrected morphological data
str(data)

## Set the name of morphological traits
trait.names <- c("gls", "cb", "nas", "bbc", "zb", "poc", "ic", "lbasi", "bps", "bam", "can", "pl", "max", "sm", "mad", "alpcor", "larcon") # name of linear measures

# Remove measures which measurement error > 20% and have sexual dimorphism in shape
trait.names <- trait.names[!trait.names %in% c("poc","bps","bam","larcon")]

## Set the name of functional traits
fun.names <- c("size", "shape", "nas.gls", "can.bam", "pl.cb", "bam.pl", "sm.max", "alpcor.mad") # name of functional traits - please, consult Material & Methods for further details on the calculations of functional traits

## Data filtering
# Only adult specimens of Marmosops incanus, with sex identified and populations defined
y <- subset(data, subset = sp == "Marmosops incanus" & idade == "adulto" & pop != "" & medido == TRUE)
y$pop <- factor(y$pop) # 'pop' as factor
populations <- levels(y$pop)  # population names


## 1.1. CALCULATE MORPHOLOGICAL VARIABILITY  ####
# Morphological variability was calculated based on the multidimensional approach described in Boucher et al. (2013)

## Standardize morphological measures
# Standardize the species' traits to mean = 0 and sd = 1, using the entire sample
y.std <- subset(y, select = c("pop", "lon", "lat", "sexo", trait.names, fun.names)) # Select only essential columns
y.std[, c(trait.names, fun.names)] <- psych::rescale(y.std[, c(trait.names, fun.names)], mean = 0, sd = 1) # standardize

## Calculate the variability extent (sensu Boucher et al. 2013)
# Create a dataframe for data input
# The variability extent is calculated based on two diferent measures
# 'raw' = linear measures of the skull corrected by size
# 'fun' = functional traits based on ratios of linear measures (e.g. can:bam, nas:gls)
df <- data.frame(population = populations, extent.raw = NA, extent.fun = NA)

## Extract the variability of each population, for each type of measures (linear and functional trait)
for(i in populations){
  
  ## 1) Calculate the variability extension of linear measures
  yi <- subset(y.std, subset = pop == i, select = c(trait.names))
  yi <- na.omit(yi) # omit NAs
  df[which(df$population == i), "n.raw"] <- nrow(yi) # sample size
  
  # var-covar matrix for linear measures
  varcor <- cov(yi)
  df[which(df$population == i), "extent.raw"] <- sum(diag(varcor)) # extract the extent value
  
  
  ## 2) Calculate the variability extension of functional traits
  yi <- subset(y.std, subset = pop == i, select = c(fun.names))
  yi <- na.omit(yi) # omit NAs
  df[which(df$population == i), "n.fun"] <- nrow(yi) # sample size
  
  # var-covar matrix for functional traits
  varcor <- cov(yi)
  df[which(df$population == i), "extent.fun"] <- sum(diag(varcor)) # extract the extent value
}


##+++++++++++++++++++++++++++##
## 2. PREDICTOR VARIABLES  ####
##+++++++++++++++++++++++++++##
## 2.1. IMPORT ENVIRONMENTAL SUITABILITY  ####
## Import niche models
# Ecological niche models were calibrated according to Braz et al. (2020) "Interspecific competition constrains local abundance in highly suitable areas". Ecography 43: 1560–1570
list <- list.files("./Model_results", pattern = "avg.asc", full.names = TRUE)
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


## 2.2. EXTRACT SUITABILITY VALUES ####
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


## 2.3. DISTANCE FROM THE RANGE EDGE ####
## Import occurrence records
# Records are filtered by a minimum distance of 10 km between them to avoid spatial bias in sample
points <- read.csv("./Data_Records_thin1.csv") # import csv
points <- SpatialPoints(points[, c("lon","lat")], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) # as spatialpoints

# Create a Minimum Convex Polygon (MCP)
limits <- readOGR("/home/alan/Documentos/Repository/GIS/shapes/borders/", layer = "TM_WORLD_BORDERS-0.3")
mcp <- mcp(points, percent = 100)
intersect <- intersect(limits, mcp)
df.points <- SpatialPointsDataFrame(df[,c("lon","lat")], df)
dist.line <- dist2Line(df.points, intersect) # calculate distance from points to MCP
df.points@data$distance <- dist.line[,"distance"]
df <- df.points@data


##++++++++++++++++++++++++++++##
## 3. STATISTICAL ANALYSES  ####
##++++++++++++++++++++++++++++##
## 3.1. TEST THE VARIABILITY-SUITABILITY RELATIONSHIP  ####
# Select predictors
pred.names <- c("domain","mahal","maxent", "distance")

# Data filtering
xy <- subset(df, subset = n.raw >= 5) # Only populations with sample size >= 5 individuals

# Something gone wrong: Catas Altas and Reserva Vale have NA values for suitabilities. Add values manually
xy[which(xy$population == "Catas Altas"), c("domain","mahal","maxent")] <- c(0.8806858, 0.9085349, 0.9675414)
xy[which(xy$population == "Reserva Vale"), c("domain","mahal","maxent")] <- c(0.9235663, 0.9540392, 0.5223794)

# Remove outliers (source from https://www.r-bloggers.com/identify-describe-plot-and-remove-the-outliers-from-the-dataset/)
# Explaratory analyses on residual distribution
source("http://goo.gl/UUyEzD")
outlierKD(xy, extent.raw) # YES: Remove 1 outlier from 26 observations


## 3.1.1 Test for mean differences and variance equality in variability between low/high populations ####
# Classify each locality as low/high suitability for each ENM
for(i in pred.names){
  
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
for(i in pred.names){
  library(car) # required for Levene's test
  
  # Characterize each sample as low/high environmental suitability
  y.low <- xy[which(xy[, paste(i, "binary", sep = "_")] == "low"),"extent.raw"]
  y.high <- xy[which(xy[, paste(i, "binary", sep = "_")] == "high"),"extent.raw"]
  x <- xy[, paste(i, "binary", sep = "_")]
  x <- as.factor(x)
  y <- xy$extent.raw
  
  # Calculate Levene's test, Welch's test (t-test with unequal variances), and ANOVA
  lev <- leveneTest(y, group = x, center= mean)
  welch <- t.test(y.low, y.high, alternative = "two.sided", var.equal = FALSE)
  anova <- aov(y ~ x); anova <- summary(anova)
  
  # Compute the results
  results[c("levene.F","levene.p"),i] <- c(lev$`F value`[1], lev$`Pr(>F)`[1])
  results[c("welch.t","welch.p"),i] <- c(welch$statistic, welch$p.value)
  results[c("anova.F","anova.p"),i] <- c(anova[[1]][["F value"]][1], anova[[1]][["Pr(>F)"]][1])
}

# Export results
write.csv(results, "./Results/Results_median groups.csv")


## 3.1.2. Test for linear relationship between variability and predictor variables ####
# Test if morphological variability data follows a gamma distribution using Kolmogorov-Smirnov test
gamma <- fitdistr(na.omit(xy$extent.raw), "gamma") 
ks.test(na.omit(xy$extent.raw), "pgamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

## Fit linear regressions 
for(i in pred.names){
  
  if(i == "distance"){
    x <- log(xy[,i]) # log-transformation for distance from the range edge
  }else{
    x <- xy[,i]
  }
  
  # Model fitting
  lm <- lm(xy$extent.raw ~ x) # fit linear models
  glm <- glm(xy$extent.raw ~ x, family = Gamma(link = "log")) # fit GLM gamma
  
  if(i != "distance"){
    lm.distance <- lm(xy$extent.raw ~ x + log(xy$distance))
    assign(paste('lm.dist', i, sep = "."), lm.distance)
    glm.distance <- glm(xy$extent.raw ~ x + xy$distance, family = Gamma(link = "log"))
    assign(paste('glm.dist', i, sep = "."), glm.distance)
  }
  
  # Assign models to an object
  assign(paste('lm', i, sep = "."), lm)
  assign(paste('glm', i, sep = "."), glm)
}

# Model selection results
lm.sel <- MuMIn::model.sel(lm.domain, lm.dist.domain, lm.mahal, lm.dist.mahal, lm.maxent, lm.dist.maxent, lm.distance, lm(xy$extent.raw ~ 1))

glm.sel <- MuMIn::model.sel(glm.domain, glm.dist.domain, glm.mahal, glm.dist.mahal, glm.maxent, glm.dist.maxent, glm.distance, glm(xy$extent.raw ~ 1))

## Test if the dependency of lm residuals to the explanatory variable (i.e. heteroskedasticity) using Breusch-Pagan test
lmtest::bptest(lm.domain)
lmtest::bptest(lm.mahal)
lmtest::bptest(lm.maxent)
lmtest::bptest(lm.distance)

## Although a visual inspection of the xy plot suggests a certain degree of heteroskedasticity, variability does not seems to follow a gamma distribution and Breusch-Pagan test does not reject the homoskedasticity
## For this reason, results were based on simple linear regressions
write.csv(lm.sel, "./Results/Model_Selection_Results.csv")


##+++++++++++++++##
## 4. FIGURES  ####
##+++++++++++++++##
# Graphic parameters
# Create a function to automate the procedures
myPlot <- function(x, x.name, y.name){
  par(mfrow = c(1, 1),
  oma = c(0, 0, 0, 0), mar = c(4, 1, 1, 4), mgp = c(1.8, 0.5, 0),
  cex.lab = 1.4, family = "sans", font.lab = 2)
  
  plot(xy$extent_raw ~ x, pch = 1, cex = 1.2,
     yaxt = "n", ylab = NA, xlab = x.name)
  axis(4)
  mtext(y.name, 4, line = 2, cex = 1.4, family = "sans", font = 2)
}

# Figure 3A - Domain X Variability
pdf("Figure 3A.pdf",  width = 5, height = 4)

myPlot(x = xy$domain_std, x.name = "Environmental suitability", y.name = "Morphological variability")
abline(lm.domain_std, col = "black", lty = 1, lwd = 2)
legend("topleft", legend = "R² = 0.16
P = 0.04
n = 25", bty = "n", text.font = 1, adj = c(0,0), cex = 1, inset = c(0, 0.1))

points(subset(xy, distance < 20000, select = c("domain_std","extent_raw")), col = "black", pch = 19)

dev.off()

# Figure 3B - Mahalanobis X Variability
pdf("Figure 3B.pdf",  width = 5, height = 4)

myPlot(x = xy$mahal_std, x.name = "Environmental suitability", y.name = "Morphological variability")
points(subset(xy, distance < 20000, select = c("mahal_std","extent_raw")), col = "black", pch = 19)

dev.off()

# Figure 3C - Maxent X Variability
pdf("Figure 3C.pdf",  width = 5, height = 4)

myPlot(x = xy$maxent_std, x.name = "Environmental suitability", y.name = "Morphological variability")
points(subset(xy, distance < 20000, select = c("maxent_std","extent_raw")), col = "black", pch = 19)

dev.off()

# Figure 4 - Distance X Variability
pdf("Figure 4.pdf",  width = 5, height = 4)

myPlot(x = log(xy$distance), x.name = expression(bold("Distance log"[10]*"(m)")), y.name = "Morphological variability")
abline(lm.distance, col = "black", lty = 1, lwd = 2)
legend("topleft", legend = "R² = 0.23
P = 0.02
n = 25", bty = "n", text.font = 1, adj = c(0,0), cex = 1, inset = c(0, 0.1))

dev.off()



## FIGURE APPENDIX - CORRELATIONS ####
# Graphic parameters
# Create a function to automate the procedures
pdf("Figure S4.pdf",  width = 10, height = 4)

par(mfrow = c(1,3), oma = c(3, 3, 0, 0), mar = c(1, 1, 1, 1), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.4, family = "sans", font.lab = 2)

myPlot <- function(x, x.name, y.name){
  plot(x ~ log(xy$distance), pch = 1, cex = 1.2,
       ylab = "", xlab = x.name)
}

# Figure 3A - Domain X Variability
myPlot(x = xy$domain_std, x.name = "", y.name = "")
legend("bottomright", legend = "Domain", bty = "n", text.font = 2, adj = c(0,0), cex = 1.5, inset = c(0.05, 0))

myPlot(x = xy$mahal_std, x.name = "", y.name = "")
legend("bottomright", legend = "Mahalanobis", bty = "n", text.font = 2, adj = c(0,0), cex = 1.5, inset = c(0.05, 0))

myPlot(x = xy$maxent_std, x.name = "", y.name = "")
legend("bottomright", legend = "Maxent", bty = "n", text.font = 2, adj = c(0,0), cex = 1.5, inset = c(0.05, 0))

title(ylab = "Environmental suitability", cex.lab = 1.4, line = 1, outer = TRUE)
title(xlab = expression(bold("Distance log"[10]*"(m)")), cex.lab = 1.4, line = 1, outer = TRUE)


dev.off()







# Graphic results
for(x in pred.names){
  
  # Import models
  lm <- get(paste("lm", x, sep = "." ))
  s <- summary(lm)
  
  # Plot variability-suitability relationship
  if(x == "distance"){
    
    plot(xy$extent_fun ~ log(xy[,x]), pch = 19, cex = 1.2,
         ylab = "", xlab = "Distance from range border (log)")
    text(xy[,x], xy$extent_fun, labels = xy$population, cex= 0.7, adj = c(0,1))
  }else{
    
    plot(xy$extent_fun ~ xy[,x], pch = 19, cex = 1.2,
       ylab = "", xlab = "Environmental suitability")
    #text(xy[,x], xy$extent_fun, labels = xy$population, cex= 0.7, adj = c(0,1))
  }
  
  # Add model name
  if(x == "mahal"){
    legend("topleft", legend = "Mahalanobis", bty = "n", text.font = 2, adj = c(0.2,0), cex = 1.8, inset = c(0, 0))
  } else if(x == "domain"){
    legend("topleft", legend = "Domain", bty = "n", text.font = 2, adj = c(0.2,0), cex = 1.8, inset = c(0, 0))
  } else if(x == "maxent"){
    legend("topleft", legend = "Maxent", bty = "n", text.font = 2, adj = c(0.2,0), cex = 1.8, inset = c(0, 0))
  }
  
  # Plot significant lines (only tau = 90 and 95 for Domain)
  if(x == "domain"){
    abline(qr.90.domain, col = "black", lty = 2, lwd = 1.5)
  }

}
title(ylab = "Morphological variability", cex.lab = 1.4, line = 0.5, outer = TRUE)

dev.off()


# QUantile regression
model <- rq(xy$extent_fun ~ xy$maxent, 0.90)
summary(model, se = "boot")


## 3.2. Spatial autocorrelation ####
library(ape); library(MuMIn); library(nlme)

df <- data.frame(lon = xy$lon, lat = xy$lat)
df <- SpatialPointsDataFrame(xy[, c("lon","lat")], xy, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

df@data <- df@data[!df@data$pop == "Mocambinho",]

lm.names <- c("lm.domain_std", "lm.mahal_std", "lm.maxent_std", "lm.distance", "lm.dist.domain_std", "lm.dist.mahal_std", "lm.dist.maxent_std")

for(i in lm.names){
  model <- get(i)
  
  df@data[,paste(i, "resid", sep = "_")] <- model$residuals
}

sample.dist <- as.matrix(dist(df@data[,c("lon","lat")]))
sample.dist <- 1/sample.dist
diag(sample.dist) <- 0

table <- data.frame(model = lm.names, observed = NA, expected = NA, sd = NA, p.value = NA)

for(i in lm.names){
  data <- df@data[, paste(i, "resid", sep = "_")]
  
  moran <- Moran.I(data, sample.dist, alternative = "greater")
  
  table[which(table$model == i), "observed"] <- moran$observed
  table[which(table$model == i), "expected"] <- moran$expected
  table[which(table$model == i), "sd"] <- moran$sd
  table[which(table$model == i), "p.value"] <- moran$p.value
}

write.csv(table, "./Results/Results_autocorr.csv")




### 3.3. FIGURE - Analyses without outliers ####
# Select predictors
pred.names <- c(paste(c("domain","mahal","maxent"), "std", sep = "_"), "distance")

# Data filtering
xy <- subset(df, subset = n_raw >= 5) # Only populations with sample size >= 5 individuals

# Something gone wrong: Catas Altas and Reserva Vale have NA values for suitabilities. Add values manually
xy[which(xy$population == "Catas Altas"), c("domain_std","mahal_std","maxent_std")] <- c(0.8806858, 0.9085349, 0.9675414)
xy[which(xy$population == "Reserva Vale"), c("domain_std","mahal_std","maxent_std")] <- c(0.9235663, 0.9540392, 0.5223794)

# Fit linear and quantile regressions
for(i in pred.names){
  
  if(i == "distance"){
    x <- log(xy[,i])
  }else{
    x <- xy[,i]
  }
  
  # Model fitting
  lm <- lm(xy$extent_raw ~ x) # fit linear models
  glm <- glm(xy$extent_raw ~ x, family = Gamma(link = "log"))
  qr.75 <- rq(xy$extent_raw ~ x, 0.75) # fit quantile regression, tau = .75
  qr.90 <- rq(xy$extent_raw ~ x, 0.90) # quantile regression, tau = .90
  qr.95 <- rq(xy$extent_raw ~ x, 0.95) # quantile regression, tau = .95
  
  if(i != "distance"){
    lm.distance <- lm(xy$extent_raw ~ x + xy$distance)
    assign(paste('lm.dist', i, sep = "."), lm.distance)
  }
  
  # Assign models to an object
  assign(paste('lm', i, sep = "."), lm)
  assign(paste('glm', i, sep = "."), glm)
  assign(paste('qr', "75", i, sep = "."), qr.75)
  assign(paste('qr', "90", i, sep = "."), qr.90)
  assign(paste('qr', "95", i, sep = "."), qr.95)
}

## FIGURE
# Graphic parameters
pdf("Figure S2.pdf",  width = 6, height = 5)

par(mfrow = c(2, 2),
    oma = c(2, 2, 0, 0), mar = c(3, 1, 1, 1), mgp = c(1.8, 0.5, 0),
    cex.lab = 1.1, family = "sans", font.lab = 2)

# Create a function to automate the procedures
myPlot <- function(x, x.name, y.name){
  plot(xy$extent_raw ~ x, pch = 1, cex = 1.2,
       ylab = NA, xlab = x.name)
  #mtext(y.name, 4, line = 2, cex = 1.4, family = "sans", font = 2)
}

# Figure 3A - Domain X Variability
myPlot(x = xy$domain_std, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Domain", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = xy$mahal_std, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Mahalanobis", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = xy$maxent_std, x.name = "Environmental suitability", y.name = "")
legend("topright", legend = "Maxent", bty = "n", text.font = 2, adj = c(0.2,-0.3), cex = 1, inset = c(0, 0.1))

myPlot(x = log(xy$distance), x.name = expression(bold("Distance log"[10]*"(m)")), y.name = "Morphological variability")

title(ylab = "Morphological variability", cex.lab = 1.1, line = 0.5, outer = TRUE)

dev.off()













## 3.1. Categorizar populações em periféricas/centrais e simpátricas/alopátricas com Marmosops paulensis ####
# Identificar localidades periféricas
# Escolhidas com base no conhecimento dos autores sobre a distribuição da espécie
peri <- c("Belo Horizonte", "Boracéia", "Caraça", "Catas Altas", "Chapada Diamantina", "Cipó", "Felício", "Intervales", "Irapé", "Mocambinho", "Peão", "Peti", "Una")

# Identificar localidades simpátricas
# Com base em registros de M. paulensis na mesma localidade
paulensis <- c("Bocaína", "Caparaó", "Serra dos Órgãos", "Intervales")


# Categorizar populaçoes
for(i in df$population){
  if(i %in% peri){
    df[which(df$population == i), "location"] <- "peripheral"
  }else{
    df[which(df$population == i), "location"] <- "central"
  }
  
  if(i %in% paulensis){
    df[which(df$population == i), "paulensis"] <- "sympatric"
  }else{
    df[which(df$population == i), "paulensis"] <- "allopatric"
    }
}

med <- median(na.omit(df$distance))
df[which(df$distance > med), "location"] <- "peripheral"
df[which(df$distance <= med), "location"] <- "central"

## 3.2. Testar diferenças entre as populações periféricas/centrais ####
# Selecionar populações com um tamanho amostral mínimo (>= 5 indivíduos na população)
xy <- subset(df, subset = n_raw >= 5)

## 3.2.1. Medidas lineares ####
# Visualizar boxplot
boxplot(extent_fun ~ location, data = xy, ylab = "Extensão da variabilidade", xlab = ""); points(factor(xy$location), xy$extent_fun, col = "black", pch = 19)
boxplot(extent_raw ~ paulensis, data = xy); points(factor(xy$paulensis), xy$extent_fun, col = "red")

summary(aov(extent_raw~ location, data = xy))




##+++++++++++++++##
## 4. FIGURAS  ####
##+++++++++++++++##

border <- shapefile("/home/alan/Documentos/Repository/GIS/shapes/borders/estados.shp")

## 4.1. Mapas ilustrativos
pdf("Figure 1.pdf",  width = 5, height = 5)

plot(sdm[["maxent"]], ext = extent(-52, -35, -27, -9), 
     col = terrain.colors(100, rev = TRUE), legend.args = list(text = 'Environmental suitability', side = 4, font = 2, line = 2.5, cex = 0.8))
lines(border, add = T)
points(ppoints, pch = 3, cex = 0.6)
points(xy[,c("lon","lat")], pch = 17, cex = 1, col = "yellow")
points(xy[,c("lon","lat")], pch = 24, cex = 1, col = "black")

dev.off()



pdf("Figure 2.pdf",  width = 5, height = 5)

plot(kernel, ext = extent(-52, -35, -27, -9), 
     col = hcl.colors(100, palette = "viridis", rev = TRUE))
lines(border, add = T)

dev.off()





