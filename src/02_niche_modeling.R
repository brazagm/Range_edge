## ---------------------------------------------------------------------------##
##                                                                            ##
##  Title: ECOLOGICAL NICHE MODELING                                          ##
##                                                                            ##
##  Description: Ecological niche modeling procedures to estimate             ##
##  environmental suitability for Marmosops incanus. Same script version      ##
##  from Braz et al. (2020) Interspecific competition constrains local        ##
##  abundance in highly suitable areas. Ecography 43: 1560–1570               ##
##                                                                            ##
##  Created by: Alan Braz (@brazagm)                                          ##
##                                                                            ##
## ---------------------------------------------------------------------------##
setwd("./")

## Load required packages
library(adehabitatHR)
library(biomod2)
library(dismo)
library(raster)
library(rgdal)
library(rgeos)
library(spThin)
library(tidyverse)
library(usdm)


# 1. OCCURRENCE RECORDS AND PRESENCES  -------------------------------------####
## Import species' occurrence records
records <- read.csv("./data/data_records.csv")

## 1.1. Occurrence data cleaning  ####
# All records' location were visually checked in QuantumGIS
# Records located within urban area (i.e. coordinates were given by municipalities, valid_georef == FALSE) were excluded from analyses due to the NDVI variable used in model calibration
records <- subset(records, valid_georef == TRUE)

# Remove duplicated observations by coordinates
records <- records[which(duplicated(records[, c("lon", "lat")]) == FALSE), ]

# Count number of occurrence records by data source
lvl <- levels(as.factor(records$institution))
df <- data.frame(n = rep(NA, 7), row.names = c(lvl, "total")) # create a dataframe for data entry

# Include the number of occurrences in dataframe (df)
for(i in lvl){
  df[i, "n"] <- nrow(subset(records, institution == i))
}
df["total", "n"] <- sum(df[c(lvl), "n"]) # total number of records

# Check dataframe (df)
df

## 1.2. Thin species' presences by a minimum distance of 10 km  ####
thin(loc.data = records, lat.col = "lat", long.col = "lon",
     spec.col = "sp", thin.par = 10, reps = 50,
     locs.thinned.list.return = FALSE,
     write.files = TRUE, max.files = 1, out.dir = "./data",
     out.base = "data_records", write.log.file = FALSE, verbose = TRUE)


# 2. PREDICTOR VARIABLES  --------------------------------------------------####
## Import thinned records
records <- read.csv("./data/data_records_thin1.csv")

## 2.1. Import variables  ####
# List variables from WorldClim 2.0
list <- list.files("./layers/WorldClim25_Global_v2", pattern = ".tif", full.names = TRUE)

# Add variables from MOD13A3 to the list
list <- append(list, list.files("./layers/MOD13A3", pattern = ".tif", full.names = TRUE))

# Set an equal extent for each variables dataset
e <- extent(raster(list[20])) # use MOD13's variables to reduce the extent

# Import as raster each variables by croping them to the same extent
for(i in list){
  if(i == list[1]){
    r <- stack(crop(raster(i), e))
  } else {
    r <- addLayer(r, crop(raster(i), e))
  }
}

# Drop layers with no or few ecological meaning for distribution modelling
r <- dropLayer(r, c("bio01","bio02","bio03","bio08","bio13","bio16","MOD13_NDVImin","MOD13_NDVImax","MOD13A3_IEVI","MOD13_EVImax","MOD13_EVImin"))
r_pred <- r

## 2.2. Set calibration area  ####
# Set CRS
proj4string(r) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Delimit calibration area as the Minimum Convex Polygon (MCP) with a 100 km buffer
# Create a MCP including all records
msk <- adehabitatHR::mcp(SpatialPoints(records[,c("lon","lat")], CRS("+proj=longlat +datum=WGS84 +no_defs")), percent = 100)
msk <- spTransform(msk, CRS("+init=epsg:32724")) # Project MCP shapefile
msk <- raster::buffer(msk, width = 100*1000, dissolve = TRUE) # Create a 200 km buffer
msk <- spTransform(msk, CRS("+proj=longlat +datum=WGS84 +no_defs")) # Back to WGS84 unprojected

# Crop predictor variables by the mask and a new extent
r <- crop(mask(r, msk), extent(-60, -30, -40, 0))
r_pred <- crop(mask(r_pred, msk), extent(-60, -30, -40, 0))

## 2.3. Variables selection ####
# Create 10,000 random points to sample climatic conditions
bg <- randomPoints(r, n = 10000, ext = extent(r[[1]]), extf = 1) 
colnames(bg) <- c('lon', 'lat') # Name columns
vals <- extract(r, bg) # Extract the climatic values

# Variation Inflation Factor (VIF)
vifstep(vals, th = 5) # Stepwise exclusion
vifcor(vals, th = 0.7) # Exclusion by correlation

# Variables selected by vifstep(th = 5)
select <- c("bio07","bio10","bio12","bio15","bio18","bio19", "MOD13A3_INDVI") # include "MOD13A3_INDVI"
r <- subset(r, select) # Rasters for calibration
r_pred <- subset(r_pred, select) # Rasters for model prediction


# 3. ECOLOGICAL NICHE MODELING  --------------------------------------------####
## Read presence records as a spatial points dataframe
pres <- SpatialPointsDataFrame(records[, c("lon","lat")], records, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

## Set names and attribute a numeric value for presences (= 1)
names(pres) <- c("M.incanus", "lon", "lat") ; pres$M.incanus <- 1

## 3.1. Create pseudo-absences & background points  ####
# Pseudo-absences were randomly generated with a 100-km minimum distance of presence records
# Create a mask with a buffer from presences
msk_pseudo <- raster::buffer(pres, width = 100*1000, dissolve = TRUE)

# Create pseudo-absences outside the mask
pseudo <- randomPoints(mask(r, msk_pseudo, inverse = TRUE), n = 1000, ext = extent(r), extf = 1)
colnames(pseudo) <- c("lon", "lat") # Set column names
pseudo <- as.data.frame((pseudo))

# Bind presences and pseudo-absences in the same dataframe
#pres_pseudo <- rbind(records[,c("lon","lat")], pseudo)
#pres_pseudo$Mincanus <- c(rep(1, nrow(records)), rep(0, nrow(pseudo)))

# Import presences/pseudo-absences dataframe as spatial points data frame
pseudo$M.incanus <- 0

# Read absences as a spatial points dataframe
pseudo <- SpatialPointsDataFrame(pseudo[, c("lon","lat")], pseudo, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Creat background points within the calibration area
background <- randomPoints(r, n = 10000, ext = extent(r), extf = 1)

## 3.2. Model calibration ####
# Presence-only and presence-background methods (i.e. Mahalanobis, Domain, and Maxent) were calibrated in dismo package
# Presence-absence (i.e. GLM and BRT) were calibrated in biomod2 package

## Automate modelling procedures, by calibrating replicates and combining them into an average model for each algorithm through:
# (3.2.1) Data preparation; (3.2.2) Model calibration; (3.2.3) Model prediction; and (3.2.4) Model evaluation

## Group each record type by kfold = 5 repetitions
pres_group <- kfold(pres, 5)
pseudo_group <- kfold(pseudo, 5)
back_group <- kfold(background, 5)
  
for(n in c(1:5)){
  ### 3.2.1 Data preparation  ####
  # For presence records
  pres_train <- pres[pres_group != n, ]
  pres_test <- pres[pres_group == n, ]
  
  # For pseudo-absences
  pseudo_train <- pseudo[pseudo_group != n, ]
  pseudo_test <- pseudo[pseudo_group == n, ]
  
  # For background points
  back_train <- background[back_group != n, ]
  back_test <- background[back_group == n, ]
  
  ## Formating data for biomod2
  resp <- c(pres_train$M.incanus, pseudo_train$M.incanus)
  xy <- rbind(pres_train[,c("lon","lat")], pseudo_train[,c("lon","lat")])
  xy <- xy@data
  data <- BIOMOD_FormatingData(resp.var = resp,
                               expl.var = r,
                               resp.xy = xy,
                               resp.name = "M.incanus",
                               na.rm = TRUE)
  options <- BIOMOD_ModelingOptions()
  
  ### 3.2.2. Model calibration  ####
  # in dismo package
  mahal <- mahal(r, pres_train)
  domain <- domain(r, pres_train)
  maxent <- maxent(r, pres_train, a = back_train, removeDuplicates = TRUE, args = "outputformat=logistic")
  
  # in biomod package
  out <- BIOMOD_Modeling(data,
                         models = c("GLM","GBM"),
                         models.options = options,
                         SaveObj = TRUE,
                         rescal.all.models = TRUE,
                         do.full.models = FALSE,
                         modeling.id = "")

  ### 3.2.3. Model prediction  ####
  # in dismo package
  pred.mahal <- predict(r_pred, mahal, progress = "")
  pred.domain <- predict(r_pred, domain, progress = "")
  pred.maxent <- predict(r_pred, maxent, progress = "")
  
  # in biomod package
  pred.biomod <- BIOMOD_Projection(modeling.output = out,
                            new.env = r_pred,
                            proj.name = "current",
                            selected.models = "all",
                            build.clamping.mask = FALSE,
                            compress = "xz",
                            output.format = ".grd",
                            do.stack = TRUE)
  pred <- stack("M.incanus/proj_current/proj_current_M.incanus.grd")
  
  # Stack raster predictions
  pred <- stack(pred.mahal, pred.domain, pred.maxent, pred)
  names(pred) <- paste(c("mahal", "domain", "maxent", "glm", "brt"), n, sep = "") # set names
  
  ### 3.2.4. Model evaluation  ####
  library(pROC)
  
  if(n == 1){
    df <- data.frame(model = names(pred), pAUC = NA, sens = NA, spec = NA, TSS = NA)
  }else{
    df <- rbind(df, data.frame(model = names(pred), pAUC = NA, sens = NA, spec = NA, TSS = NA))
  }
  
  p <- extract(pred, pres_test)
  a <- extract(pred, back_test)
  p_train <- extract(pred, pres_train)
  a_train <- extract(pred, back_train)
  
  for(i in names(pred)){
    # Calculate partial AUC
    table <- rbind(cbind(pred = p[,i], obs = 1), cbind(pred = a[,i], obs = 0))
    curve <- roc(as.factor(table[,"obs"]), table[,"pred"])
    pauc <- auc(curve, partial.auc = c(min(curve$sensitivities[curve$sensitivities > 0]), 1), partial.auc.focus = "sens", partial.auc.correct = FALSE)
    df[which(df$model == i), "pAUC"] <- pauc 
    
    # Get threshold value that maximizes TSS score
    values <- c(p_train[,i], a_train[,i])
    table <- data.frame(threshold = values, sens = NA, spec = NA, TSS = NA)
    
    for(x in unique(values)){
      binary <- pred[[i]] > x
      
      pres.values <- extract(binary, pres_test)
      pres.values <- na.omit(pres.values)
      true.pres <- length(which(pres.values == 1))
      false.pres <- length(which(pres.values == 0))
      
      abse.values <- extract(binary, pseudo_test)
      abse.values <- na.omit(abse.values)
      true.abse <- length(which(abse.values == 0))
      false.abse <- length(which(abse.values == 1))
      
      sensibility <- true.pres/(true.pres + false.pres)
      specificity <- true.abse/(true.abse + false.abse)
      tss <- (sensibility + specificity) - 1
      
      table[which(table$threshold == x), "sens"] <- sensibility
      table[which(table$threshold == x), "spec"] <- specificity
      table[which(table$threshold == x), "TSS"] <- tss
    }
    
    table <- na.omit(table)
    
    df[which(df$model == i), c("sens", "spec", "TSS")] <- unique(table[which(table$TSS == max(table$TSS)), c("sens", "spec", "TSS")])
    
  }
  
  # Export asc file
  writeRaster(pred, filename = "./results/niche_models/", format = "ascii", bylayer = TRUE, suffix = names(pred), overwrite = TRUE)
  
  
}

write.csv(df, "./results/niche_models/01_replicates_eval.csv")

## 3.3. Average model ####
list <- list.files("./results/niche_models", pattern = ".asc", full.names = TRUE)
model.names <- c("mahal","domain","maxent","glm","brt")

df.avg <- data.frame(model = model.names, n = NA, pAUC_mean = NA, pAUC_sd = NA, TSS_mean = NA, TSS_sd = NA)

for(i in model.names){
  sub.list <- grep(i, list, value = TRUE)
  model <- stack(sub.list)
  names(model) <- paste(i, c(1:5), sep = "")
  
  good <- df[which(df$TSS >= 0.5), "model"]
  good <- grep(i, good, value = TRUE)
  
  if(length(good) >= 1){
    model <- subset(model, good, drop = FALSE)
    
    avg <- mean(model)
    
    df.avg[which(df.avg$model == i), "n"] <- length(good)
    df.avg[which(df.avg$model == i), "pAUC_mean"] <- mean(df[which(df$model %in% good), "pAUC"])
    df.avg[which(df.avg$model == i), "pAUC_sd"] <- sd(df[which(df$model %in% good), "pAUC"])
    df.avg[which(df.avg$model == i), "TSS_mean"] <- mean(df[which(df$model %in% good), "TSS"])
    df.avg[which(df.avg$model == i), "TSS_sd"] <- sd(df[which(df$model %in% good), "TSS"])
    
    writeRaster(avg, paste("./Model_results/", i, "_avg.asc", sep = ""), overwrite = TRUE)
  }
  
}

write.csv(df.avg, "./results/niche_models/02_avg_eval.csv")

