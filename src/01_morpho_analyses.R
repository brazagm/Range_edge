## ---------------------------------------------------------------------------##
##                                                                            ##
##  Title: MORPHOLOGICAL ANALYSES                                             ##
##                                                                            ##
##  Description: Morphological data preparation and analyses. Import          ##
##  occurrence records, populations and trait measures of the Gray slender    ##
##  opossum, Marmosops incanus. Export measurement errors for each trait      ##
##  test sexual dimorphism and calculate trait measures corrected by size.    ##
##                                                                            ##
##  Created by: Alan Braz (@brazagm)                                          ##
##                                                                            ##
## ---------------------------------------------------------------------------##
setwd("./")

library(ICC)
library(smatr)


# 0. DATA INPUT AND OUTPUT  ------------------------------------------------####
## Import data (i.e.occurrence, populations and traits of Marmosops incanus)
data <- read.csv("./data/data_morphology.csv")
str(data)

## Set the name of morphological traits
trait_names <- c("gls", "cb", "nas", "bbc", "zb", "poc", "ic", "lbasi", "bps", "bam", "can", "pl", "max", "sm", "mad", "alpcor", "larcon") # name of linear measures

## Set the name of functional traits (consult Material & Methods for further details on the calculations of functional traits)
func_names <- c("size", "shape", "nas.gls", "can.bam", "pl.cb", "bam.pl", "sm.max", "alpcor.mad")


# 1. MEASUREMENT ERRORS  ---------------------------------------------------####
## Subset data by selecting only specimens which skull has repeated measurements
df <- subset(data, sp == "Marmosops incanus" & idade == "adulto" & replicas == TRUE)

## Calculate sufficient sample size (number of specimens) to estimate measurement errors
# "h" = a hypothetical test, no preliminary results
# w = confidence interval, this value has to be smaller as possible
# ICC = ICC value which I wish to estimate, +/- the confidence interval
# k = repetition number to be done
# alpha = p-value for test significance
(n <- Nest("h", w = 0.1, ICC = 0.9, k = 5, alpha = 0.05))

## 1.1. Calculate the percentage of measurement error (%ME)  ####
## Create a matrix for data input
tab <- matrix(data = NA, nrow = 3, ncol = length(trait_names))
rownames(tab) <- c("p", "%", "ICC")
colnames(tab) <- trait_names

## Automate the calculation for each trait
for (x in trait_names) {

  # Import replicates ('r') of the measurement ('x') for each specimen
  m <- cbind(
    Trait = c(
      df[, x],
      df[, paste(x, ".1", sep = "")],
      df[, paste(x, ".2", sep = "")],
      df[, paste(x, ".3", sep = "")],
      df[, paste(x, ".4", sep = "")]
    ),
    Ind = rep(1:length(df[, x]), 5)
  )
  m <- as.data.frame(m)

  # Test if the individual as a factor in ANOVA
  mod <- summary(aov(Trait ~ as.factor(Ind), data = m))
  mod

  ## Calculate the percentage of measurement error (%ME)
  # The percentage of measurement error is defined as the ratio of the within
  # measurement component of variance on the sum of the within- and among-measurement component
  s2within <- MSwithin <- mod[[1]][2, 3] # variance within the same specimen is equal to the ANOVA's mean squares
  MSamong <- mod[[1]][1, 3] # mean squares among specimens
  s2among <- (MSamong - MSwithin) / 4 # difference among the mean squares divided by the number of replicates

  tab["%", x] <- s2within / (s2within + s2among) * 100 # calculate the percentage of mesaurement error
  tab["p", x] <- mod[[1]][1, 5]
  tab["ICC", x] <- s2among / (s2within + s2among) # calculate ICC
}

## Export results
write.csv(t(tab), "./results/01_measur_error.csv")


# 2. SEXUAL DIMORPHISM  ----------------------------------------------------####
## Select only adult specimens with sex identification and skull measures
df <- subset(data, sp == "Marmosops incanus" & sexo %in% c("male", "female") & idade == "adulto" & medido == TRUE & pop != "")

## Read 'sex' column as factor
df$sexo <- factor(df$sexo)

## Summary
table(df$sexo)
summary(df)

## 2.1. COMPARE SEXES WITH BOXPLOTS  ####
# Create boxplots for each trait between sexes, comparing the mean value, sd, and ci 95 and 99%
par(mfrow = c(4, 10), mar = c(2, 2, 2, 2))

for (i in trait_names) {
  x <- df[, c(i, "sexo")]
  x <- x[complete.cases(x), ] # remove NAs

  ## Create a boxplot based on the median, quartis and range
  ju <- boxplot(x[, i] ~ x$sexo, range = 0, notch = TRUE, main = i, border = c("blue", "red"))

  ## Create a boxplot based on the mean, standard deviation and 95% en 99% confidence intervals
  # For that purpose, the values of some components of the boxplot 'ju' are modified
  mean <- tapply(x[, i], x$sexo, mean) # mean for each group
  sd <- tapply(x[, i], x$sexo, sd) # standard deviation
  se <- sd / sqrt(as.numeric(table(x$sexo))) # standard error

  # Graphical output
  ju$stats[3, ] <- mean
  ju$stats[c(2, 4), ] <- rbind(mean - qnorm(0.995) * se, mean + qnorm(0.995) * se)
  ju$stats[c(1, 5), ] <- rbind(mean - sd, mean + sd)
  ju$conf <- rbind(mean - qnorm(0.975) * se, mean + qnorm(0.975) * se)
  bxp(ju, notch = TRUE, main = i, border = c("blue", "red"))
}

## 2.2. TEST SEXUAL DIMORPHISM USING ANOVA  ####
## Use ANOVA for test the differences between mean and variance of male/females
# Create a matrix for the result
matrix <- matrix(NA, 1, length(trait_names))
colnames(matrix) <- trait_names
rownames(matrix) <- "Pr(>F)"

# Are traits mean and variance different among sexes?
for (i in trait_names) {

  # ANOVA and F-test
  matrix[1, i] <- summary(aov(df[, i] ~ df$sexo))[[1]][1, 5]
}

# p-values for the ANOVA for each trait
matrix

# Export the results
write.csv(matrix, "./results/02_sex_dimorph_anova.csv")


## 2.3. TEST SEXUAL DIMORPHISM ON SKULL SHAPE ####
## Is dimorphism only a size difference? Or sexes also differ in skull shape?
# Here, size is defined as the geometric mean of CB and ZB. However, size could be the geometric mean of all the skull traits... Size was already calculated by a function in LibreOffice calc
# Sexual dimorphism on skull shape was tested using traits mean of each sex per population
## Comparisons of the allometry among sexes were performed using the 'sma' function from 'smatr' package
# Create a matrix for the result
matrix <- data.frame(trait = trait_names)

df$pop <- as.factor(df$pop)
df <- droplevels(df)
pop <- data.frame(sex_pop = paste("female", levels(df$pop), sep = "_"))
pop <- rbind(pop, data.frame(sex_pop = paste("male", levels(df$pop), sep = "_")))

# Extract the variability of each population, for each type of measures (linear and functional trait)
for (i in levels(df$pop)) {
  for (t in trait_names) {
    y_female <- subset(df, pop == i & sexo == "female")
    y_male <- subset(df, pop == i & sexo == "male")

    pop[which(pop$sex_pop == paste("female", i, sep = "_")), t] <- mean(y_female[, t])
    pop[which(pop$sex_pop == paste("female", i, sep = "_")), "sexo"] <- "female"
    pop[which(pop$sex_pop == paste("male", i, sep = "_")), t] <- mean(y_male[, t])
    pop[which(pop$sex_pop == paste("male", i, sep = "_")), "sexo"] <- "male"
  }
  pop[which(pop$sex_pop == paste("female", i, sep = "_")), "size"] <- mean(y_female[, "size"])
  pop[which(pop$sex_pop == paste("male", i, sep = "_")), "size"] <- mean(y_male[, "size"])
}

# Test shift and slope differences in allometry of male and females
for (i in trait_names) {
  slope <- sma(pop[, i] ~ pop$size * as.factor(pop$sexo)) # test differences in slopes
  shift <- sma(pop[, i] ~ pop$size + as.factor(pop$sexo), type = "shift") # test for shifts

  # Compute data in result matrix
  matrix[which(matrix$trait == i), "male-slope"] <- slope$coef$male[2, 1]
  matrix[which(matrix$trait == i), "male-r²"] <- slope$r2$male[1]
  matrix[which(matrix$trait == i), "male-p"] <- slope$pval$male[1]
  matrix[which(matrix$trait == i), "female-slope"] <- slope$coef$female[2, 1]
  matrix[which(matrix$trait == i), "female-r²"] <- slope$r2$female[1]
  matrix[which(matrix$trait == i), "female-p"] <- slope$pval$female[1]
  matrix[which(matrix$trait == i), "slope-test"] <- slope$commoncoef$p[1]
  matrix[which(matrix$trait == i), "shift-test"] <- shift$gtr$p[1, 1]
}

# Export results
write.csv(matrix, file = "./results/03_sex_dimorph_shape.csv")


## 2.4. FIGURE: SEXUAL DIMORPHISM ####
source("./src/function_corner_text.R")

## Plot examples of size-only dimorphism and shape of the skull
par(mfrow = c(1, 2), oma = c(2, 2, 2, 0))

plot(sma(df$bps ~ df$size * df$sexo),
  xlab = "Skull size", ylab = "NAS"
)
Corner_text("A", location = "bottomright")

plot(sma(df$sm ~ df$size * df$sexo),
  xlab = "Skull size", ylab = "BPS"
)
Corner_text("B", location = "bottomright")


# 3. SIZE CORRECTION  ------------------------------------------------------####
## Subset data - only measured and adult specimens of M. incanus, with sex identified
df <- subset(data, sp == "Marmosops incanus" & idade == "adulto")

# Extract regression residuals for each trait
for (i in trait_names) {
  m <- sma(df[, i] ~ df$size)
  df[which(is.na(df[, i]) == FALSE & is.na(df$size) == FALSE), i] <- residuals(m)
  df[which(is.na(df$size) == TRUE), i] <- NA
}

# Export size-corrected data
write.csv(df, "./data/data_size-corrected.csv")


# 4. EXPLORATORY ANALYSES  -------------------------------------------------####
## Subset data: size-corrected measures of M. incanus' adult specimens...
df <- read.csv("./data/data_size-corrected.csv")

## 4.1. Are traits normally distributed?  ####
# Linear traits
par(mfrow = c(5, 4), mar = c(2, 2, 1, 1), oma = c(0.5, 0.5, 0.5, 0.5))
for (i in trait_names) {
  t <- na.omit(df[, i])

  hist(t,
    breaks = 10, prob = TRUE, col = "grey",
    main = paste(i), xlab = "", ylab = "", axes = TRUE
  )
  lines(density(t), col = "red", lwd = 2)
}

# Functional traits
par(mfrow = c(2, 4), mar = c(1, 1, 1, 1), oma = c(0.5, 0.5, 0.5, 0.5))
for (i in func_names) {
  t <- na.omit(df[, i])

  hist(t,
    breaks = 10, prob = TRUE, col = "grey", xaxt = "n", yaxt = "n",
    main = paste(i), xlab = "", ylab = "", axes = TRUE
  )
  lines(density(t), col = "red", lwd = 2)
}


## 4.2. Are outliers present? Check QQ-plot  ####
# Linear traits
par(mfrow = c(5, 4), mar = c(2, 2, 1, 1), oma = c(0.5, 0.5, 0.5, 0.5))
for (i in trait_names) {
  t <- na.omit(df[, i])

  qqnorm(t, main = paste(i))
  qqline(t)
}

# Functional traits
par(mfrow = c(2, 4), mar = c(2, 2, 1, 1), oma = c(0.5, 0.5, 0.5, 0.5))
for (i in func_names) {
  t <- na.omit(df[, i])

  qqnorm(t, main = paste(i))
  qqline(t)
}

## 4.3. Shapiro-Wilk's test of normality  ####
for (i in trait_names) {
  t <- na.omit(df[, i])


  if (shapiro.test(t)$p.value < 0.05) {
    print(paste(i, "'s distribution differs from normal distribution (P < 0.05)", sep = ""))
  }
}

