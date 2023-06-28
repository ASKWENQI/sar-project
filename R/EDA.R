rm(list=ls())
gc()

# loading and filtering
library(phyloseq)
neon_dob <- readRDS("phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
neon <- subset_samples(neon, !is.na(horizon))
#rm(neon_dob)

# neon dataset
d <- sample_data(neon) # sample data data frame
d$collectYear <- as.factor(d$collectYear)
aggregate(d$collectYear,list(d$siteID),function(x) length(unique(x)))
a <- unique(d$Site) # get all the sites

d_site <- d[which(d$siteID==a[2]),]
boxplot(d_site$soilInCaClpH ~ d_site$collectYear)
d_site <- d[which(d$siteID==a[3]),]
boxplot(d_site$soilInCaClpH ~ d_site$collectYear)
d_site <- d[which(d$siteID==a[4]),]
boxplot(d_site$soilInCaClpH ~ d_site$collectYear)
d_site <- d[which(d$siteID==a[42]),]
boxplot(d_site$soilInCaClpH ~ d_site$collectYear)
# delete 2014 data
# MLBS no 2017 pH data

d_site <- d[which(d$siteID==a[1]),]
boxplot(d_site$soilMoisture ~ d_site$collectYear)
d_site <- d[which(d$siteID==a[2]),]
boxplot(d_site$soilMoisture ~ d_site$collectYear)
d_site <- d[which(d$siteID==a[4]),]
boxplot(d_site$soilMoisture ~ d_site$collectYear)
d_site <- d[which(d$siteID==a[39]),]
boxplot(d_site$soilMoisture ~ d_site$collectYear)


d_site <- d[which(d$siteID==a[2]),]
boxplot(d_site$soilInCaClpH ~ d_site$sampleTiming)
d_site <- d[which(d$siteID==a[3]),]
boxplot(d_site$soilInCaClpH ~ d_site$sampleTiming)
d_site <- d[which(d$siteID==a[4]),]
boxplot(d_site$soilInCaClpH ~ d_site$sampleTiming)
d_site <- d[which(d$siteID==a[40]),]
boxplot(d_site$soilInCaClpH ~ d_site$sampleTiming)
# delete 2014 data
summary(aov(d_site$soilInCaClpH ~ d_site$collectYear))

d_site <- d[which(d$siteID==a[1]),]
boxplot(d_site$soilMoisture ~ d_site$sampleTiming)
d_site <- d[which(d$siteID==a[2]),]
boxplot(d_site$soilMoisture ~ d_site$sampleTiming)
d_site <- d[which(d$siteID==a[4]),]
boxplot(d_site$soilMoisture ~ d_site$sampleTiming)
d_site <- d[which(d$siteID==a[42]),]
boxplot(d_site$soilMoisture ~ d_site$sampleTiming)

p_value <- vector(length = 45)
for (i in c(1:41,43:45)){
  d_site <- d[which(d$siteID==a[i]),]
  
  if (length(unique(d_site$collectYear)) > 1){
    p_value[i] <- summary(aov(d_site$soilInCaClpH ~ d_site$collectYear))[[1]][["Pr(>F)"]][1]
  }else{
    p_value[i] <- NA
  }
}
p_value[42] <- NA
sum(p_value < 0.05, na.rm = T)

p_value <- vector(length = 45)
for (i in c(1:45)){
  d_site <- d[which(d$siteID==a[i]),]
  
  if (length(unique(d_site$sampleTiming)) > 1){
    p_value[i] <- summary(aov(d_site$soilInCaClpH ~ d_site$sampleTiming))[[1]][["Pr(>F)"]][1]
  }else{
    p_value[i] <- NA
  }
}
p_value
sum(p_value < 0.05, na.rm = T)

# soil moisture
p_value <- vector(length = 45)
for (i in c(1:45)){
  d_site <- d[which(d$siteID==a[i]),]
  
  if (length(unique(d_site$collectYear)) > 1){
    p_value[i] <- summary(aov(d_site$soilMoisture ~ d_site$collectYear))[[1]][["Pr(>F)"]][1]
  }else{
    p_value[i] <- NA
  }
}
p_value
sum(p_value < 0.05, na.rm = T)

p_value <- vector(length = 45)
for (i in c(1:45)){
  d_site <- d[which(d$siteID==a[i]),]
  
  if (length(unique(d_site$sampleTiming)) > 1){
    p_value[i] <- summary(aov(d_site$soilMoisture ~ d_site$sampleTiming))[[1]][["Pr(>F)"]][1]
  }else{
    p_value[i] <- NA
  }
}
p_value
sum(p_value < 0.05, na.rm = T)

# sampleTiming * collectYear
p_value <- vector(length = 45)
for (i in c(1:45)){
  d_site <- d[which(d$siteID==a[i]),]
  
  if (length(unique(d_site$sampleTiming)) > 1){
    p_value[i] <- summary(aov(d_site$soilMoisture ~ d_site$sampleTiming * d_site$collectYear))[[1]][["Pr(>F)"]][1]
  }else{
    p_value[i] <- NA
  }
}
p_value
sum(p_value < 0.05, na.rm = T)

p_value <- vector(length = 45)
for (i in c(1:45)){
  d_site <- d[which(d$siteID==a[i]),]
  b <- paste(d_site$collectYear,d_site$sampleTiming, sep = "_")
  boxplot(d_site$soilMoisture ~ b)
  if (length(unique(b)) > 1){
    p_value[i] <- summary(aov(d_site$soilMoisture ~ b))[[1]][["Pr(>F)"]][1]
  }else{
    p_value[i] <- NA
  }
}
p_value
sum(p_value < 0.05, na.rm = T)

p_value <- vector(length = 45)
for (i in c(1:45)){
  d_site <- d[which(d$siteID==a[i]),]
  b <- paste(d_site$collectYear,d_site$sampleTiming, sep = "_")
  boxplot(d_site$soilInCaClpH ~ b)
  if (length(unique(b)) > 1){
    p_value[i] <- summary(aov(d_site$soilInCaClpH ~ b))[[1]][["Pr(>F)"]][1]
  }else{
    p_value[i] <- NA
  }
}
p_value
sum(p_value < 0.05, na.rm = T)
