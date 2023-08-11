rm(list=ls())
library(phyloseq)
library(dplyr)
neon_dob_agg <- readRDS("../data/neon_dob_agg_v2.1.Rds")
neon_dob <- readRDS("../data/phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)

xgrid <- 1/6*seq(1,600)-160 # Based on longitudinal extent of sample data
ygrid <- 1/6*seq(1,360)+15 # Based on latitudinal extent of sample data

sample_data(neon_dob) %>%
  as("data.frame") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(lon, xgrid)], 2),
                            round(ygrid[findInterval(lat, ygrid)], 2), sep="_")) ->
  sample_data(neon_dob)

# take out sample info
dfSampleData <- sample_data(neon_dob)
dfSampleDataAgg <- sample_data(neon_dob_agg)

site <- unique(dfSampleData$Site)
# turn NaN into NA 
dfSampleDataAgg$nitrogenPercent[is.nan(dfSampleDataAgg$nitrogenPercent)] <- NA
dfSampleDataAgg$soilMoisture[is.nan(dfSampleDataAgg$soilMoisture)] <- NA

# create a matrix to store 
# [1] "soilInCaClpH"     "nitrogenPercent"  "organicCPercent"  "soilMoisture"    
# [5] "mat_celsius"      "map_mm"           "temp_seasonality" "prec_seasonality"
matSite <- matrix(ncol = 8,  dimnames = list(NA,c(colnames(dfSampleDataAgg)[c(6,8,9,17,20,21,22,23)])))

# Use Clara's aggregate data to fill in NA
for (i in 1:length(site)){
  temp <- dfSampleData[which(dfSampleData$Site== site[i])]
  vecSiteCoords <- unique(temp$coordsSite)
  if (length(vecSiteCoords) == 1){
    # when one site covers only one grid, take out the value of the grid
    matSite <- rbind(matSite, dfSampleDataAgg[vecSiteCoords,c(6,8,9,17,20,21,22,23)])
  }else if (length(vecSiteCoords) > 1){
    # when covering more than one grid, take the mean of all the grids that are covered
    temp1 <- dfSampleDataAgg[vecSiteCoords,c(6,8,9,17,20,21,22,23)]
    matSite <- rbind(matSite, apply(temp1,2,function(x)mean(x,na.rm = T)))
  }
}

matSite <- matSite[-1,] # delete the first NA row
matSite$nitrogenPercent[is.nan(matSite$nitrogenPercent)] <- NA
matSite$soilMoisture[is.nan(matSite$soilMoisture)] <- NA
rownames(matSite) <- site

# check the number of NA
for (i in 1:8){
  print(sum(is.na(matSite[,i])))
}

# use the sample info (sample_data(neon_dob)) in each site to fill in
for (i in 1:length(site)){
  NAvalue <- is.na(matSite[i,]) # find the NA variable in each site
  if (sum(NAvalue) > 0){
    temp <- dfSampleData[which(dfSampleData$Site== site[i])] # take out the sample info 
    
    # take mean of the environmental variable of all the samples in the site
    matSite[i, colnames(matSite)[which(NAvalue == TRUE)]] <- 
      apply(temp[, colnames(matSite)[which(NAvalue == TRUE)]], 2, function(x) mean(x, na.rm = T))
  }
}

for (i in 1:8){
  print(sum(is.na(matSite[,i])))
}

# some sites have perc.soil.moisture variable instead of soilMoisture
# fill in the soilMoisture variable in these sites
NAmoisture <- which(is.na(matSite$soilMoisture))
for (i in NAmoisture){
  temp <- dfSampleData[which(dfSampleData$Site== site[i])]
  temp$Perc.Soil.Moisture <- as.numeric(temp$Perc.Soil.Moisture)
  matSite[i, "soilMoisture"] <- mean(temp$Perc.Soil.Moisture, na.rm = T)
}

for (i in 1:8){
  print(sum(is.na(matSite[,i])))
}

write.csv(matSite, "FilledInSiteEnv.csv")
# number of NA value
#[1] "soilInCaClpH"     "nitrogenPercent"  "organicCPercent"  "soilMoisture"    
#[5] "mat_celsius"      "map_mm"           "temp_seasonality" "prec_seasonality"
# 0 4 0 6
