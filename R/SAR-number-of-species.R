rm(list=ls())
gc()
library(sars)
# loading and filtering
library(phyloseq)
neon_dob <- readRDS("../data/phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
neon <- subset_samples(neon, !is.na(horizon))
#rm(neon_dob)

# neon dataset
d <- sample_data(neon) # sample data data frame

a <- unique(d$Site) # get all the sites
# create matrix to save power
power.c <- power.z <- matrix(nrow = length(a), ncol = 4,
         dimnames = list(a, c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )))

models <- vector(length = 20)

# species_sum <- vector(length = length(a))
# loga.c <- matrix(nrow = 45, ncol = 4)
# loga.z <- matrix(nrow = 45, ncol = 4)
# computing the relationship between number of species and number of samples
for (i in 1:length(a)){
  # take out one site
    cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  for (j in 1:(dim1[1])){ 
    
    # randomly sample j samples in the site 
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(neon_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
    
    # compute number of species
    species[j] <- sum(otu_table(temp)["TRUE"] > 0)
  }
  ex <- as.data.frame(cbind("A"=c(1:dim1[1]),species))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  power.c[i,] <- temp[1,]
  power.z[i,] <- temp[2,]
  
  #species_sum[i] <- sum(colSums(otu_table(neon_sub)) > 0)
  #species[length(species)] == species_sum[i]
  
  # temp <- summary(nls(species ~ c + z * log(A),ex,start = list(c=1,z=1)))[["coefficients"]]
  # loga.c[i,] <- temp[1,]
  # loga.z[i,] <- temp[2,]
  # temp <- sar_multi(ex)
  # for (k in 1:20){
  #   aic[k] <- temp[[k]][["AIC"]]
  # }
  # models[i] <- names(temp)[which.min(aic)]
}

power.z
power.c
write.table(power.z,"neon.z.txt")
write.table(power.c,"neon.c.txt")
# merge samples and linear regression
# neon_site <- merge_samples(neon, "Site")
#d_site <- sample_data(neon_site)
d$nitrogenPercent[is.nan(d$nitrogenPercent)] = NA
d$organicCPercent[is.nan(d$organicCPercent)] = NA
d_site <- aggregate(cbind(d$soilInCaClpH,d$nitrogenPercent,d$organicCPercent,
                          d$soilMoisture, d$map_mm,d$prec_seasonality,d$mat_celsius,
                          d$temp_seasonality),list(d$Site), function(x) mean(x,na.rm = T))
d_site$V2[is.nan(d_site$V2)] = NA
d_site$V3[is.nan(d_site$V3)] = NA
row.names(d_site) <- d_site$Group.1
d_site <- d_site[,-1]
colnames(d_site) <- c("soilInCaClpH","nitrogenPercent","organicCPercent",
                      "soilMoisture", "map_mm",'prec_seasonality',"mat_celsius",
                      "temp_seasonality")
# standardization seems to have no effect on significance
d_site <- d_site[a,]
d_site <- as.data.frame(scale(d_site))
summary(lm(power.z[,1] ~ d_site$soilInCaClpH + d_site$soilMoisture + d_site$mat_celsius +
     d_site$map_mm + d_site$temp_seasonality))

# data with no missing values
summary(lm(scale(power.z[,1]) ~ d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality + d_site$prec_seasonality))

library(ggplot2)
library(ggtrendline)
library(gridExtra)
p1 <- ggtrendline(d_site$soilInCaClpH, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilInCaClpH, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilInCaClpH")+ylab("z")

p2 <- ggtrendline(d_site$soilMoisture, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilMoisture, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilMoisture")+ylab("z")

p3 <- ggtrendline(d_site$mat_celsius, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$mat_celsius, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("mat_celsius")+ylab("z")

p4 <- ggtrendline(d_site$map_mm, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$map_mm, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("map_mm")+ylab("z")

p5 <- ggtrendline(d_site$temp_seasonality, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$temp_seasonality, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("temp_seasonality")+ylab("z") 

p6 <- ggtrendline(d_site$prec_seasonality, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$prec_seasonality, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("prec_seasonality")+ylab("z") 

grid.arrange(p1,p2,p3,p4,p5,p6,ncol = 3)

# dob
dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")
dob <- subset_samples(dob, !is.na(Site))
# dob <- subset_samples(dob, horizon %in% c("O", "M"))
rm(neon_dob)

d1 <- sample_data(dob)

a1 <- unique(d1$Site)# be noted that the powering site orders may not be identical to that of generated by the 'merge_samples' function you used at line 128
power.c <- power.z <- matrix(nrow = length(a1), ncol = 4,
                             dimnames = list(a1, c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )))
species_sum <- vector(length = length(a1))
for (i in 1:length(a1)){
  # take out one site
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))
  dob_sub <- subset_samples(dob, Site==a1[i])
  dim1 <- dim(otu_table(dob_sub)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  for (j in 1:(dim1[1])){ 
    
    # randomly sample j samples in the site 
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(dob_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
    
    # compute number of species
    species[j] <- sum(otu_table(temp)["TRUE"] > 0)
  }
  ex <- as.data.frame(cbind(species, "A"=c(1:dim1[1])))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  power.c[i,] <- temp[1,]
  power.z[i,] <- temp[2,]
  
  species_sum[i] <- species[length(species)]
  # temp <- summary(nls(species ~ c + z * log(A),ex,start = list(c=1,z=1)))[["coefficients"]]
  # loga.c[i,] <- temp[1,]
  # loga.z[i,] <- temp[2,]
  # temp <- sar_multi(ex)
  # aic <- vector(length = 20)
  # for (k in 1:20){
  #   aic[k] <- temp[[k]][["AIC"]]
  # }
  # models[i] <- names(temp)[which.min(aic)]
}

colnames(power.z) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
power.z[,4] > 0.05
hist(power.z[,1])
summary(power.z[,1])

# merge samples and linear regression
write.table(power.z,"dob.z.txt")
write.table(power.c,"dob.c.txt")
# merge samples and linear regression
# dob_site <- merge_samples(dob, "Site")
#d_site <- sample_data(dob_site)
d1$nitrogenPercent[is.nan(d1$nitrogenPercent)] = NA
d1$organicCPercent[is.nan(d1$organicCPercent)] = NA
d1$map_mm[is.nan(d1$map_mm)] = NA

d_site <- aggregate(cbind(d1$soilInCaClpH,d1$nitrogenPercent,d1$organicCPercent,
                          d1$soilMoisture, d1$map_mm,d1$prec_seasonality,d1$mat_celsius,
                          d1$temp_seasonality),list(d1$Site), function(x) mean(x,na.rm = T))
d_site$V2[is.nan(d_site$V2)] = NA
d_site$V3[is.nan(d_site$V3)] = NA
d_site$V4[is.nan(d_site$V4)] = NA
row.names(d_site) <- d_site$Group.1
d_site <- data.frame(scale(d_site[,-1]))
colnames(d_site) <- c("soilInCaClpH","nitrogenPercent","organicCPercent",
                      "soilMoisture", "map_mm",'prec_seasonality',"mat_celsius",
                      "temp_seasonality")
# standardization seems to have no effect on significance
d_site <- d_site[a1,]
summary(lm(power.z[,1] ~ d_site$soilInCaClpH + d_site$soilMoisture + d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality))

# data with no missing values
summary(lm(scale(power.z[,1]) ~ d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality + d_site$prec_seasonality))

summary(lm(power.z[,1] ~ d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality + d_site$prec_seasonality))

p1 <- ggtrendline(d_site$soilInCaClpH, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilInCaClpH, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilInCaClpH")+ylab("z")

p2 <- ggtrendline(d_site$soilMoisture, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilMoisture, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilMoisture")+ylab("z")

p3 <- ggtrendline(d_site$mat_celsius, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$mat_celsius, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("mat_celsius")+ylab("z")

p4 <- ggtrendline(d_site$map_mm, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$map_mm, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("map_mm")+ylab("z")

p5 <- ggtrendline(d_site$temp_seasonality, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$temp_seasonality, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("temp_seasonality")+ylab("z") 

p6 <- ggtrendline(d_site$prec_seasonality, power.z[,1],linecolor="red")+
  geom_point(aes(d_site$prec_seasonality, power.z[,1]),color="black",pch=21,fill='white')+
  xlab("prec_seasonality")+ylab("z") 

grid.arrange(p1,p2,p3,p4,p5,p6,ncol = 3)
