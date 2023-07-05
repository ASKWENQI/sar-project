rm(list=ls())
gc()
library(sars)
library(phyloseq)
library(ggplot2)
library(ggtrendline)
library(gridExtra)
library(permute)
library(reshape2)
set.seed(10001)

neon_dob <- readRDS("../data/phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
neon <- subset_samples(neon, !is.na(horizon))
rm(neon_dob)

## 1            ##
## neon dataset ##
d <- sample_data(neon) # sample data data frame

a <- unique(d$Site) # get all the sites
# create matrix to save power
power.increase.c <- power.increase.z <- power.decrease.c <- power.decrease.z <- 
  power.random.c <- power.random.z <- 
  matrix(nrow = 45, ncol = 4,dimnames = list(a, c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )))

models <- vector(length = 20)
# loga.c <- matrix(nrow = 45, ncol = 4)
# loga.z <- matrix(nrow = 45, ncol = 4)
# computing the relationship between number of species and number of samples
for (i in 1:length(a)){
  # take out one site
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, "/", length(a), collapse = ''))# informs the processing
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  otu_tab <- otu_table(neon_sub)
  species_increase <- sort(apply(otu_tab, 1, function(x) sum(x > 0)))
  species_increase_sum <- species_decrease_sum <- species_random <- vector(length = dim1[1])
  species_increase_sum[1] <- species_increase[1]
  species_decrease_sum[1] <- species_increase[dim1[1]]
  sample_seq <- shuffle(c(1:dim1[1]))
  species_random[1] <- species_increase[rownames(otu_tab)[sample_seq[1]]]
  
  for (j in 2:(dim1[1])){
    # take out samples as the sequence in sample_seq
    temp <- colSums(otu_tab[c(names(species_increase)[1:j]),])
    # count species
    species_increase_sum[j] <- sum(temp > 0)
    
    temp <- colSums(otu_tab[c(names(species_increase)[(dim1[1] - j + 1) : dim1[1]]),])
    species_decrease_sum[j] <- sum(temp > 0)
    
    temp <- colSums(otu_tab[c(sample_seq[1:j]),])
    species_random[j] <- sum(temp > 0)
  }
    
  ex_increase <- as.data.frame(cbind("A"=c(1:dim1[1]), "species" = species_increase_sum))
  ex_decrease <- as.data.frame(cbind("A"=c(1:dim1[1]), "species" = species_decrease_sum))
  ex_random <- as.data.frame(cbind("A"=c(1:dim1[1]), "species" = species_random))
  
  temp <- summary(nls(species~c*A^z,ex_increase,start = list(c=1,z=1)))[["coefficients"]]
  power.increase.c[i,] <- temp[1,]
  power.increase.z[i,] <- temp[2,]
  
  temp <- summary(nls(species~c*A^z,ex_decrease,start = list(c=1,z=1)))[["coefficients"]]
  power.decrease.c[i,] <- temp[1,]
  power.decrease.z[i,] <- temp[2,]
  
  temp <- summary(nls(species~c*A^z,ex_random,start = list(c=1,z=1)))[["coefficients"]]
  power.random.c[i,] <- temp[1,]
  power.random.z[i,] <- temp[2,]
  # temp <- summary(nls(species ~ c + z * log(A),ex,start = list(c=1,z=1)))[["coefficients"]]
  # loga.c[i,] <- temp[1,]
  # loga.z[i,] <- temp[2,]
  # temp <- sar_multi(ex)
  # for (k in 1:20){
  #   aic[k] <- temp[[k]][["AIC"]]
  # }
  # models[i] <- names(temp)[which.min(aic)]
}
z_df <- cbind(power.increase.z[,1], power.decrease.z[,1], power.random.z[,1])
colnames(z_df) <- c("increase","decrease","random")
z_df_long <- melt(z_df)
plot(x = z_df_long$Var2, y = z_df_long$value, xlab = "the order of pooling",
     ylab = "z value", 
     main = "Box plot of the z values calculated with different pooling method")
write.table(power.increase.z,"NEONPowerInZ.txt")
write.table(power.increase.z,"NEONPowerDeZ.txt")

# merge samples and linear regression
# neon_site <- merge_samples(neon, "Site")
#d_site <- sample_data(neon_site)
sample_data_neon$nitrogenPercent[is.nan(sample_data_neon$nitrogenPercent)] = NA
sample_data_neon$organicCPercent[is.nan(sample_data_neon$organicCPercent)] = NA
d_site <- aggregate(cbind(sample_data_neon$soilInCaClpH,sample_data_neon$nitrogenPercent,sample_data_neon$organicCPercent,
                          sample_data_neon$soilMoisture, sample_data_neon$map_mm,sample_data_neon$prec_seasonality,sample_data_neon$mat_celsius,
                          sample_data_neon$temp_seasonality),list(sample_data_neon$Site), function(x) mean(x,na.rm = T))
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


## 2   ##
## dob ##
dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")
dob <- subset_samples(dob, !is.na(Site))
rm(neon_dob)

d1 <- sample_data(dob)

a1 <- unique(d1$Site)# be noted that the powering site orders may not be identical to that of generated by the 'merge_samples' function you used at line 128
power.increase.c <- power.increase.z <- power.decrease.c <- power.decrease.z <- 
  power.random.c <- power.random.z <- 
  matrix(nrow = length(a1), ncol = 4,dimnames = list(a, c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )))

for (i in 1:length(a1)){
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, "/", length(a), collapse = ''))# informs the processing
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  otu_tab <- otu_table(neon_sub)
  species_increase <- sort(apply(otu_tab, 1, function(x) sum(x > 0)))
  species_increase_sum <- species_decrease_sum <- species_random <- vector(length = dim1[1])
  species_increase_sum[1] <- species_increase[1]
  species_decrease_sum[1] <- species_increase[dim1[1]]
  sample_seq <- shuffle(c(1:dim1[1]))
  species_random[1] <- species_increase[rownames(otu_tab)[sample_seq[1]]]
  
  for (j in 2:(dim1[1])){
    # take out samples as the sequence in sample_seq
    temp <- colSums(otu_tab[c(names(species_increase)[1:j]),])
    # count species
    species_increase_sum[j] <- sum(temp > 0)
    
    temp <- colSums(otu_tab[c(names(species_increase)[(dim1[1] - j + 1) : dim1[1]]),])
    species_decrease_sum[j] <- sum(temp > 0)
    
    temp <- colSums(otu_tab[c(sample_seq[1:j]),])
    species_random[j] <- sum(temp > 0)
  }
  
  ex_increase <- as.data.frame(cbind("A"=c(1:dim1[1]), "species" = species_increase_sum))
  ex_decrease <- as.data.frame(cbind("A"=c(1:dim1[1]), "species" = species_decrease_sum))
  ex_random <- as.data.frame(cbind("A"=c(1:dim1[1]), "species" = species_random))
  
  temp <- summary(nls(species~c*A^z,ex_increase,start = list(c=1,z=1)))[["coefficients"]]
  power.increase.c[i,] <- temp[1,]
  power.increase.z[i,] <- temp[2,]
  
  temp <- summary(nls(species~c*A^z,ex_decrease,start = list(c=1,z=1)))[["coefficients"]]
  power.decrease.c[i,] <- temp[1,]
  power.decrease.z[i,] <- temp[2,]
  
  temp <- summary(nls(species~c*A^z,ex_random,start = list(c=1,z=1)))[["coefficients"]]
  power.random.c[i,] <- temp[1,]
  power.random.z[i,] <- temp[2,]
  # temp <- summary(nls(species ~ c + z * log(A),ex,start = list(c=1,z=1)))[["coefficients"]]
  # loga.c[i,] <- temp[1,]
  # loga.z[i,] <- temp[2,]
  # temp <- sar_multi(ex)
  # for (k in 1:20){
  #   aic[k] <- temp[[k]][["AIC"]]
  # }
  # models[i] <- names(temp)[which.min(aic)]
}

z_df <- cbind(power.increase.z[,1], power.decrease.z[,1], power.random.z[,1])
colnames(z_df) <- c("increase","decrease","random")
z_df_long <- melt(z_df)
plot(x = z_df_long$Var2, y = z_df_long$value, xlab = "the order of pooling",
     ylab = "z value", 
     main = "Box plot of the z values calculated with different pooling method")

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
