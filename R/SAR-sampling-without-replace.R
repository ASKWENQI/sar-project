library(permute)
set.seed(100001)
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

a <- unique(d$Site) # get all the sites
# create matrix to save result
result.c <- matrix(nrow = 45, ncol = 4)
result.z <- matrix(nrow = 45, ncol = 4)
for (i in 1:length(a)){
  # take out one site
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  
  # shuffled sequence e.g. sample_seq = c(4,2,3,1,5) for 5 samples in one site
  # take out samples number 4,2,3,1,5 and add them one by one
  sample_seq <- shuffle(c(1:dim1[1]))
  
  # take out otu_tab for each site
  otu_tab <- otu_table(neon_sub)
  otu_tab <- matrix(otu_tab, nrow = dim1[1], byrow = TRUE)
  
  # j = 1 as a special case because line 23
  # 'temp <- colSums(otu_tab[c(sample_seq[1:j]),])' 
  # would give error
  species[1] <- sum(otu_tab[sample_seq[1],] > 0)
  
  for (j in 2:(dim1[1])){ 
    # take out samples as the sequence in sample_seq
    temp <- colSums(otu_tab[c(sample_seq[1:j]),])
    # count species
    species[j] <- sum(temp > 0)
  }
  ex <- as.data.frame(cbind(species, "A"=c(1:dim1[1])))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  result.c[i,] <- temp[1,]
  result.z[i,] <- temp[2,]
}
colnames(result.z) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
colnames(result.c) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
result.z
result.c

hist(result.z[,1]) # many z bigger than 1

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
summary(lm(result.z[,1] ~ d_site$soilInCaClpH + d_site$soilMoisture + d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality + d_site$prec_seasonality))

# data with no missing values
summary(lm(scale(result.z[,1]) ~ d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality + d_site$prec_seasonality))

library(ggplot2)
library(ggtrendline)
library(gridExtra)
p1 <- ggtrendline(d_site$soilInCaClpH, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilInCaClpH, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilInCaClpH")+ylab("z")

p2 <- ggtrendline(d_site$soilMoisture, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilMoisture, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilMoisture")+ylab("z")

p3 <- ggtrendline(d_site$mat_celsius, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$mat_celsius, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("mat_celsius")+ylab("z")

p4 <- ggtrendline(d_site$map_mm, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$map_mm, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("map_mm")+ylab("z")

p5 <- ggtrendline(d_site$temp_seasonality, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$temp_seasonality, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("temp_seasonality")+ylab("z") 

p6 <- ggtrendline(d_site$prec_seasonality, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$prec_seasonality, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("prec_seasonality")+ylab("z") 

grid.arrange(p1,p2,p3,p4,p5,p6,ncol = 3)

#####
#Dob#

dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")
dob <- subset_samples(dob, !is.na(Site))
rm(neon_dob)

d1 <- sample_data(dob)

a1 <- unique(d1$Site)
result.c <- matrix(nrow = 68, ncol = 4)
result.z <- matrix(nrow = 68, ncol = 4)

for (i in 1:length(a1)){
  # take out one site
  dob_sub <- subset_samples(dob, Site==a1[i])
  dim1 <- dim(otu_table(dob_sub)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  
  # shuffled sequence e.g. sample_seq = c(4,2,3,1,5) for 5 samples in one site
  # take out samples number 4,2,3,1,5 and add them one by one
  sample_seq <- shuffle(c(1:dim1[1]))
  
  # take out otu_tab for each site
  otu_tab <- otu_table(dob_sub)
  otu_tab <- matrix(otu_tab, nrow = dim1[1], byrow = TRUE)
  
  # j = 1 as a special case because line 23
  # 'temp <- colSums(otu_tab[c(sample_seq[1:j]),])' 
  # would give error
  species[1] <- sum(otu_tab[sample_seq[1],] > 0)
  
  for (j in 2:(dim1[1])){ 
    # take out samples as the sequence in sample_seq
    temp <- colSums(otu_tab[c(sample_seq[1:j]),])
    # count species
    species[j] <- sum(temp > 0)
  }
  ex <- as.data.frame(cbind(species, "A"=c(1:dim1[1])))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  result.c[i,] <- temp[1,]
  result.z[i,] <- temp[2,]
}
colnames(result.z) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
colnames(result.c) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
result.z
result.c

hist(result.z[,1]) # many z bigger than 1

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
summary(lm(result.z[,1] ~ d_site$soilInCaClpH + d_site$soilMoisture + d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality))

# data with no missing values

summary(lm(result.z[,1] ~ d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality + d_site$prec_seasonality))

p1 <- ggtrendline(d_site$soilInCaClpH, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilInCaClpH, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilInCaClpH")+ylab("z")

p2 <- ggtrendline(d_site$soilMoisture, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$soilMoisture, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("soilMoisture")+ylab("z")

p3 <- ggtrendline(d_site$mat_celsius, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$mat_celsius, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("mat_celsius")+ylab("z")

p4 <- ggtrendline(d_site$map_mm, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$map_mm, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("map_mm")+ylab("z")

p5 <- ggtrendline(d_site$temp_seasonality, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$temp_seasonality, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("temp_seasonality")+ylab("z") 

p6 <- ggtrendline(d_site$prec_seasonality, result.z[,1],linecolor="red")+
  geom_point(aes(d_site$prec_seasonality, result.z[,1]),color="black",pch=21,fill='white')+
  xlab("prec_seasonality")+ylab("z") 

grid.arrange(p1,p2,p3,p4,p5,p6,ncol = 3)
