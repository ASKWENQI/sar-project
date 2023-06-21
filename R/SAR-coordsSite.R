rm(list=ls())
library(dplyr)
gc()

# loading and filtering
library(phyloseq)
neon_dob <- readRDS("../Data/phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, !is.na(Site))
d <- sample_data(neon_dob)[,c(6,8,9,17,20:23)]
# d1 <- d[which(d$Project=="DoB"),c(5,6,8,9,17,20:23)]
# d1 <- d1[which(d1$horizon=="O"),-1]
ggcorrplot(cor(d, use = "pairwise.complete.obs"),type = "upper",lab = T, p.mat = cor_pmat(d1))

xgrid <- 1/6*seq(1,600)-160 # Based on longitudinal extent of sample data
ygrid <- 1/6*seq(1,360)+15 # Based on latitudinal extent of sample data
sample_data(neon_dob) %>%
  as("data.frame") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(lon, xgrid)], 2),
                            round(ygrid[findInterval(lat, ygrid)], 2), sep="_")) ->
  sample_data(neon_dob)
#rm(neon_dob)

# neon dataset
d <- sample_data(neon_dob) # sample data data frame

a <- unique(d$coordsSite) # get all the sites
# create matrix to save result
result.c <- matrix(nrow = 128, ncol = 4)
result.z <- matrix(nrow = 128, ncol = 4)

# computing the relationship between number of species and number of samples
for (i in 109:length(a)){
  # take out one site
  neon_dob_sub <- subset_samples(neon_dob, coordsSite==a[i])
  dim1 <- dim(otu_table(neon_dob_sub)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  for (j in 1:(dim1[1])){ 
    
    # randomly sample j samples in the site 
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(neon_dob_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
    
    # compute number of species
    species[j] <- sum(otu_table(temp)["TRUE"] > 0)
  }
  ex <- as.data.frame(cbind(species, "A"=c(1:dim1[1])))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  result.c[i,] <- temp[1,]
  result.z[i,] <- temp[2,]
}
colnames(result.z) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
colnames(result.c) <- colnames(summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]])
rownames(result.c) <- a
rownames(result.z) <- a

result.z <- na.omit(result.z)
result.c <- na.omit(result.c)
result.z

# merge samples and linear regression
neon_dob_agg <- readRDS("neon_dob_agg_v2.1.Rds")
d_site <- sample_data(neon_dob_agg)
d_site$nitrogenPercent[is.nan(d_site$nitrogenPercent)]<-NA
d_site <- d_site[row.names(result.z),c(6,8,9,17,20:23)]
d_site <- data.frame(scale(d_site))# standardization seems to have no effect on significance

summary(lm(result.z[,1] ~ d_site$soilInCaClpH + d_site$organicCPercent + d_site$soilMoisture + d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality))

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

