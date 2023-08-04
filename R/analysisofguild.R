listPowerZEM <- readRDS("E:/power.z.guild.per.neon.RDS")
checkPValue <- function(list_of_z){
  flag <- vector()
  for (i in 1:length(list_of_z)){
    if (sum(list_of_z[[i]][2, ] > 0.05) > 0)
      flag <- c(flag, i)
  }
  if (length(flag) == 0){
    print("all z values are significant")
  }else{
    return(flag)
  }
}
checkPValue(listPowerZEM)

matPowerZEM <- matrix(nrow = length(listPowerZEM), ncol = 25, dimnames = list(names(listPowerZEM), NULL))
for (i in 1:length(listPowerZEM))
  matPowerZEM[i,] <- listPowerZEM[[i]][1, ]
powerZEMMean <- apply(matPowerZEM, 1, mean)
deviation <- apply(matPowerZEM, 1, function(x){sum((x - mean(x))^2)})
boxplot(t(matPowerZEM))

# listPowerZ <- readRDS("E:/power.z.dob.per.RDS")
# 
# #rm(neon_dob)
# 
# checkPValue(listPowerZ)
# 
# matPowerZ <- matrix(nrow = length(listPowerZ), ncol = 25, dimnames = list(names(listPowerZ), NULL))
# for (i in 1:length(listPowerZ))
#   matPowerZ[i,] <- listPowerZ[[i]][1, ]
# power_z_mean <- apply(matPowerZ, 1, mean)
# deviation <- apply(matPowerZ, 1, function(x){sum((x - mean(x))^2)})

matPowerZ <- read.table("E:/power.z.per.txt")

pValue <- vector(length = length(listPowerZEM))
for (i in 1:length(listPowerZEM)){
  pValue[i] <- t.test(matPowerZ[i,], matPowerZEM[i,])[["p.value"]]
}
pValue < 0.05
sum(pValue < 0.05)

d <- read.table("E:/sample_data_neon_dob.txt")
d <- d[which(d$Project == "NEON"),]

a <- unique(d$Site) # get all the sites
# create matrix to save power
names(deviation) <- a

species_sum <- vector(length = length(a))
names(species_sum) <- a
EMSum <- vector(length = length(a))
for (i in 1:length(a)){
  # take out one site
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  species_sum[i] <- sum(colSums(otu_table(neon_sub)) > 0)
  
  temp <- apply(tax_table(neon_sub), 1, function(x) {
    if (x[["genus"]] %in% ft$genus){
      return(ft$guild2[which(ft$genus == x[["genus"]])])
    }else{
      return (NA)
    }
  })
  tax_table(neon_sub) <- cbind(tax_table(neon_sub), "guild2" = temp)
  neon_guild <- subset_taxa(neon_sub, guild2 == "EM")
  EMSum[i] <- sum(colSums(otu_table(neon_guild)) > 0)
}

sample_sum <- data.frame(aggregate(d$geneticSampleID, list(d$Site), function(x) length(unique(x))), row.names = 1)
sample_sum <- sample_sum[a,]

library(ggplot2)
library(ggtrendline)
library(gridExtra)
p1 <- ggtrendline(EMSum, powerZEMMean,linecolor="red")+
  geom_point(aes(EMSum, powerZEMMean),color="black",pch=21,fill='white')+
  xlab("number of species")+ylab("z")

p2 <- ggtrendline(sample_sum, powerZEMMean,linecolor="red")+
  geom_point(aes(sample_sum, powerZEMMean),color="black",pch=21,fill='white')+
  xlab("number of samples")+ylab("z")

grid.arrange(p1,p2,ncol = 2)

library(ggplot2)
library(ggtrendline)
library(gridExtra)
p1 <- ggtrendline(EMSum[1:43], powerZEMMean[1:43],linecolor="red")+
  geom_point(aes(EMSum[1:43], powerZEMMean[1:43]),color="black",pch=21,fill='white')+
  xlab("number of species")+ylab("z of O & M horizon")

p3 <- ggtrendline(EMSum[44:68], powerZEMMean[44:68],linecolor="red")+
  geom_point(aes(EMSum[44:68], powerZEMMean[44:68]),color="black",pch=21,fill='white')+
  xlab("number of species")+ylab("z of Oh & Ah horizon")

ggtrendline(species_sum[44:68], powerZEMMean[44:68],linecolor="red", model = "line3P")+
  geom_point(aes(species_sum[44:68], powerZEMMean[44:68]),color="black",pch=21,fill='white')+
  xlab("number of species")+ylab("z of O & M horizon")

p2 <- ggtrendline(sample_sum[1:43], powerZEMMean[1:43],linecolor="red")+
  geom_point(aes(sample_sum[1:43], powerZEMMean[1:43]),color="black",pch=21,fill='white')+
  xlab("number of samples")+ylab("z of O & M horizon")

p4 <- ggtrendline(sample_sum[44:68], powerZEMMean[44:68],linecolor="red")+
  geom_point(aes(sample_sum[44:68], powerZEMMean[44:68]),color="black",pch=21,fill='white')+
  xlab("number of samples")+ylab("z of Oh & Ah horizon")

grid.arrange(p1,p2,p3,p4,ncol = 2)

powerZMeanCombine <- cbind(c(power_z_mean, powerZEMMean),
                           "group" = c(rep("all species", times = 68), rep("EM"),times = 68)

                           
df <- data.frame(cbind(sample_sum, EMSum, c(rep("O & M", times = 43), rep("Oh & Ah", times = 68-43))))
names(df) <- c("sample", "species", "horizon")
df$horizon[which(df$horizon == TRUE)] = "O & M"
df$horizon[which(df$horizon == 0)] = "Oh & Ah"
df$horizon <- as.factor(df$horizon)
df$sample <- as.numeric(df$sample)
df$species <- as.numeric(df$species)
ggplot(df, aes(x = sample, y = species, color = horizon))+
  geom_point()+
  geom_smooth(method =  "lm", show.legend = TRUE)

power_z_mean_OM <- powerZEMMean[1:43]
power_z_mean_OA <- powerZEMMean[44:68]


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
d_site_OM <- d_site[1:43,]
d_site_OA <- d_site[44:68,]

d_site <- d_site_OA
temp <- powerZEMMean

d_site <- as.data.frame(scale(d_site))
summary(lm(temp ~ d_site$soilInCaClpH + d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality + d_site$prec_seasonality))

# data with no missing values
summary(lm(scale(temp) ~ d_site$mat_celsius +
             d_site$map_mm + d_site$temp_seasonality + d_site$prec_seasonality))

library(ggplot2)
library(ggtrendline)
library(gridExtra)
p1 <- ggtrendline(d_site$soilInCaClpH, temp,linecolor="red")+
  geom_point(aes(d_site$soilInCaClpH, temp),color="black",pch=21,fill='white')+
  xlab("soilInCaClpH")+ylab("z")

p2 <- ggtrendline(d_site$soilMoisture, temp,linecolor="red")+
  geom_point(aes(d_site$soilMoisture, temp),color="black",pch=21,fill='white')+
  xlab("soilMoisture")+ylab("z")

p3 <- ggtrendline(d_site$mat_celsius, temp,linecolor="red")+
  geom_point(aes(d_site$mat_celsius, temp),color="black",pch=21,fill='white')+
  xlab("mat_celsius")+ylab("z")

p4 <- ggtrendline(d_site$map_mm, temp,linecolor="red")+
  geom_point(aes(d_site$map_mm, temp),color="black",pch=21,fill='white')+
  xlab("map_mm")+ylab("z")

p5 <- ggtrendline(d_site$temp_seasonality, temp,linecolor="red")+
  geom_point(aes(d_site$temp_seasonality, temp),color="black",pch=21,fill='white')+
  xlab("temp_seasonality")+ylab("z") 

p6 <- ggtrendline(d_site$prec_seasonality, temp,linecolor="red")+
  geom_point(aes(d_site$prec_seasonality, temp),color="black",pch=21,fill='white')+
  xlab("prec_seasonality")+ylab("z") 

grid.arrange(p1,p2,p3,p4,p5,p6,ncol = 3)


p2 <- ggtrendline(d_site$soilMoisture[which(d_site$soilMoisture<2)], temp[which(d_site$soilMoisture<2)],linecolor="red")+
  geom_point(aes(d_site$soilMoisture[which(d_site$soilMoisture<2)], temp[which(d_site$soilMoisture<2)]),color="black",pch=21,fill='white')+
  xlab("soilMoisture")+ylab("z")

p3 <- ggtrendline(d_site$map_mm[which(d_site$map_mm < 2)], temp[which(d_site$map_mm < 2)],linecolor="red")+
  geom_point(aes(d_site$map_mm[which(d_site$map_mm < 2)], temp[which(d_site$map_mm < 2)]),color="black",pch=21,fill='white')+
  xlab("map_mm")+ylab("z")

grid.arrange(p2,p3,ncol = 2)

df <- cbind(temp, d_site)
library(ggplot2)
p1 <- ggtrendline(d_site$mat_celsius, temp,linecolor="red", model = "line3P")+
  geom_point(aes(d_site$mat_celsius, temp),color="black",pch=21,fill='white')+
  xlab("mat_celsius")+ylab("z") 
p2 <- ggtrendline(d_site$map_mm, temp,linecolor="red", model = "line3P")+
  geom_point(aes(d_site$prec_seasonality, temp),color="black",pch=21,fill='white')+
  xlab("prec_seasonality")+ylab("z") 

grid.arrange(p1,p2,ncol = 2)

# 50 times bootstraps
times = 25
library(parallel)
foreach(i = c(1:20), .packages = c("glmnet",))