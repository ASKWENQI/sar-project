# 3. Prepare guild data ------------------------------------------------------
library(dplyr)
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
# First see if genus is unique, or if higher taxonomic ranks are necessary for matching
tax_table(neon_dob) %>%
  as("matrix") %>% as.data.frame() %>%
  group_by(genus) %>%
  summarise(in_n_families = n_distinct(family)) %>%
  arrange(desc(in_n_families))
# Genus is unambiguous!

# Load FungalTraits
ft <- read.csv("../data/FungalTraits_1.2_ver_16Dec_2020.csv", na.strings="")
ft <- ft[which(!duplicated(ft)),]
names(ft) <- tolower(names(ft))

# Yes! Thus, safe to remove duplicates.
ft <- ft[-which(duplicated(ft$genus)),]

# Create new column(s), representing broader groupings
guild2 <- rep(NA, nrow(ft))
guild2[grep("ectomycorrhizal", ft$primary_lifestyle)] <- "EM"
guild2[grep("arbuscular_mycorrhizal", ft$primary_lifestyle)] <- "AM"
guild2[grep("saprotroph", ft$primary_lifestyle)] <- "saprotroph"
guild2[grep("pathogen", ft$primary_lifestyle)] <- "pathogen"
guild2[grep("symbiotroph", ft$primary_lifestyle)] <- "symbiotroph"
guild2[grep("symbiont", ft$primary_lifestyle)] <- "symbiont"
guild2[grep("parasite", ft$primary_lifestyle)] <- "parasite"
guild2[grep("lichenized", ft$primary_lifestyle)] <- "lichenized"
guild2[grep("epiphyte", ft$primary_lifestyle)] <- "epiphyte"
guild2[grep("endophyte", ft$primary_lifestyle)] <- "endophyte"
guild2[is.na(guild2) & ft$primary_lifestyle!="unspecified"] <- "other"
ft$guild2 <- guild2
# table(ft$guild2, useNA="ifany")

# Use a subset of trait data for this analysis
ft <- ft %>%
  dplyr::select(
    class:genus, primary_lifestyle, secondary_lifestyle,
    ectomycorrhiza_exploration_type_template,
    growth_form_template, guild2) %>%
  dplyr::filter(!is.na(primary_lifestyle))

d <- sample_data(neon) # sample data data frame

a <- unique(d$Site) # get all the sites
# create matrix to save power
power.c <- power.z <- matrix(nrow = length(a), ncol = 4,
                             dimnames = list(a, c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )))

species_sum <- vector(length = length(a))
# loga.c <- matrix(nrow = 45, ncol = 4)
# loga.z <- matrix(nrow = 45, ncol = 4)
# computing the relationship between number of species and number of samples
for (i in 1:length(a)){
  # take out one site
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  neon_sub <- subset_samples(neon_dob, Site==a[i])
  temp <- apply(tax_table(neon_sub), 1, function(x) {
    if (x[["genus"]] %in% ft$genus){
      return(ft$guild2[which(ft$genus == x[["genus"]])])
    }else{
      return (NA)
    }
  })
  tax_table(neon_sub) <- cbind(tax_table(neon_sub), "guild2" = temp)
  neon_guild <- subset_taxa(neon_sub, guild2 == "EM")
  # neon_guild <- subset_taxa(neon_sub, guild2 == "saprotroph")
  
  dim1 <- dim(otu_table(neon_guild)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  for (j in 1:(dim1[1])){ 
    
    # randomly sample j samples in the site 
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(neon_guild, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
    
    # compute number of species
    species[j] <- sum(otu_table(temp)["TRUE"] > 0)
  }
  ex <- as.data.frame(cbind("A"=c(1:dim1[1]),species))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  power.c[i,] <- temp[1,]
  power.z[i,] <- temp[2,]
}

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

#DoB
dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")
dob <- subset_samples(dob, !is.na(Site))
dob <- subset_samples(dob, horizon %in% c("O", "M"))
rm(neon_dob)

d1 <- sample_data(dob)

a1 <- unique(d1$Site)# be noted that the powering site orders may not be identical to that of generated by the 'merge_samples' function you used at line 128
power.c <- power.z <- matrix(nrow = length(a1), ncol = 4,
                             dimnames = list(a1, c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )))
species_sum <- vector(length = length(a1))

for (i in 1:length(a1)){
  # take out one site
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  dob_sub <- subset_samples(dob, Site==a1[i])
  temp <- apply(tax_table(dob_sub), 1, function(x) {
    if (x[["genus"]] %in% ft$genus){
      return(ft$guild2[which(ft$genus == x[["genus"]])])
    }else{
      return (NA)
    }
  })
  tax_table(dob_sub) <- cbind(tax_table(dob_sub), "guild2" = temp)
  dob_guild <- subset_taxa(dob_sub, guild2 == "EM")
  # dob_guild <- subset_taxa(dob_sub, guild2 == "saprotroph")
  
  dim1 <- dim(otu_table(dob_guild)) # the number of samples in one site
  species <- vector(length = dim1[1]) # create a vector to save diversity
  for (j in 1:(dim1[1])){ 
    
    # randomly sample j samples in the site 
    flag <- rep(FALSE, dim1[1])
    flag[sample(1:dim1[1], j)] <- TRUE
    temp <- merge_samples(dob_guild, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
    
    # compute number of species
    species[j] <- sum(otu_table(temp)["TRUE"] > 0)
  }
  ex <- as.data.frame(cbind("A"=c(1:dim1[1]),species))
  temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
  power.c[i,] <- temp[1,]
  power.z[i,] <- temp[2,]
  
  species_sum[i] <- species[length(species)]
}

write.table(power.z, "../result/sapr.power.z2.txt")
plot(power.z[,1], species_sum)
summary(lm(power.z[,1] ~ species_sum + I(species_sum ^ 2)))
