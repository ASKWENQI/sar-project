library(permute)
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

#####
#Dob#

dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")
dob <- subset_samples(dob, !is.na(Site))
rm(neon_dob)

d1 <- sample_data(dob)

a1 <- unique(d1$Site)
result.c <- matrix(nrow = 68, ncol = 4)
result.z <- matrix(nrow = 68, ncol = 4)

for (i in 1:length(a)){
  # take out one site
  dob_sub <- subset_samples(dob, Site==a[i])
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