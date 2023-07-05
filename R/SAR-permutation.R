rm(list=ls())
library(phyloseq)
library(doParallel)
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

times = 25
power.c <- power.z <- matrix(nrow = 45, ncol = times,
       dimnames = list(a, c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )))

# computing the relationship between number of species and number of samples
otu_tab <- otu_table(neon)

for (i in 1:length(a)){
  # take out one site
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  cl <- makeCluster(3)
  registerDoParallel(cl)
  power.z[i,] <- foreach(i = 1:times, .combine = "c") %dopar% {
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
    return(temp[2,1])
  }
  stopCluster(cl)
}

write.table(power.z,"permutation.power.z.txt")