rm(list=ls())
library(phyloseq)
library(doParallel)

neon_dob <- readRDS("../data/phylo_V3.1.RDS")
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
rm(neon_dob)
# neon <- subset_samples(neon, horizon == "O")
# neon <- subset_samples(neon, horizon == "M")
neon <- subset_samples(neon, !is.na(lon) & !is.na(lat))
neon <- subset_taxa(neon, taxa_sums(neon) > 0)
neon <- subset_samples(neon, sample_sums(neon) > 0)
neon <- subset_samples(neon, !is.na(horizon))

# neon dataset
d <- sample_data(neon) # sample data data frame

a <- unique(d$Site) # get all the sites
# create matrix to save power

times = 30
power.z <- vector("list", length(a))

# computing the relationship between number of species and number of samples
otu_tab <- otu_table(neon)

for (i in 1:length(a)){
  # take out one site
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  
  neon_sub <- subset_samples(neon, Site==a[i])
  dim1 <- dim(otu_table(neon_sub)) # the number of samples in one site
  if (dim1[1] >= 30){
    cl <- makeCluster(3)
    registerDoParallel(cl)
    
    power.z[[i]] <- foreach(i = 1:times, .combine = "cbind", .packages = c("phyloseq")) %dopar% {
      species <- vector(length = 30) # create a vector to save diversity
      for (j in 1:30){ 
        
        # randomly sample j samples in the site 
        flag <- rep(FALSE, dim1[1])
        flag[sample(1:dim1[1], j)] <- TRUE
        temp <- merge_samples(neon_sub, flag, function(x) mean(x, na.rm = TRUE)) # the j samples aggregated by mean
        
        # compute number of species
        species[j] <- sum(otu_table(temp)["TRUE"] > 0)
      }
      ex <- as.data.frame(cbind("A"=c(1:30),species))
      temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
      
      return(c(temp[2,1], temp[2,4]))
    }
    stopCluster(cl)
  }
}

power.z.matrix <- matrix(nrow = length(power.z), ncol = times, dimnames = list(a, NULL)) 
pvalue.matrix <- matrix(nrow = length(power.z), ncol = times, dimnames = list(a, NULL)) 

for (i in 1:length(power.z)){
  if (!is.null(power.z[[i]])){
    power.z.matrix[i,] <- power.z[[i]][1,]
    pvalue.matrix[i,] <- power.z[[i]][2,]
  }
}

write.table(power.z.matrix,"D_power.z.NEON.txt")
write.table(pvalue.matrix,"D_pvalue.NEON.txt")

write.table(power.z.matrix,"D_power.z.NEON.O.txt")
write.table(pvalue.matrix,"D_pvalue.NEON.O.txt")

write.table(power.z.matrix,"D_power.z.NEON.M.txt")
write.table(pvalue.matrix,"D_pvalue.NEON.M.txt")