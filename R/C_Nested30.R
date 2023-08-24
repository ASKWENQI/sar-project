rm(list=ls())
library(phyloseq)
library(foreach)
library(doParallel)
library(permute)
neon_dob <- readRDS("../data/phylo_V3.1.RDS")
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
rm(neon_dob)
# neon <- subset_samples(neon, horizon == "O")
# neon <- subset_samples(neon, horizon == "M")
neon <- subset_samples(neon, !is.na(lon) & !is.na(lat))
neon <- subset_taxa(neon, taxa_sums(neon) > 0)
neon <- subset_samples(neon, sample_sums(neon) > 0)
neon <- subset_samples(neon, !is.na(horizon))
gc()
# neon dataset
d <- sample_data(neon) # sample data data frame

a <- unique(d$Site) # get all the sites
# create matrix to save power

times = 30
power.z <- vector("list", length(a))

# computing the relationship between number of species and number of samples
otu_tab <- otu_table(neon)

for (i in 28:length(a)){
  # take out one site
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, collapse = ''))# informs the processing
  
  otu_tab <- otu_table(subset_samples(neon, Site==a[i]))
  dim1 <- dim(otu_tab) # the number of samples in one site
  if (dim1[1] >= 30){
    cl <- makeCluster(3)
    registerDoParallel(cl)
    
    power.z[[i]] <- foreach(k = 1:times, .combine = "cbind", .packages = c("phyloseq", "permute")) %dopar% {
      species <- vector(length = 30) # create a vector to save diversity
      
      # shuffled sequence e.g. sample_seq = c(4,2,3,1,5) for 5 samples in one site
      # take out samples number 4,2,3,1,5 and add them one by one
      sample_seq <- shuffle(c(1:dim1[1]))
      
      # take out otu_tab for each site
      otu_tab <- matrix(otu_tab, nrow = dim1[1], byrow = TRUE)
      
      # j = 1 as a special case because line 23
      # 'temp <- colSums(otu_tab[c(sample_seq[1:j]),])' 
      # would give error
      species[1] <- sum(otu_tab[sample_seq[1],] > 0)
      
      for (j in 2:30){ 
        # take out samples as the sequence in sample_seq
        temp <- colSums(otu_tab[c(sample_seq[1:j]),])
        # count species
        species[j] <- sum(temp > 0)
      }
      ex <- as.data.frame(cbind(species, "A"=c(1:30)))
      temp <- summary(nls(species~c*A^z,ex,start = list(c=1,z=1)))[["coefficients"]]
      return(c(temp[2,1], temp[2,4]))
    }
    stopCluster(cl)
    gc()
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

write.table(power.z.matrix,"C_power.z.NEON.txt")
write.table(pvalue.matrix,"C_pvalue.NEON.txt")

write.table(power.z.matrix,"C_power.z.NEON.O.txt")
write.table(pvalue.matrix,"C_pvalue.NEON.O.txt")

write.table(power.z.matrix,"C_power.z.NEON.M.txt")
write.table(pvalue.matrix,"C_pvalue.NEON.M.txt")


