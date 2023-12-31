rm(list=ls())
gc()
library(sars)
library(phyloseq)
library(ggplot2)
library(ggtrendline)
library(gridExtra)
library(permute)
library(reshape2)
library(ggplot2)
set.seed(10001)

# 0.function to compare sequence #
longest_common_subsequence <- function(seq1, seq2) {
  m <- length(seq1)
  n <- length(seq2)
  
  # Create a matrix to store the lengths of common subsequences
  lcs_lengths <- matrix(0, nrow = m + 1, ncol = n + 1)
  
  # Iterate through the sequences to find the longest common subsequence
  for (i in 1:m) {
    for (j in 1:n) {
      if (seq1[i] == seq2[j]) {
        lcs_lengths[i + 1, j + 1] <- lcs_lengths[i, j] + 1
      } else {
        lcs_lengths[i + 1, j + 1] <- max(lcs_lengths[i + 1, j], lcs_lengths[i, j + 1])
      }
    }
  }
  
  # Backtrack to reconstruct the longest common subsequence
  lcs <- character(0)
  i <- m
  j <- n
  while (i > 0 && j > 0) {
    if (seq1[i] == seq2[j]) {
      lcs <- c(seq1[i], lcs)
      i <- i - 1
      j <- j - 1
    } else if (lcs_lengths[i + 1, j] > lcs_lengths[i, j + 1]) {
      j <- j - 1
    } else {
      i <- i - 1
    }
  }
  
  return(lcs)
}

neon_dob <- readRDS("../data/phylo_V3.1.RDS")
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
neon <- subset_samples(neon, !is.na(horizon))
rm(neon_dob)

## 1.neon dataset ##
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
}

z_df <- cbind(power.increase.z[,1], power.decrease.z[,1], power.random.z[,1])
colnames(z_df) <- c("increase","decrease","random")
z_df_long <- melt(z_df)
plot(x = z_df_long$Var2, y = z_df_long$value, xlab = "the order of pooling",
     ylab = "z value", 
     main = "Box plot of the z values calculated with different pooling method")
write.table(power.increase.z,"NEONPowerInZ.txt")
write.table(power.increase.z,"NEONPowerDeZ.txt")

length(longest_common_subsequence(names(sort(power.increase.z[,1])),names(sort(power.decrease.z[,1])))) # 
length(longest_common_subsequence(names(sort(power.increase.z[,1])),names(sort(power.random.z[,1])))) # 
length(longest_common_subsequence(names(sort(power.decrease.z[,1])),names(sort(power.random.z[,1])))) # 

## 2   ##
## dob ##
dob <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="DoB")
dob <- subset_samples(dob, !is.na(Site))
rm(neon_dob)

d1 <- sample_data(dob)

a1 <- unique(d1$Site)# be noted that the powering site orders may not be identical to that of generated by the 'merge_samples' function you used at line 128
power.increase.c <- power.increase.z <- power.decrease.c <- power.decrease.z <- 
  power.random.c <- power.random.z <- 
  matrix(nrow = length(a1), ncol = 4,dimnames = list(a1, c("Estimate", "Std. Error", "t value", "Pr(>|t|)" )))

for (i in 1:length(a1)){
  cat('\r',paste(paste0(rep("*", round(i/ 1, 0)), collapse = ''), i, "/", length(a1), collapse = ''))# informs the processing
  dob_sub <- subset_samples(dob, Site==a1[i])
  dim1 <- dim(otu_table(dob_sub)) # the number of samples in one site
  otu_tab <- otu_table(dob_sub)
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
ggplot(data = z_df_long)+
  geom_density(aes(x = value, color = Var2, alpha = 0.2, fill = Var2))

write.table(power.increase.z,"DOBPowerInZ.txt")
write.table(power.increase.z,"DOBPowerDeZ.txt")

length(longest_common_subsequence(names(sort(power.increase.z[,1])),names(sort(power.decrease.z[,1])))) # 18 out of all 68
length(longest_common_subsequence(names(sort(power.increase.z[,1])),names(sort(power.random.z[,1])))) # 23
length(longest_common_subsequence(names(sort(power.decrease.z[,1])),names(sort(power.random.z[,1])))) # 22
cor(power.increase.z[,1], power.decrease.z[,1], method = "kendall") # 0.55
cor(power.increase.z[,1], power.random.z[,1], method = "kendall") # 0.64
cor(power.decrease.z[,1], power.random.z[,1], method = "kendall") # 0.59
