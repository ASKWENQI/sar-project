rm(list=ls())
neon_dob <- readRDS("/data/ZHULAB/soil/NEON_DOB/phylo_V3.1.RDS")

# We can't use any sites that lack coordinates
neon_dob <- subset_samples(neon_dob, !is.na(lon) & !is.na(lat))

# No need to agglomerate to a standard taxonomic level; all OTUs
# are already defined at 97% sequence similarity.

# Remove zero-abundance taxa and samples
neon_dob <- subset_taxa(neon_dob, taxa_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, sample_sums(neon_dob) > 0)
neon_dob <- subset_samples(neon_dob, !is.na(Site))
neon <- subset_samples(neon_dob, get_variable(neon_dob, "Project")=="NEON")
neon <- subset_samples(neon, !is.na(horizon))

# Aggregate sites within 10-minute of a degree; this is the resolution of the
# climate data
xgrid <- 1/6*seq(1,600)-160 # Based on longitudinal extent of sample data
ygrid <- 1/6*seq(1,360)+15 # Based on latitudinal extent of sample data
sample_data(neon) %>%
  as("data.frame") %>%
  mutate(coordsSite = paste(round(xgrid[findInterval(lon, xgrid)], 2),
                            round(ygrid[findInterval(lat, ygrid)], 2), sep="_")) ->
  sample_data(neon)

neon_agg <- merge_samples(neon, group="coordsSite")
