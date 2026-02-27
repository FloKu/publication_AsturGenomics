# *************************************************************
# Post hoc analyses of genomic accipiter data
# *************************************************************
# This script reads vcf file from STACKS output, identifies missing percentages, applies data checks, filters samples and calculates DAPC.
# author: Florian Kunz
# last changed: 18.02.26

{library("dartRverse")
  library("adegenet")
  library("vcfR")
  library("ggplot2")
  library("ggsci")
  library("mmod")
}
  
path <- ""


# 1) read in data
# **************************

# Full data: 64 individuals, 40116 SNPs, 28.62% missing data
gl_full <- dartR::gl.read.vcf(paste0(path, "analysis/final_data_vcf/FINAL_63inds.vcf"))
gl_full$other$ind.metrics <- read.csv(paste0(path, "analysis/population.txt"))

pop(gl_full) <- gl_full$other$ind.metrics$pop
plyr::count(as.vector(gl_full$pop))

# 2) check SNP data
# **************************
#----
dartR::gl.compliance.check(gl_full)
gl_full
gl_full$ind.names
gl_full$pop

# inds per pop
gi_full <- dartR::gl2gi(gl_full)
summary(gi_full)$n.by.pop

dartR::gl.report.monomorphs(gl_full)
# no monomorphic loci

dartR::gl.report.callrate(gl_full, by_pop=F)
# about 75% callrate on average - most SNPs pretty good, but a some have very low call rate

dartR::gl.report.callrate(gl_full, by_pop=T)
# looking at per pop, I believe these SNPs with the low call rates are the ones which are unique to the populations, so they are important for differentiation

# overall individuals
dartR::gl.report.callrate(gl_full, by_pop=F, method="ind")

# 2.1.) check callrates per individual
# **************************
# plot
glPlot(gl_full)
dartR::gl.smearplot(gl_full)

# missing data per individual
missing_percentage <- rowSums(is.na(as.matrix(gl_full)))
total_loci <- nLoc(gl_full)
missing <- as.data.frame(missing_percentage) %>% 
  mutate(missing_percentage=(missing_percentage / total_loci) * 100) %>% 
  mutate(pop = gl_full$pop)
missing

# subset of NAM, HEN and MEL
subset(missing, pop == "NAM" | pop == "HEN" | pop == "MEL")

# excluding all above a certain percentage
subset(missing, missing_percentage > 66)

# create subset
filtered <- rownames(subset(missing, missing_percentage < 66))
individuals_to_keep <- indNames(gl_full) %in% filtered # subset gnlight needs atomic vector
gl_filtered <- gl_full[individuals_to_keep, ]

# filtering using dartR, for checks
gl_filtered_dR <- dartR::gl.filter.callrate(gl_full, method="ind", threshold=0.34, recalc=T)

# Check loci again in the filtered data
dartR::gl.report.monomorphs(gl_full) # monomorphic loci
dartR::gl.filter.allna(gl_filtered, by.pop = FALSE, recalc = T, verbose = 5)

# 2.2.) check filtered data
# **************************
dartR::gl.compliance.check(gl_filtered)
gl_filtered$pop

dartR::gl.report.monomorphs(gl_filtered)
# 213 monomorphic loci

gl_filtered <- dartR::gl.filter.monomorphs(gl_filtered)
# 39903 SNPs remaining

dartR::gl.report.callrate(gl_filtered, by_pop=F)
# about 83% callrate on average - most SNPs pretty good, but a some have very low call rate

# overall individuals
dartR::gl.report.callrate(gl_filtered, by_pop=F, method="ind")

# check callrates per individual
glPlot(gl_filtered)
dartR::gl.smearplot(gl_filtered)

# missing data per individual
missing_percentage <- rowSums(is.na(as.matrix(gl_filtered)))
total_loci <- nLoc(gl_filtered)
missing <- as.data.frame(missing_percentage) %>% 
  mutate(missing_percentage=(missing_percentage / total_loci) * 100) %>% 
  mutate(pop = gl_filtered$pop)
missing

# last filter step
dartR::gl.filter.allna(gl_filtered, by.pop = FALSE, recalc = T, verbose = 5)

# save cleaned version
dartR::gl.save(gl_filtered, paste0(path, "analysis/final_data_vcf/FINAL_52inds.rds"), verbose=5)
dartR::gl2snapp(gl_filtered, outfile = "FINAL_52inds_dartR.nexus", outpath = paste0(path, "analysis/final_data_nexus/"), verbose = NULL) # sadly, this nexus is not working in IQ-TREE as the datatype cannot be intergerdata

dartR::gl2vcf(gl_full, plink_path="C:/Users/floriankunz/plink_win64_20250819", outfile="test_full", outpath=paste0(getwd(), "/final_data_vcf/dartR")) # not sure if this works, it looks different
#----

# Read in filtered data: 53 individuals, 39903 SNPs, 16.89% missing data
gl_filtered <- dartR::gl.load(paste0(path, "analysis/final_data_vcf/FINAL_52inds.rds"))


# 3) DAPC
# **************************
#----
gl_filtered_subset <- gl_filtered[pop(gl_filtered) != "OUT", ]
gl_filtered_subset@pop

# change names of subpops
pop(gl_filtered_subset) <- ifelse(pop(gl_filtered_subset) == "EUR", "A. gentilis", ifelse(pop(gl_filtered_subset) == "HEN", "A. henstii", ifelse(pop(gl_filtered_subset) == "MEL", "A. melanoleucus", "A. atricapillus")))

# cross-validation
# split
mat <- tab(gl_filtered_subset, NA.method = "mean")
grp <- pop(gl_filtered_subset)

# evaluate
xval <- xvalDapc(mat, grp, n.pca.max = 100, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval[2:7]

# The model with the lowest MSE is recommended in the tutorial
xval[6]

# extracting the best DAPC as chosen by lowest MSE
dapc <- xval$DAPC

# Explore Results
dapc
summary(dapc)

dapc$n.pca # retained PCA axes
dapc$n.da # retained DA axes
dapc$var # proportion of variance covered by PCA principal components

# How many percet are correclty assgined per biological group?
summary(dapc)$assign.per.pop

assignplot(dapc)
compoplot(dapc)

# Plot DAPC
col <- c("cornflowerblue", "darkorange", "turquoise3", "springgreen3")

png("C:/Users/floriankunz/2020_projects/2020_Habicht_FollowUp/manuscript/graphs/DAPC.png", 16000, 8000, res=1000, type='cairo-png', antialias=c("subpixel")) 
par(mfrow = c(1, 2), cex.axis = 1.5, cex.lab = 1.5)
scatter(dapc, xax=1, yax=1, col=col, legend=F, posi.leg="topright", scree.pca = F, posi.pca="topright", cleg=1.5)
labs <- levels(dapc$grp)              # group labels used by scatter.dapc
legend("topright", legend = labs, col = col, pch = 20, pt.cex = 2.5, cex = 1.5,
       text.font = 3,  # 1=plain, 2=bold, 3=italic, 4=bold-italic
       bty = "o",       # draw box ("n" = none)
       box.lwd = 1.5,   # box line width
       box.lty = 1,     # box line type
       box.col = "black",
       bg = "white")
scatter(dapc, xax=2, yax=2, col=col, legend=F, posi.leg="topright", scree.pca = T, posi.pca="topleft", ratio.pca=0.3)
par(mfrow = c(1, 1))
dev.off()

#----


# 4) G''ST and DJost
# **************************
#----
# Package calculates G''ST (Hedrick 2011), not G'ST (Hedrick 2005) as indicated in the documentation, confirmed with the package author by Florian Kunz

# Package diveRsity does bias correction on CI automatically, method explained in help manual - mmod does not!
# So recreated that method using mmod:
# 1) Pairwise indices (point estimates) are calculated
# 2) 1000 bootstraped resamples are created with mmod - this function creates alist of genind objects that are boostrapped samples from the original genind, and the boostrapping is thereby accounting for sample size of the populations
# 3) Then a loop is initiated, which calculates and saves the pairwise genetic distance of each genind object within the list of boostrapped samples
# 4) Table is then analysed to calculate mean
# 5) Bias corrected CI are calulated by substracting the difference between mean and point estimate from the CI interval, basically re-centering it
#   See guide_genetic_distance_metrics.docx for a detailed explanation of the bias correction.

gl_filtered_subset <- gl_filtered[pop(gl_filtered) != "OUT", ]
gl_filtered_subset@pop

# change names of subpops
pop(gl_filtered_subset) <- ifelse(pop(gl_filtered_subset) == "EUR", "A. gentilis", ifelse(pop(gl_filtered_subset) == "HEN", "A. henstii", ifelse(pop(gl_filtered_subset) == "MEL", "A. melanoleucus", "A. atricapillus")))

# create genind
genind.data <- dartR::gl2gi(gl_filtered_subset)
genind.data
genind.data@pop

pwD <- pairwise_D(genind.data, linearized = FALSE, hsht_mean = "arithmetic")
pwGst <- pairwise_Gst_Hedrick(genind.data, linearized = FALSE)

rep <- 1000 # Number of bootstrap replicates
boots <- chao_bootstrap(genind.data, nreps = rep) # Generating new genind objects with bootstrapping

pairwise_D_all <- matrix(0, ncol=6, nrow=rep) # Generates an emtpy matrix that will be filed with resultsdata
colnames(pairwise_D_all) <- c("A. atricapillus/A. gentilis", "A. atricapillus/A. henstii", "A. atricapillus/A. melanoleucus", "A. gentilis/A. henstii", "A. gentilis/A. melanoleucus", "A. henstii/A. melanoleucus")
for (i in 1:rep) {
  print(i)
  pw<-pairwise_D(boots$BS[[i]], linearized = FALSE, hsht_mean = "arithmetic") 
  pairwise_D_all[i,]  <- pw
}

pwD_sum <- matrix(0, ncol=7, nrow=6, dimnames = list(c("A. atricapillus/A. gentilis", "A. atricapillus/A. henstii", "A. atricapillus/A. melanoleucus", "A. gentilis/A. henstii", "A. gentilis/A. melanoleucus", "A. henstii/A. melanoleucus"), c("comparison", "point estimate", "mean", "sd", "CI_lower_bc", "CI_upper_bc", "scenario")))
pwD_sum[,2] <- pwD
pwD_sum[,3] <- apply(pairwise_D_all, 2, mean)
pwD_sum[,4] <- apply(pairwise_D_all, 2, sd)
biasD <- pwD_sum[,3]-pwD_sum[,2]
pwD_sum[,5] <- (qnorm(0.025, pwD_sum[,3], pwD_sum[,4])-biasD)
pwD_sum[,6] <- (qnorm(0.975, pwD_sum[,3], pwD_sum[,4])-biasD)
pwD_sum[,1] <- c("A. atricapillus/A. gentilis", "A. atricapillus/A. henstii", "A. atricapillus/A. melanoleucus", "A. gentilis/A. henstii", "A. gentilis/A. melanoleucus", "A. henstii/A. melanoleucus") # It is important to insert text as a last
pwD_sum[,7] <- c("pwD")

write.csv(pwD_sum, "C:/Users/floriankunz/2020_projects/2020_Habicht_FollowUp/analysis/pairwiseD_mmod_results_1000bt.csv")

pairwise_Gst_all <- matrix(0, ncol=6, nrow=rep) # Generates an emtpy matrix that will be filed with results
colnames(pairwise_Gst_all) <- c("A. atricapillus/A. gentilis", "A. atricapillus/A. henstii", "A. atricapillus/A. melanoleucus", "A. gentilis/A. henstii", "A. gentilis/A. melanoleucus", "A. henstii/A. melanoleucus")
for (i in 1:rep) {
  print(i)
  pw<-pairwise_Gst_Hedrick(boots$BS[[i]], linearized = FALSE) 
  pairwise_Gst_all[i,]  <- pw
}

pwGst_sum <- matrix(0, ncol=7, nrow=6, dimnames = list(c("A. atricapillus/A. gentilis", "A. atricapillus/A. henstii", "A. atricapillus/A. melanoleucus", "A. gentilis/A. henstii", "A. gentilis/A. melanoleucus", "A. henstii/A. melanoleucus"), c("comparison", "point estimate", "mean", "sd", "CI_lower_bc", "CI_upper_bc", "scenario")))
pwGst_sum[,2] <- pwGst
pwGst_sum[,3] <- apply(pairwise_Gst_all, 2, mean)
pwGst_sum[,4] <- apply(pairwise_Gst_all, 2, sd)
biasGst <- pwGst_sum[,3]-pwGst_sum[,2]
pwGst_sum[,5] <- (qnorm(0.025, pwGst_sum[,3], pwGst_sum[,4])-biasGst)
pwGst_sum[,6] <- (qnorm(0.975, pwGst_sum[,3], pwGst_sum[,4])-biasGst)
pwGst_sum[,1] <- c("A. atricapillus/A. gentilis", "A. atricapillus/A. henstii", "A. atricapillus/A. melanoleucus", "A. gentilis/A. henstii", "A. gentilis/A. melanoleucus", "A. henstii/A. melanoleucus") # It is important to insert text as a last
pwGst_sum[,7] <- c("pwG2st")

write.csv(pwGst_sum, "C:/Users/floriankunz/2020_projects/2020_Habicht_FollowUp/analysis/pairwiseG2st_mmod_results_1000bt.csv")
#----