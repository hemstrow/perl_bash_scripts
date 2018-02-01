### This excersise has been modified from the original script provided by Matteo Fumagalli

# load all R functions we need
source("~/scripts/sfs_demography_functions.R")

## 1) READ AND STORE OBSERVED GENETIC DATA (2D-SFS)

# the file "Data/polar.brown.sfs" includes the joint (2D) site frequency spectrum (SFS) between polar bears (on the rows) and brown bears (on the columns)
# if you want to see this file type "cat Data/polar.brown.sfs" in your terminal

# read this file and store the the joint SFS of polar VS brown bears into a matrix
COLR.HUMB.2dsf = read.table("~/longfin/pop_gen/demographic_modeling_sfs/COLR_HUMB/COLR.HUMB.2dsfs",header=FALSE)
COLR.HUMB.2dsf = as.numeric(COLR.HUMB.2dsf[1,])
dim(COLR.HUMB.2dsf) = c(17,13)
# this matrix represents the joint (2D) site frequency spectrum (SFS)

#

# QUESTION: what is the sample size for this dataset?

# ANSWER: this is the unfolded spectrum which means that each population has 2n+1 entries in its spectrum, with n being the number of individuals
COLR.nr_ind = (ncol(COLR.HUMB.2dsf)-1)/2
HUMB.nr_ind = (nrow(COLR.HUMB.2dsf)-1)/2
# on the other hand the number of chromosomes can be retrieved as
COLR.nr_chrom = ncol(COLR.HUMB.2dsf)-1
HUMB.nr_chrom = nrow(COLR.HUMB.2dsf)-1

# -----

# QUESTION: how many sites (and SNPs) are there in this 2D-SFS?

# ANSWER: the number of analysed sites is simply the sum of all entries in the SFS
nr_sites = sum(COLR.HUMB.2dsf)
# whereas the number of polymorphic sites is equal to the number of sites minus the count of sites with joint allele frequency (0,0) or (2*polar.nr_ind; 2*brown.nr_ind)
nr_snps = as.numeric(nr_sites-COLR.HUMB.2dsf[1,1]-COLR.HUMB.2dsf[(2*HUMB.nr_ind+1), (2*COLR.nr_ind+1)])

# -----

# for convenience, let's set to NA entries in the matrix where the sites is not a SNP
#polar.brown.sfs[1,1]=NA
#polar.brown.sfs[(2*polar.nr_ind+1), (2*brown.nr_ind+1)]=NA

# plot the spectrum
#plot2DSFS(polar.brown.sfs, xlab="Polar", ylab="Brown", main="2D-SFS")

# this is the unfolded spectrum BUT it has been generated using a reference allele to polarise the spectrum; this means that the polarisation is arbitrary and we need to use the folded spectrum instead
COLR.HUMB.2dsf = fold2DSFS(COLR.HUMB.2dsf)

# plot the folded spectrum
#plot2DSFS(polar.brown.sfs, xlab="Polar", ylab="Brown", main="folded 2D-SFS")

# as an illustration, compute the FST value for the comparison polar vs brown
COLR.HUMB.fst = doFST(COLR.HUMB.2dsf)

# -----

## 2) SIMULATE GENETIC DATA UNDER DIFFERENT VALUES OF THE PARAMETER TO BE ESTIMATED

# define how many simulations we want to perform (ideally a lot)
nr_simul = 1000

# define the prior distribution of our parameter to be estimated (divergence time)
# use a uniform prior bounded at realistic range values
#tdiv_min<-300000 # 300k years ago
#tdiv_max<-700000 # 700k years ago
gen_time = 0.33
mu = 2e-8
ref_theta = 0.004468 ### HUMB ###
ref_Nef = ref_theta/(4*mu*gen_time)
# convert time in coalescent units
# we need generation time and the reference effective population size we used for writing our model in ms
# T(years)=T(coal)*4*Nref
#ref_pop_size <- 68000
min_M12 = 0
max_M12 = 500
min_M21 = 0
max_M21 = 500

# -----

# initialise output files
cat("", file="~/longfin/pop_gen/demographic_modeling_sfs/COLR_HUMB/results/Theta_M_distance_1K.txt")

debug<-FALSE # set this to TRUE if you want to plot intermediate results

# set the directory wher you installed "ms" software
ms_dir = "~/bin/ms" # this is my specific case, yours could be different

# iterate across all simulations
for (i in 1:nr_simul) # run this "for loop" nr_simul times
{

	# print
        cat("\n",i, "/", nr_simul)

	# pick a divergence time randomly from the prior distribution defined above
  #tdiv_random_coal<-runif(1, min=tdiv_min_coal, max=tdiv_max_coal)
  M12_random = runif(1, min = min_M12, max = max_M12)
  M21_random = runif(1, min = min_M21, max = max_M21)
#  m12_random = M12_random/ref_Nef
#  m21_random = M21_random/ref_Nef 

   # convert this value into years-units
	#tdiv_random<-round(tdiv_random_coal*gen_time*4*ref_pop_size)
	# record the sampled value of divergence time in a file
	cat(M12_random, "\t", M21_random, "\t", file="~/longfin/pop_gen/demographic_modeling_sfs/COLR_HUMB/results/Theta_M_distance_1K.txt", append=TRUE)
	cat("\t", M12_random, "\t", M21_random)

	# simulate genetic data of polar and brown bears under this sampled divergence time
	# use "ms" to simulate data

	# initalise ms output file (create an empty file)
	cat("", file="ms.txt")        

	# QUESTION: how many SNPs do we need to simulate?

	# ANSWER: this should match the observed value, stored in "nr_snps" 

	ms.command <- paste(ms_dir, "28", nr_snps, "-s 1 -I 2 16 12 -n 1 1 -n 2 1.8 -m 1 2", M12_random," -m 2 1 ", M21_random," > ms.txt")
	#ms 15 1000 -t 2.0 -I 3 10 4 1 5.0 -m 1 2 10.0 -m 2 1 9.0
	# run ms
	system(ms.command)

	# read ms output file and compute the 2D-SFS

	simulated.2dsfs<-fromMStoSFSwith1site("ms.txt", nr_snps, HUMB.nr_chrom, COLR.nr_chrom)

	# fold it
	simulated.2dsfs<-fold2DSFS(simulated.2dsfs)
	simulated.2dsfs[1,1]=NA

	# if you want to plot it
  #	plot2DSFS(simulated.sfs, xlab="Polar", ylab="Brown", main="Simulated 2D-SFS")
	# and calculate FST
	cat("\t", doFST(simulated.2dsfs))

	# compute the Eucledian distance between observed and simulated SFS (our summary statistics)
	# sqrt((OBS-SIM)^2)
	eucl_dist<-sum(sqrt((COLR.HUMB.2dsf-simulated.2dsfs)^2), na.rm=T)
	# store it in a file
	cat(eucl_dist, "\n", file="~/longfin/pop_gen/demographic_modeling_sfs/COLR_HUMB/results/Theta_M_distance_1K.txt", append=T)

	cat("\t", eucl_dist) # print in std output

}
