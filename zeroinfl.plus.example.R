# Title: Artifical simulations to illustrate our omnibus test for differential abundance
# analysis of microbiome-seq data. Only one taxon are simulated for illustration purpose.
# Authors: Jun Chen (chen.jun2@mayo.edu)
# Date: 2017/02/07

source('~/Dropbox/Workspace/MayoClinic/Methodology/2014_06_01_Omnibus_Test/zeroinfl.plus.github.R')

ns <- 100      # Number of sample
zp <- 0.4      # Zero propability
prp <- 0.02    # OTU relative abundance
iter <- 50     # Iteration number
grp <- as.numeric(gl(2, ns/2)) # Variable of interest
grp0 <- grp
#####################################
# Simulation 1: Type I error control
#####################################
pv.lrt1 <- pv.lrt2 <- NULL
for (i in 1:iter) {
	if (i %% 10 == 0) cat(".")
	# Generate library size
	size.factor <- rnbinom(ns, mu = 2000, size = 2.5)
	# Generate structural zeros
	zind <- rbinom(ns, 1, zp)
	y <- numeric(ns)
	# Generate counts for nonzero part
	y[zind != 0] <-  rnbinom(sum(zind != 0), mu = size.factor[zind != 0] * prp, size = 1.0)
	data <- data.frame(y = y, grp = grp)
	
	# ominbus test for prevalence/abundance/dispersion
	zinb.obj1 <- zinb.lrt(formula.H1 = y ~ grp + offset(log(size.factor)) | grp | grp, 
			formula.H0 = y ~ 1 + offset(log(size.factor)) | 1 | 1, data=data,
			control = zinb.control(trace=TRUE))
	
	# test for differential prevalence/abundance, allowing dispersion depending on the group
	zinb.obj2 <- zinb.lrt(formula.H1=y ~ grp + offset(log(size.factor)) | grp | grp, 
			formula.H0 = y ~ 1 + offset(log(size.factor)) | 1 | grp, data=data,
			control = zinb.control(trace=FALSE))
	
	pv.lrt1 <- c(pv.lrt1, zinb.obj1$p.value)
	pv.lrt2 <- c(pv.lrt2, zinb.obj2$p.value)
}
par(mfrow=c(1, 2))
hist(pv.lrt1)
hist(pv.lrt2)

###########################################
# Simulation 1: Power - prevalence change
###########################################
grp <- as.numeric(grp0)
grp <- grp - mean(grp)
pv.lrt1 <- pv.lrt2 <- NULL
for (i in 1:iter) {
	if (i %% 10 == 0) cat(".")
	size.factor <- rnbinom(ns, mu = 2000, size = 2.5)

	zind <- rbinom(ns, 1, exp(grp) / (1 + exp(grp)))
	
	y <- numeric(ns)
	y[zind != 0] <-  rnbinom(sum(zind != 0), mu = size.factor[zind != 0] * prp, size = 1.0)
	data <- data.frame(y = y, grp = grp)
	# ominbus test for prevalence/abundance/dispersion
	zinb.obj1 <- zinb.lrt(formula.H1 = y ~ grp + offset(log(size.factor)) | grp | grp, 
			formula.H0 = y ~ 1 + offset(log(size.factor)) | 1 | 1, data=data,
			control = zinb.control(trace=FALSE))
	
	# test for differential prevalence/abundance, allowing dispersion depending on the group
	zinb.obj2 <- zinb.lrt(formula.H1=y ~ grp + offset(log(size.factor)) | grp | grp, 
			formula.H0 = y ~ 1 + offset(log(size.factor)) | 1 | grp, data=data,
			control = zinb.control(trace=FALSE))
	pv.lrt1 <- c(pv.lrt1, zinb.obj1$p.value)
	pv.lrt2 <- c(pv.lrt2, zinb.obj2$p.value)
}
par(mfrow=c(1, 2))
hist(pv.lrt1)
hist(pv.lrt2)

###########################################
# Simulation 2: Power - dispersion change
###########################################
grp <- as.numeric(grp0)
grp <- grp - mean(grp)
pv.lrt1 <- pv.lrt2 <- NULL
for (i in 1:iter) {
	if (i %% 10 == 0) cat(".")
	size.factor <- rnbinom(ns, mu = 2000, size = 2.5)
	zind <- rbinom(ns, 1, zp)
	y <- numeric(ns)
	y[zind != 0] <-  rnbinom(sum(zind != 0), mu = size.factor[zind != 0] * prp, size = grp[zind != 0] + 1.0)
	data <- data.frame(y = y, grp = grp)
	# ominbus test for prevalence/abundance/dispersion
	zinb.obj1 <- zinb.lrt(formula.H1 = y ~ grp + offset(log(size.factor)) | grp | grp, 
			formula.H0 = y ~ 1 + offset(log(size.factor)) | 1 | 1, data=data,
			control = zinb.control(trace=FALSE))
	
	# test for differential prevalence/abundance, allowing dispersion depending on the group
	zinb.obj2 <- zinb.lrt(formula.H1=y ~ grp + offset(log(size.factor)) | grp | grp, 
			formula.H0 = y ~ 1 + offset(log(size.factor)) | 1 | grp, data=data,
			control = zinb.control(trace=FALSE))
	pv.lrt1 <- c(pv.lrt1, zinb.obj1$p.value)
	pv.lrt2 <- c(pv.lrt2, zinb.obj2$p.value)
}
par(mfrow=c(1, 2))
hist(pv.lrt1)
hist(pv.lrt2)

