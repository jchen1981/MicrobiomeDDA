# Here we use a real throat microbiome data set (Charlson, 2010, PLoS ONE) to illustrate the proposed method.
# We aim to identify bacterial OTUs that differentiates between smokers and nonsmokers.
# The basic process for real data is as follows:
# 1. Estimate the libary size using GMPR
# 2. Replace outliers using winsorization (97% quantile)
# 3  Filtering
# 4. Run the omnibus test

#setwd('~/Dropbox/Workspace/MayoClinic/Methodology/2014_06_01_Omnibus_Test/')
source('zeroinfl.plus.github.R')
source('zeroinfl.plus.daa.R')

require(GUniFrac)

data(throat.otu.tab)
data(throat.meta)

otu.tab <- t(throat.otu.tab)

###############################################
# Ominbus test
meta.dat <- throat.meta
obj <- ZISeq(otu.tab, meta.dat, 
		size.factor = NULL,                                  # Normalization (GMPR)
		winsor = TRUE, winsor.qt = 0.97,                     # Winsorization
		grp.name = 'SmokingStatus', adj.name = NULL, 
		method = 'omnibus',
		filter = TRUE, prev.filter = 0.1, ct.cutoff = 10     # Filter
	) 
obj$result


# Prevalence/abundance test - dispersion still depends on covariate
meta.dat <- throat.meta
obj <- ZISeq(otu.tab, meta.dat, 
		size.factor = NULL,                                  # Normalization (GMPR)
		winsor = TRUE, winsor.qt = 0.97,                     # Winsorization
		grp.name = 'SmokingStatus', adj.name = NULL, 
		method = 'prev.abund1',
		filter = TRUE, prev.filter = 0.1, ct.cutoff = 10     # Filter
) 
obj$result

# Prevalence/abundance test - common dispersion
meta.dat <- throat.meta
obj <- ZISeq(otu.tab, meta.dat, 
		size.factor = NULL,                                  # Normalization (GMPR)
		winsor = TRUE, winsor.qt = 0.97,                     # Winsorization
		grp.name = 'SmokingStatus', adj.name = NULL, 
		method = 'prev.abund2',
		filter = TRUE, prev.filter = 0.1, ct.cutoff = 10     # Filter
) 
obj$result

###############################################
