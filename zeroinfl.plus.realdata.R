# Here we use a real throat microbiome data set (Charlson, 2010, PLoS ONE) to illustrate the proposed method.
# We aim to identify bacterial OTUs that differentiates between smokers and nonsmokers.
# The basic process for real data is as follows:
# 1. Estimate the libary size using GMPR
# 2. Replace outliers using winsorization (97% quantile)
# 3  [Filtering / model choice - will not investigate here]
# 4. Run the omnibus test

source('zeroinfl.plus.github.R')

require(GUniFrac)
require(matrixStats)

data(throat.otu.tab)
data(throat.meta)

otu.tab <- t(throat.otu.tab)

# Step 1: estimate the library size
size.factor <- GMPR(otu.tab)
hist(size.factor)

# Step 2: winsorization at 97% quantile
otu.tab.p <- t(t(otu.tab) / size.factor)
otu.tab.p <- apply(otu.tab.p, 1, function(x) {
			cutoff <- quantile(x, 0.97)
			x[x >= cutoff] <- cutoff
			x
		}
)

otu.tab.win <- t(round(otu.tab.p * size.factor))
range(rowMeans(otu.tab.win != 0))
hist(rowMeans(otu.tab.win != 0))

# Step 3-4: For demonstration purpose, we just run on a subset of OTUs with many 0's,
# In practice, we may need to make decision what OTUs will be tested and which model
# will use (goodness-of-test). 

otu.tab.win.sel <- otu.tab.win[rowMeans(otu.tab.win != 0) > 0.20 & 
				rowMeans(otu.tab.win != 0) < 0.80, ]
res <- sapply(1:nrow(otu.tab.win.sel),  function (i){
			if (i %% 10 == 0) cat('.')
			y <- otu.tab.win.sel[i, ]
			error <- try(
					zinb.obj <- zinb.lrt(
							formula.H1 = y ~ SmokingStatus + offset(log(size.factor)) | SmokingStatus | SmokingStatus, 
							formula.H0 = y ~ 1 + offset(log(size.factor)) | 1 | 1, data = throat.meta,
							control = zinb.control(trace=FALSE)
					)
			)
			if (class(error) == 'try-error') {
				return(NA)
			} else {
				zinb.obj$p.value
			}
			
		})
hist(res)
sort(res)
