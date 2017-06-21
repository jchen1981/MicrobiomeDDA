# Title: An omnibus test for different distribution analysis of zeroinflated seq data
# Version: 0.0.1
# Authors: Jun Chen (chen.jun2@mayo.edu)
# Description: A user-friendly interface to perform the omnibus test. The input
# are simply the OTU table (row - OTUs, column - samples) and the meta data file. 
# Note: Due to many parameters and asymptotics invovled, it is not recommended
# to run on a small smaple size (e.g. < 50).
# Date: 2017/06/12

# New: 2017_06_20

ZISeq <- function (otu.tab, meta.dat, 
		size.factor = NULL,                                  # Normalization (GMPR)
		winsor = TRUE, winsor.qt = 0.97,                     # Winsorization
		grp.name, adj.name = NULL, 
		method = c('omnibus', 'prev.abund1', 'prev.abund2', 'general'), 
		formula.H1 = NULL, formula.H0 = NULL, 
		filter = TRUE, prev.filter = 0.1, ct.cutoff = 10,    # Filter
		fail.method = 'permutation')    
{                       
	
	# Args:
	#		otu.tab: OTU count table, row - OTUs, column - samples
	#       meta.dat: a data.frame recording the information of the samples.  Make sure the 
	#			   variables are cast into the appropriate data types.
	#		size.factor: a vector of the size factors to address variable sequencing depth. 
	#              The default is NULL and GMPR normalization will be used.
	#		winsor: a logical value indicating whether winsorization should be performed to
	#			   replace outliers.  The default is TRUE.
	#       winsor.qt: the winsorization quantile, above which the counts will be replaced
	#       grp.name: the name of the co-variate of interest. It could be numeric or categorical
	#       adj.name: a vector of the names of the covariates to be adjusted.  
	#       method:  the method to be used. The default is the 'omnibus' test. A joint test of
	#			  prevalence and abundance parameter ('prev.abund1') is also given, where the
	#			  dispersion is still allowed to depend on the covariates. 'prev.abund2' assumes
	#             a common dispersion parameter (i.e. does not depend on the covariates). The 'general'
	#             method allows the user to specify its H0 and H1 hypotheses (see below).            
	#       formula.H1, formula.H0: characters, the formulae for the H0 and H1 hypotheses. These 
	#		are for advanced use when method = 'general'. Example:
	#			  formula.H1 = "y ~ Smoking + Batch + offset(log(size.factor)) | Batch + Smoking | Batch", 
	#			  formula.H0 = "y ~ 1 + Batch + offset(log(size.factor)) | Batch | Batch"
	#             which tests the smoking effects on prevalence/abundance while adjusting batch effects. The
	#             batch variable usually affects both the scale and dispersion parameter so we let it 
	#             depend on all the parameters.
	#   	filter: a logical value indicating whether filtering should be performed
	#       prev.filter, ct.cutoff: filter criteria,  taxa with at least 'ct.cuoff' counts
	#              in at least 'prev.filter' samples are tested.
	#       fail.method: the method used for failed taxa (due to very low prevalence/abundance, or
	#              very high prevalence/abundance)
	#
	# Returns:
	# 	result: a data.frame containing the p-values, test statistic, degree of freedoms,
	#        baseline prevalence, abundance, dispersion, estimates of the parameters, their 
	#        standard errors and the method used.
	#	size.factor: the size.factor used for normalization
	#   test.ind: a vector of logical values indicating which taxa are tested
        #
	#
	if (sum(colnames(otu.tab) != rownames(meta.dat)) != 0) {
		stop("The sample names for the otu table and the meta data table do not agree!\n")
	}
	
	method <- match.arg(method)
	
# Step 1: estimate the library size
	if(is.null(size.factor)) {
		cat('Start GMPR normalization ...\n')
		size.factor <- GMPR(otu.tab)
		size.factor <- size.factor * median(colSums(otu.tab))
	} else {
		if (length(size.factor) != ncol(otu.tab)) {
			stop("The length of size factor differs from the sample size!\n")
		}
	}
	
	
# Step 2: winsorization at 97% quantile
	if (winsor == TRUE) {
		cat('Start Winsorization ...\n')
		otu.tab.p <- t(t(otu.tab) / size.factor)
		otu.tab.p <- apply(otu.tab.p, 1, function(x) {
					cutoff <- quantile(x, winsor.qt)
					x[x >= cutoff] <- cutoff
					x
				}
		)
		otu.tab.win <- t(round(otu.tab.p * size.factor))
	} else {
		otu.tab.win <- otu.tab
	}
	
	
# Step 3: Filtering
	if (filter == TRUE) {
		cat('Perform filtering ...\n')
		test.ind <- rowMeans(otu.tab.win >= ct.cutoff) > prev.filter
		otu.tab.win.sel <- otu.tab.win[test.ind, ]
		cat('--A total of ', nrow(otu.tab.win.sel), ' taxa will be tested with a sample size of', nrow(meta.dat), '!\n')
	} else {
		cat('No filtering will be performed.\n')
		test.ind <- rep(TRUE, nrow(otu.tab.win))
		otu.tab.win.sel <- otu.tab.win
		cat('--A total of ', nrow(otu.tab.win.sel), ' taxa will be tested with a sample size of', nrow(meta.dat), '!\n')
	}		
	
	
	
# Step 4: Create formula and perform test
	if (method == 'omnibus') {
		cat('--Omnibus test is selected!\n--Dispersion is treated as a parameter of interest!\n')
		if (!is.null(adj.name)) {
			formula.H1 <- (paste0('y ~ ', grp.name, '+', paste(adj.name, collapse = ' + '), ' + offset(log(size.factor)) | ',
							grp.name, ' + ', paste(adj.name, collapse = ' + '), ' | ', 
							grp.name, ' + ', paste(adj.name, collapse = ' + ')))
			formula.H0 <- (paste0('y ~ ', paste(adj.name, collapse = ' + '), ' + offset(log(size.factor)) | ',
							paste(adj.name, collapse = ' + '), ' | ', 
							paste(adj.name, collapse = ' + ')))
		} else {
			formula.H1 <- (paste0('y ~ ', grp.name, ' + offset(log(size.factor)) | ',
							grp.name,  ' | ', 
							grp.name))
			formula.H0 <- ('y ~  1 + offset(log(size.factor)) | 1 | 1 ')
		}
	}
	
	if (method == 'prev.abund1') {
		cat('--Prev.abund test is selected!\n--Dispersion is treated as a nuisance parameter but it is allowed to vary with covariates!\n')
		if (!is.null(adj.name)) {
			formula.H1 <- (paste0('y ~ ', grp.name, '+', paste(adj.name, collapse = ' + '), ' + offset(log(size.factor)) | ',
							grp.name, ' + ', paste(adj.name, collapse = ' + '), ' | ', 
							grp.name, ' + ', paste(adj.name, collapse = ' + ')))
			formula.H0 <- (paste0('y ~ ', paste(adj.name, collapse = ' + '), ' + offset(log(size.factor)) | ',
							paste(adj.name, collapse = ' + '), ' | ', 
							grp.name, ' + ',paste(adj.name, collapse = ' + ')))
		} else {
			formula.H1 <- (paste0('y ~ ', grp.name, ' + offset(log(size.factor)) | ',
							grp.name,  ' | ', 
							grp.name))
			formula.H0 <- (paste0('y ~  1 + offset(log(size.factor)) | 1 | ', 
							grp.name))
		}
	}
	
	if (method == 'prev.abund2') {
		cat('--Prev.abund test is selected!\n--Dispersion is treated as a nuisance parameter and is common for all samples!\n')
		if (!is.null(adj.name)) {
			formula.H1 <- (paste0('y ~ ', grp.name, '+', paste(adj.name, collapse = ' + '), ' + offset(log(size.factor)) | ',
								grp.name, ' + ', paste(adj.name, collapse = ' + '), ' | ', 
								'1'))
			formula.H0 <- (paste0('y ~ ', paste(adj.name, collapse = ' + '), ' + offset(log(size.factor)) | ',
								paste(adj.name, collapse = ' + '), ' | ', 
								'1'))
		} else {
			formula.H1 <- (paste0('y ~ ', grp.name, ' + offset(log(size.factor)) | ',
								grp.name,  ' | ', 
								'1'))
			formula.H0 <- (paste0('y ~  1 + offset(log(size.factor)) | 1 | ', 
								'1'))
		}
	}
	
	if (method == 'general') {
		if (is.null(formula.H1) | is.null(formula.H0)) {
			stop("Method 'general' requires specification of 'formula.H1' and 'formula.H0'!\n")
		}
	}
	
	cat('Start testing ...\n')
	unit <- ceiling(nrow(otu.tab.win.sel) / 10)
	col.names <- NULL
	res <- lapply(1:nrow(otu.tab.win.sel),  function (i){
				if (i %% unit == 0) cat((i %/% unit * 10), '%\n')
				y <- otu.tab.win.sel[i, ]
				error <- try(
						zinb.obj <- zinb.lrt(
								formula.H1 = as.formula(formula.H1), formula.H0 = as.formula(formula.H0), data = meta.dat,
								control = zinb.control(trace=FALSE)
						), silent = TRUE
				)
				
				if (class(error) == 'try-error') {
					return(NA)
				} else {
					if (!is.na(zinb.obj$p.value)) {
						
						tab <- t(zinb.obj$coef[, 1:2])
						tab.names <- as.vector(t(outer(colnames(tab), c('.est', '.se'), paste0)))
						tab.names <- gsub('count.', 'count.LFC.', tab.names)
						tab.names <- gsub('zero.', 'zero.LOD.', tab.names)
						tab.names <- gsub('dispersion.', 'dispersion.LFC.', tab.names)
						
						tab <- as.vector(tab)
						names(tab) <- tab.names
						
						baseline <- sapply(zinb.obj$mod.H1$coefficients, function(x) x[1])
						names(baseline) <- gsub('\\(Intercept\\)', 'baseline', names(baseline))
						baseline[1] <- exp(baseline[1])
						baseline[2] <- exp(baseline[2]) / (1 + exp(baseline[2]))
						baseline[3] <- exp(baseline[3])
						
						ret <- c(p.value=zinb.obj$p.value, chi.stat=zinb.obj$chi.stat, df=zinb.obj$dof, baseline, tab)	
						col.names <<- names(ret)
						return(ret)
					} else {
						return(NA)
					}
					
				}
				
			})
	cat('100%!\n')
	n.col <- max(sapply(res, length))
	res <- lapply(res, function (x) {
				if (length(x) == 1) {
					return(rep(NA, n.col)) 
				} else {
					return(x)
				}
			} )
	res <- matrix(unlist(res), length(res), n.col, byrow=TRUE)
	col.names <- gsub('count.', 'abund.', col.names)
	col.names <- gsub('zero.', 'prev.', col.names)
	rownames(res) <- rownames(otu.tab.win.sel)
	colnames(res) <- col.names
	res <- data.frame(res, method = method)
	res$method <- as.character(res$method)
# Step 5: Error handling

   if (sum(is.na(res$p.value)) != 0) {
	   if (!is.null(fail.method)) {
		   if (fail.method == 'permutation') {
			   if (!is.null(grp.name)) {
				   cat('Handle failed taxa using permutation test!\n')
				   prop <- t(t(otu.tab.win.sel) / size.factor)
				   prop <- prop[is.na(res$p.value), , drop = FALSE]
				   obj <- permute_differential_analysis(meta.dat, prop, grp.name=grp.name, adj.name = adj.name)
				   res[is.na(res$p.value), 'method'] <- 'permutation(n=1000)'
				   res[is.na(res$p.value), 'p.value'] <- obj$p.raw
			   } else {
				   warning('Handle failed taxa requires grp.name and adj.name!\n')
			   }
		   } else {
			   warning('The specified method for failed taxa is not supported currently!\n')
		   }
	   }
   }
  cat('Completed!\n')
	return(list(call = match.call(), result = res, size.factor=size.factor, test.ind = test.ind))
	
}
