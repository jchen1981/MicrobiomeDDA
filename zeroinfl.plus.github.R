# Title: An omnibust for different abundance analysis of microbiome-seq data
# Version: 0.1
# Authors: Jun Chen (chen.jun2@mayo.edu)
# Description: It implements a general regression framework allowing the prevalance, 
# abundance, and dispersion to depend on covariates. Existing packages does not allowing 
# covariate-dependent dispersion, which could lead to either reduced power or inflated
# type I error if the heterogenety is not taken into account
# Date: 2017/02/07

#####################################################################################
#                              Instructions
# The code was adapted on the function 'zeroinfl' from 'pscl' package.
# It has a similar interface as 'zeroinfl' for the main functions 'zinb.reg' 
# Please see 'zeroinfl' for more details.
# In this implementation, we allow covariate-dependent dispersion, which can be used to
# model data heterogeneity. 
# Thus the only difference is in the specification of the formula, which allows the model for 
# dispersion component. A typical use is as follows:
#
# Y ~ Group + offset(offsetx) | Group + offset(offsetz) | Batch + offset(offsetm),#
#
# where we specify models for count (abundance), zero (prevalence) and dispersion part. 
# We also provide likelihood ratio test (zinb.lrt) for different models.
#####################################################################################

model_offset_2 <- function (x, terms = NULL, offset = TRUE)  {
	if (is.null(terms)) 
		terms <- attr(x, "terms")
	offsets <- attr(terms, "offset")
	if (length(offsets) > 0) {
		ans <- if (offset) 
					x$"(offset)"
				else NULL
		if (is.null(ans)) 
			ans <- 0
		for (i in offsets) ans <- ans + x[[deparse(attr(terms, 
									"variables")[[i + 1]])]]
		ans
	} else {
		ans <- if (offset) 
					x$"(offset)"
				else NULL
	}
	if (!is.null(ans) && !is.numeric(ans)) 
		stop("'offset' must be numeric")
	ans
}

zinb.control <- function (maxIter = 200, tol = 1e-5, trace = TRUE, start = NULL) {
	rval <- list(maxIter = maxIter, trace = trace, tol = tol, start = start)
	rval
}

# Zeroinflated negative binomial regression with covariate-dependent dispersion
# Formula can be specified for three components (count, zero, and dispersion) respectively
# e.g. Y ~ Group + offset(offsetx) | Group + offset(offsetz) | Batch + offset(offsetm).
# The code was adapted on the skeleton of the 'zeroinfl' function; Please see 'zeroinfl' usage for details!
zinb.reg <- function (formula, data, subset, na.action, weights, offset, 
		link = c("logit", "probit", "cloglog", "cauchit", "log"), control = zinb.control(), 
		model = TRUE, y = TRUE, x = FALSE, ...) {
	
	# Augmented likelihood
	loglik <- function (paras, p.ind, pi.ind, phi.ind, X.p, X.pi, X.phi, y, a, offsetz, offsetx, offsetm) {
		# Make sure all the coefficients are in matrix format
		beta.p <- paras[p.ind]
		beta.pi <- paras[pi.ind]
		beta.phi <- paras[phi.ind]
		
		eta.p <- as.vector(X.p %*% beta.p + offsetz)
		eta.pi <- as.vector(X.pi %*% beta.pi + offsetx)
		eta.phi <- as.vector(X.phi %*% beta.phi + offsetm)
		
		p <- linkinv(eta.p)
		mu <- exp(eta.pi)
		phi <- exp(eta.phi)
		
		- sum(a * log(1 - p) + (1 - a) * log(p) + (1 - a) * dnbinom(y, mu = mu, size = phi, log = TRUE))
	}
	
	# Original likelihood
	loglik0 <- function (paras, p.ind, pi.ind, phi.ind, X.p, X.pi, X.phi, y, offsetz, offsetx, offsetm) {
		# Make sure all the coefficients are in matrix format
		beta.p <- paras[p.ind]
		beta.pi <- paras[pi.ind]
		beta.phi <- paras[phi.ind]
		
		eta.p <- as.vector(X.p %*% beta.p + offsetz)
		eta.pi <- as.vector(X.pi %*% beta.pi + offsetx)
		eta.phi <- as.vector(X.phi %*% beta.phi + offsetm)
		
		p <- linkinv(eta.p)
		mu <- exp(eta.pi) 
		phi <- exp(eta.phi)
		
		I <- as.numeric(y == 0)
		sum(log((1 - p) * I + p * dnbinom(y, mu = mu, size = phi)))
	}
	
	
	Estep <- function (paras, p.ind, pi.ind, phi.ind, X.p, X.pi, X.phi, y, offsetz, offsetx, offsetm) {
		# Make sure all the coefficients are in matrix format
		# Make sure y not all zeros or nonzeros
		beta.p <- paras[p.ind]
		beta.pi <- paras[pi.ind]
		beta.phi <- paras[phi.ind]
		
		eta.p <- as.vector(X.p %*% beta.p + offsetz)
		eta.pi <- as.vector(X.pi %*% beta.pi + offsetx)
		eta.phi <- as.vector(X.phi %*% beta.phi + offsetm)
		
		# Use more generic link function
		p <- linkinv(eta.p)
		mu <- exp(eta.pi)
		phi <- exp(eta.phi)
		
		a <- numeric(length(y))
		ind <- y == 0
		
		y <- y[ind]
		p <- p[ind]
		mu <- mu[ind]
		phi <- phi[ind]
		
		a[ind] <- (1 - p)  / (1 - p  + p * dnbinom(y, mu = mu, size = phi))
		a
	}
	
	Mstep <- function (obj.func, beta0, p.ind, pi.ind, phi.ind, X.p, X.pi, X.phi, y, a, offsetz, offsetx, offsetm, ...) {
		# Need to add error control for nlm
		nlm.obj <- nlm(obj.func, p=beta0, p.ind=p.ind, pi.ind=pi.ind, phi.ind=phi.ind, 
				X.p=X.p, X.pi=X.pi, X.phi=X.phi, y=y, a=a, offsetz=offsetz, offsetx=offsetx, offsetm=offsetm, ...)
		nlm.obj$Q1 <- loglik(nlm.obj$estimate, p.ind, pi.ind, phi.ind, X.p, X.pi, X.phi, y, a, offsetz, offsetx, offsetm)
		nlm.obj
	}
	
	linkstr <- match.arg(link)
	linkobj <- make.link(linkstr)
	linkinv <- linkobj$linkinv
	
	if (control$trace) {
		cat("Zero-inflated Count Model\n", paste("count model: negative binomial model",
						"with log link\n"), paste("zero-inflation model: binomial with", 
						linkstr, "link\n"), paste("dispersion model: log link\n"), sep = "")
	}

	cl <- match.call()
	
	if (missing(data)) {
		data <- environment(formula)
	}
	
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], 
			as.name("|"))) {
		ff <- formula
		formula[[3]][1] <- call("+")
		formula[[3]][[2]][1] <- call("+")
		mf$formula <- formula
		ffc <- . ~ .
		ffz <- ~.
		ffm <- ~.
		ffc[[2]] <- ff[[2]]
		ffc[[3]] <- ff[[3]][[2]][[2]]
		ffz[[3]] <- ff[[3]][[2]][[3]]
		ffz[[2]] <- NULL
		ffm[[3]] <- ff[[3]][[3]]
		ffm[[2]] <- NULL
	} else {
		ffz <- ffc <- ffm <- ff <- formula
		ffz[[2]] <- ffm[[2]] <- NULL
	}
	
	if (inherits(try(terms(ffz), silent = TRUE), "try-error")) {
		ffz <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])), 
								deparse(ffz))))
	}
	if (inherits(try(terms(ffm), silent = TRUE), "try-error")) {
		ffm <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])), 
								deparse(ffm))))
	}
	
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	
	mtX <- terms(ffc, data = data)
	X <- model.matrix(mtX, mf)
	
	mtZ <- terms(ffz, data = data)
	mtZ <- terms(update(mtZ, ~.), data = data)
	Z <- model.matrix(mtZ, mf)
	
	mtM <- terms(ffm, data = data)
	mtM <- terms(update(mtM, ~.), data = data)
	M <- model.matrix(mtM, mf)
	
	Y <- model.response(mf, "numeric")
	
	if (length(Y) < 1) {
		stop("empty model")
	}
	if (all(Y > 0)) {
		stop("invalid dependent variable, minimum count is not zero")
	}
	
	if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001))))) {
		stop("invalid dependent variable, non-integer values")
	}
	
	Y <- as.integer(round(Y + 0.001))
	
	if (any(Y < 0)) {
		stop("invalid dependent variable, negative counts")
	}
	
	if (sum(Y == 0) == 0) {
		stop("invalid dependent variable, no zeros")
	}
	
	if (control$trace) {
		cat("dependent variable:\n")
		tab <- table(Y)
		names(dimnames(tab)) <- NULL
		print(tab)
	}
	
	n <- length(Y)
	kx <- NCOL(X)
	kz <- NCOL(Z)
	km <- NCOL(M)
	Y0 <- Y <= 0
	Y1 <- Y > 0
	
	
	weights <- model.weights(mf)
	
	if (is.null(weights)) {
		weights <- 1
	}
	if (length(weights) == 1) {
		weights <- rep.int(weights, n)
	}
	weights <- as.vector(weights)
	names(weights) <- rownames(mf)
	
	offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE)
	if (is.null(offsetx)) {
		offsetx <- 0
	}
	if (length(offsetx) == 1) {
		offsetx <- rep.int(offsetx, n)
	}
	offsetx <- as.vector(offsetx)
	
	
	offsetz <- model_offset_2(mf, terms = mtZ, offset = FALSE)
	if (is.null(offsetz)) {
		offsetz <- 0
	}
	if (length(offsetz) == 1) {
		offsetz <- rep.int(offsetz, n)
	}	
	offsetz <- as.vector(offsetz)
	
	offsetm <- model_offset_2(mf, terms = mtM, offset = FALSE)
	if (is.null(offsetm)) {
		offsetm <- 0
	}
	if (length(offsetm) == 1) {
		offsetm <- rep.int(offsetm, n)
	}	
	offsetm <- as.vector(offsetm)
	
	
	start <- control$start
	if (!is.null(start)) {
		valid <- TRUE
		if (!("count" %in% names(start))) {
			valid <- FALSE
			warning("invalid starting values, count model coefficients not specified")
			start$count <- rep.int(0, kx)
		}
		if (!("zero" %in% names(start))) {
			valid <- FALSE
			warning("invalid starting values, zero-inflation model coefficients not specified")
			start$zero <- rep.int(0, kz)
		}
		if (!("theta" %in% names(start))) {
			valid <- FALSE
			warning("invalid starting values, dispersion model coefficients not specified")
			start$zero <- rep.int(0, km)
		}
		if (length(start$count) != kx) {
			valid <- FALSE
			warning("invalid starting values, wrong number of count model coefficients")
		}
		if (length(start$zero) != kz) {
			valid <- FALSE
			warning("invalid starting values, wrong number of zero-inflation model coefficients")
		}
		if (length(start$theta) != km) {
			valid <- FALSE
			warning("invalid starting values, wrong number of dispersion coefficients")
		}
		
		if (!valid) {
			start <- NULL
		}	
	}
	
	if (is.null(start)) {
		if (control$trace) {
			cat("generating starting values...\n")
		}
		# Using zeroinfl 
		m0 <- zeroinfl(Y ~ 1 | 1 , offset = offsetx, dist = "negbin", EM = TRUE)
		start$zero <- c(m0$coefficients$zero, rep(0, kz - 1))
		start$count <- c(m0$coefficients$count, rep(0, kx - 1))
		start$theta <- c(log(m0$theta), rep(0, km - 1))
		
	}
	if (control$trace) {
		cat("EM estimation:\n")
	}
	
	p.ind <- 1:kz
	pi.ind <- (kz + 1) : (kz + kx)
	phi.ind <- (kz + kx + 1) : (kz + kx + km)
	
	beta0 <- c(start$zero, start$count, start$theta)
	a0 <- 1 - linkinv(start$zero[1] + offsetz)
	
	a0[Y != 0] <- 0
	Q0 <- loglik(beta0, p.ind, pi.ind, phi.ind, Z, X, M, Y, a = a0, offsetz, offsetx, offsetm)
	
	nIter <- 0
	code.all <- NULL
	while (TRUE) {
		if (control$trace) {
			cat(Q0, '\n')
		}
		nIter <- nIter + 1
		a1 <- Estep(beta0, p.ind, pi.ind, phi.ind, Z, X, M, Y, offsetz, offsetx, offsetm)
		M.obj <- Mstep(loglik, beta0, p.ind, pi.ind, phi.ind, Z, X, M, Y, a1, offsetz, offsetx, offsetm, hessian=FALSE, ...)
		code.all <- union(code.all, M.obj$code)
		beta1 <- M.obj$estimate
		Q1 <- M.obj$Q1
		
		if (abs(Q1 - Q0) / abs(Q0) < control$tol | nIter >= control$maxIter) break
		
		Q0 <- Q1
		beta0 <- beta1
	}
	if (control$trace) {
		cat('Finished!\n')
	}
	
	# obtain the hessian
	M.obj <- Mstep(loglik, beta0, p.ind, pi.ind, phi.ind, Z, X, M, Y, a1, offsetz, offsetx, offsetm, hessian=TRUE, ...)
	
	fit <- M.obj
	fit$nIter <- nIter
	fit$code.all <- code.all
	fit$loglik <- loglik0(beta1, p.ind, pi.ind, phi.ind, Z, X, M, Y, offsetz, offsetx, offsetm)

	if (sum(!(fit$code.all %in% c(1, 2)))) {
		warning("optimization failed to converge in some iterations!\n")
		fit$converged <- FALSE
	} else {
		fit$converged <- TRUE
	}
	
	coefz <- fit$estimate[1:kz]
	names(coefz) <- names(start$zero) <-colnames(Z)
	
	coefc <- fit$estimate[(kz + 1):(kx + kz)]
	names(coefc) <- names(start$count) <- colnames(X)
	
	coefm <- fit$estimate[(kx + kz + 1):(kx + kz + km)]
	names(coefm) <- names(start$theta) <- colnames(M)
	
	vc <- solve(as.matrix(fit$hessian))
	
	colnames(vc) <- rownames(vc) <- c(paste("zero", colnames(Z), sep = "."),
			paste("count", colnames(X), sep = "."),
			paste("dispersion", colnames(M), sep = "."))
	
	renames <- c(paste("count", colnames(Z), sep = "."),
			paste("zero", colnames(X), sep = "."),
			paste("dispersion", colnames(M), sep = "."))
	vc <- vc[renames, renames]
	
	mu <- exp(X %*% coefc + offsetx)[, 1]
	phi <- linkinv(Z %*% coefz + offsetz)[, 1]
	theta <- exp(M %*% coefm + offsetm)[, 1]
	
	Yhat <- (1 - phi) * mu
	res <- sqrt(weights) * (Y - Yhat)
	nobs <- sum(weights > 0)
	rval <- list(coefficients = list(count = coefc, zero = coefz, dispersion = coefm), 
			residuals = res, fitted.values = Yhat, nlmfit = fit, 
			control = control, start = start, weights = if (identical(as.vector(weights), 
							rep.int(1L, n))) NULL else weights, offset = list(count = if (identical(offsetx, 
									rep.int(0, n))) NULL else offsetx, zero = if (identical(offsetz, 
									rep.int(0, n))) NULL else offsetz, dispersion = if (identical(offsetm, 
									rep.int(0, n))) NULL else offsetm), n = nobs, 
			df.null = nobs - 3, df.residual = nobs - (kx + kz + km), 
			terms = list(count = mtX, zero = mtZ, dispersion = mtM, full = mt),  loglik = fit$loglik, vcov = vc, 
			link = linkstr, linkinv = linkinv, converged = fit$converged, call = cl, formula = ff, levels = .getXlevels(mt, 
					mf), contrasts = list(count = attr(X, "contrasts"), 
					zero = attr(Z, "contrasts"), dispersion = attr(M, "contrasts")))
	if (model) {
		rval$model <- mf
	}
	
	if (y) {
		rval$y <- Y
	}
	
	if (x) {
		rval$x <- list(count = X, zero = Z, dispersion = M)
	}
	class(rval) <- "zinb.reg"
	return(rval)
}

nb.control <- function (maxIter = 200, tol = 1e-5, trace = TRUE, start = NULL) {
	rval <- list(maxIter = maxIter, trace = trace, tol = tol, start = start)
	rval
}

# Negative binomial regression with covariate-dependent dispersion
# Formula can be specified for two components (count and dispersion) respectively
# e.g. Y ~ Group + offset(offsetx) | Batch + offset(offsetm).
nb.reg <- function (formula, data, subset, na.action, weights, offset, 
		link = c("logit", "probit", "cloglog", "cauchit", "log"), control = nb.control(), 
		model = TRUE, y = TRUE, x = FALSE, ...) {
	
	# Augmented likelihood
	loglik <- function (paras, pi.ind, phi.ind, X.pi, X.phi, y, offsetx, offsetm) {

		beta.pi <- paras[pi.ind]
		beta.phi <- paras[phi.ind]
		
		eta.pi <- as.vector(X.pi %*% beta.pi + offsetx)
		eta.phi <- as.vector(X.phi %*% beta.phi + offsetm)

		mu <- exp(eta.pi)
		phi <- exp(eta.phi)
		
		- sum(dnbinom(y, mu = mu, size = phi, log = TRUE))
	}
	
	if (control$trace) {
		cat("Zero-inflated Count Model\n", paste("count model: negative binomial model",
						"with log link\n"), paste("dispersion model: log link\n"), sep = "")
	}
	
	cl <- match.call()
	
	if (missing(data)) {
		data <- environment(formula)
	}
	
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
	mf <- mf[c(1, m)]
	mf$drop.unused.levels <- TRUE
	
	if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], 
			as.name("|"))) {
		ff <- formula
		formula[[3]][1] <- call("+")
		mf$formula <- formula
		ffc <- . ~ .
		ffm <- ~.
		ffc[[2]] <- ff[[2]]
		ffc[[3]] <- ff[[3]][[2]]
		ffm[[3]] <- ff[[3]][[3]]
		ffm[[2]] <- NULL
	} else {
		ffm <- ffc <- ff <- formula
		ffm[[2]] <- NULL
	}
	
	if (inherits(try(terms(ffm), silent = TRUE), "try-error")) {
		ffm <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])), 
								deparse(ffm))))
	}
	
	mf[[1]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	
	mtX <- terms(ffc, data = data)
	X <- model.matrix(mtX, mf)
	
	mtM <- terms(ffm, data = data)
	mtM <- terms(update(mtM, ~.), data = data)
	M <- model.matrix(mtM, mf)
	
	Y <- model.response(mf, "numeric")
	
	if (length(Y) < 1) {
		stop("empty model")
	}

	if (!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001))))) {
		stop("invalid dependent variable, non-integer values")
	}
	
	Y <- as.integer(round(Y + 0.001))
	
	if (any(Y < 0)) {
		stop("invalid dependent variable, negative counts")
	}
	
	if (control$trace) {
		cat("dependent variable:\n")
		tab <- table(Y)
		names(dimnames(tab)) <- NULL
		print(tab)
	}
	
	n <- length(Y)
	kx <- NCOL(X)
	km <- NCOL(M)

	weights <- model.weights(mf)
	
	if (is.null(weights)) {
		weights <- 1
	}
	if (length(weights) == 1) {
		weights <- rep.int(weights, n)
	}
	weights <- as.vector(weights)
	names(weights) <- rownames(mf)
	
	offsetx <- model_offset_2(mf, terms = mtX, offset = TRUE)
	if (is.null(offsetx)) {
		offsetx <- 0
	}
	if (length(offsetx) == 1) {
		offsetx <- rep.int(offsetx, n)
	}
	offsetx <- as.vector(offsetx)

	
	offsetm <- model_offset_2(mf, terms = mtM, offset = FALSE)
	if (is.null(offsetm)) {
		offsetm <- 0
	}
	if (length(offsetm) == 1) {
		offsetm <- rep.int(offsetm, n)
	}	
	offsetm <- as.vector(offsetm)
	
	
	start <- control$start
	if (!is.null(start)) {
		valid <- TRUE
		if (!("count" %in% names(start))) {
			valid <- FALSE
			warning("invalid starting values, count model coefficients not specified")
			start$count <- rep.int(0, kx)
		}

		if (!("theta" %in% names(start))) {
			valid <- FALSE
			warning("invalid starting values, dispersion model coefficients not specified")
			start$zero <- rep.int(0, km)
		}
		if (length(start$count) != kx) {
			valid <- FALSE
			warning("invalid starting values, wrong number of count model coefficients")
		}

		if (length(start$theta) != km) {
			valid <- FALSE
			warning("invalid starting values, wrong number of dispersion coefficients")
		}
		
		if (!valid) {
			start <- NULL
		}	
	}
	
	if (is.null(start)) {
		if (control$trace) {
			cat("generating starting values...\n")
		}
		# offset confusion, to circumvent, we use it as a covariate
		m0 <- glm.nb(Y ~ 1 + offsetx)
		start$count <- c(m0$coefficients[1], rep(0, kx - 1))
		start$theta <- c(log(m0$theta), rep(0, km - 1))
		
	}
	if (control$trace) {
		cat("begin model fitting...\n")
	}
	pi.ind <- 1 : kx
	phi.ind <- (kx + 1) : (kx + km)
	
	beta0 <- c(start$count, start$theta)

	nlm.obj <- nlm(loglik, p=beta0,  pi.ind=pi.ind, phi.ind=phi.ind, X.pi=X,  X.phi=M, y=Y, 
			offsetx=offsetx, offsetm=offsetm, hessian=TRUE, ...)

	fit <- nlm.obj
	fit$loglik <- - loglik(nlm.obj$estimate, pi.ind, phi.ind, X, M, Y, offsetx, offsetm)
	
	if (sum(!(fit$code %in% c(1, 2)))) {
		warning("optimization failed to converge in some iterations!\n")
		fit$converged <- FALSE
	} else {
		fit$converged <- TRUE
	}
	if (control$trace) {
		if (fit$converged) {
			cat("model converged...\n")
		} else {
			cat("model not converged...\n")
		}
		
	}
	coefc <- fit$estimate[(1):(kx)]
	names(coefc) <- names(start$count) <- colnames(X)
	
	coefm <- fit$estimate[(kx + 1):(kx + km)]
	names(coefm) <- names(start$theta) <- colnames(M)
	vc <- solve(as.matrix(fit$hessian))
	
	colnames(vc) <- rownames(vc) <- 
			c(paste("count", colnames(X), sep = "."),
			paste("dispersion", colnames(M), sep = "."))
	
	mu <- exp(X %*% coefc + offsetx)[, 1]
	theta <- exp(M %*% coefm + offsetm)[, 1]
	
	Yhat <- mu
	res <- sqrt(weights) * (Y - Yhat)
	nobs <- sum(weights > 0)
	rval <- list(coefficients = list(count = coefc, dispersion = coefm), 
			residuals = res, fitted.values = Yhat, nlmfit = fit, 
			control = control, start = start, weights = if (identical(as.vector(weights), 
							rep.int(1L, n))) NULL else weights, offset = list(count = if (identical(offsetx, 
									rep.int(0, n))) NULL else offsetx,  dispersion = if (identical(offsetm, 
									rep.int(0, n))) NULL else offsetm), n = nobs, 
			df.null = nobs - 2, df.residual = nobs - (kx + km), 
			terms = list(count = mtX, dispersion = mtM, full = mt),  loglik = fit$loglik, vcov = vc, 
			converged = fit$converged, call = cl, formula = ff, levels = .getXlevels(mt, 
					mf), contrasts = list(count = attr(X, "contrasts"), dispersion = attr(M, "contrasts")))
	if (model) {
		rval$model <- mf
	}
	
	if (y) {
		rval$y <- Y
	}
	
	if (x) {
		rval$x <- list(count = X, dispersion = M)
	}
	class(rval) <- "nb.reg"
	return(rval)
}


# LRT for prevalence and/or abundance and/or dispersion 
# Wald test has inflated type I error if dispersion is tested
zinb.lrt <- function (formula.H1, formula.H0, data, ...) {
	
	mod1 <- zinb.reg(formula=formula.H1, data=data, ...)
	mod0 <- zinb.reg(formula=formula.H0, data=data, ...)
	p.value <- pchisq(2 * (mod1$loglik - mod0$loglik), df = mod0$df.residual - mod1$df.residual, lower.tail = FALSE)
	# Extract variables of interest
	var.name <- setdiff(colnames(mod1$vcov), colnames(mod0$vcov))
	coef <- unlist(mod1$coefficients)[var.name]
	se <- sqrt(diag(mod1$vcov[var.name, var.name, drop=FALSE]))
    se[is.nan(se)] <- NA
	zstat <- coef / se
#	pval <- 2 * pnorm(-abs(zstat))
	coef <- cbind(coef, se, coef - 1.96 * se, coef + 1.96 * se, zstat)
	colnames(coef) <- c("Estimate", "Std. Error", "2.5%", "97.5%",  "z value")
	list(mod.H1 = mod1, mod.H0 = mod0, p.value = p.value, coef = coef)
}

# LRT for prevalence and/or abundance and/or dispersion 
# Wald test has inflated type I error if dispersion is tested
nb.lrt <- function (formula.H1, formula.H0, data, ...) {
	
	mod1 <- nb.reg(formula=formula.H1, data=data, ...)
	mod0 <- nb.reg(formula=formula.H0, data=data, ...)
	p.value <- pchisq(2 * (mod1$loglik - mod0$loglik), df = mod0$df.residual - mod1$df.residual, lower.tail = FALSE)
	# Extract variables of interest
	var.name <- setdiff(colnames(mod1$vcov), colnames(mod0$vcov))
	coef <- unlist(mod1$coefficients)[var.name]
	se <- sqrt(diag(mod1$vcov[var.name, var.name, drop=FALSE]))
	se[is.nan(se)] <- NA
	zstat <- coef / se
#	pval <- 2 * pnorm(-abs(zstat))
	coef <- cbind(coef, se, coef - 1.96 * se, coef + 1.96 * se, zstat)
	colnames(coef) <- c("Estimate", "Std. Error", "2.5%", "97.5%",  "z value")
	list(mod.H1 = mod1, mod.H0 = mod0, p.value = p.value, coef = coef)
}

# Bootstrap test for zeroinflation
# If the data are not truly zeroinflated, using ZINB can be less powerful than using NB;  vice versa.
# Based on simulation studies, both NB/ZINB have conservative type I error control under model under/over-specification
# Therefore, selecting a more appropriate model is more about increasing the power
# The current test can be used as a screening step
zeroinfl.test <- function (y, size.factor, B=999) {
	zero.obs <- sum(y == 0)
	if (zero.obs == 0) return(1)
	n <- length(y)
	obj <- glm.nb(y ~ 1 + offset(log(size.factor)))
	# Not converged usually indicate severe zeroinflation
	if (!obj$converged) {
		return(1 / (B + 1))
	} 
	mu <- exp(obj$coefficient + log(size.factor))
	size <- obj$theta
	y.bt <- matrix(rnbinom(n * B, mu = mu, size = size), ncol=B)
	return(mean(c(colSums(y.bt == 0) >= zero.obs, 1)))
}

