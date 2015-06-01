curveFit <- function(x, expr, eq = c("Hill", "Weibull", "Logit", "BCW", "BCL", "GL"),
			param, 	effv, sigLev = 0.05, noec = TRUE, fig = TRUE,
			algo = "default"){
	# NLS curve fitting for sigmoidal equations
	# x is a vector 
	# y is a vector or matrix
	
	## Jacobian matrix calculation
	jacobian <- function(eq, x, paraHat){
		n <- length(x)
		mpara <- length(paraHat)
		Alpha <- paraHat[1]
		Beta <- paraHat[2]
		if (mpara == 3) Gamma <- paraHat[3]
		jac <- matrix(rep(0, n * mpara), n, mpara)
		
		jacFun <- switch(eq,
			Hill = c('x / (Beta + x)', '-Alpha * x / (Beta + x)^2'),
			Weibull = c('exp(Alpha + Beta * log(x) / log(10)) * exp(-exp(Alpha + Beta * log(x) / log(10)))', 
						'log(x) / log(10) * exp(Alpha + Beta * log(x) / log(10)) * exp(-exp(Alpha + Beta * log(x) / log(10)))'),
			Logit = c('1 / (1 + exp(-Alpha - Beta * log(x) / log(10)))^2 * exp(-Alpha - Beta * log(x) / log(10))', 
						'1 / (1 + exp(-Alpha - Beta * log(x) / log(10)))^2 * log(x) / log(10) * exp(-Alpha - Beta * log(x) / log(10))'),
			BCW = c('exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))', 
					'(x^Gamma - 1) / Gamma * exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))', 
					'(Beta * x^Gamma * log(x) / Gamma - Beta * (x^Gamma - 1) / Gamma^2) * exp(Alpha + Beta * (x^Gamma - 1) / Gamma) * exp(-exp(Alpha + Beta * (x^Gamma - 1) / Gamma))'),
			BCL = c('1 /(1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma))^2 * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)', 
					'1 / (1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma))^2 * (x^Gamma - 1) / Gamma * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)', 
					'-1 / (1 + exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)) ^ 2 * (-Beta * x^Gamma * log(x) / Gamma + Beta * (x^Gamma - 1) / Gamma^2) * exp(-Alpha - Beta * (x^Gamma - 1) / Gamma)'),
			GL = c('1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * Gamma * exp(-Alpha - Beta * log(x) / log(10)) / (1 + exp(-Alpha - Beta * log(x) / log(10)))', 
					'1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * Gamma * log(x) / log(10) * exp(-Alpha - Beta * log(x) / log(10)) / (1 + exp(-Alpha - Beta * log(x) / log(10)))', 
					'-1 / ((1 + exp(-Alpha - Beta * log(x) / log(10)))^Gamma) * log(1 + exp(-Alpha - Beta * log(x) / log(10)))')
		)
		
		for (i in seq(mpara)) jac[, i] <- eval(parse(text = jacFun[i]))
		return(jac)
	}
	
	#############################################################
	ecxCI <- function(ciInfo, effv){
		# effect concentration and associated confidence intervals calculation
		oci.up <- spline(ciInfo[, 3], ciInfo[, 1], method = 'fmm', xout = effv)$y # effect concentration upper bound
		oci.low <- spline(ciInfo[, 4], ciInfo[, 1], method = 'fmm', xout = effv)$y # effect concentration lower bound
		oci.low[which(oci.low < 0)] = 0
		
		fci.up <- spline(ciInfo[, 5], ciInfo[, 1], method = 'fmm', xout = effv)$y # effect concentration upper bound
		fci.low <- spline(ciInfo[, 6], ciInfo[, 1], method = 'fmm', xout = effv)$y # effect concentration lower bound
		fci.low[which(fci.low < 0)] = 0
		
		ec.CI <- cbind(oci.low, oci.up, fci.low, fci.up)
		colnames(ec.CI) <- c('OCI.low', 'OCI.up', 'FCI.low', 'FCI.up')
		return(ec.CI)
	}
	
	#############################################################
	## confidence intervals for effect
	effvCI <- function(ciInfo, effv, ecx){
		# confidence interval for effect based on spline interpolation
		eoci.low <- spline(ciInfo[, 1], ciInfo[, 3], method = 'fmm', xout = ecx)$y
		eoci.up <- spline(ciInfo[, 1], ciInfo[, 4], method = 'fmm', xout = ecx)$y
		eoci.low[which(eoci.low < 0)] = 0
		
		efci.up <- spline(ciInfo[, 1], ciInfo[, 6], method = 'fmm', xout = ecx)$y # effect concentration upper bound
		efci.low <- spline(ciInfo[, 1], ciInfo[, 5], method = 'fmm', xout = ecx)$y # effect concentration lower bound
		efci.low[which(efci.low < 0)] = 0
		
		effv.CI <- cbind(eoci.low, eoci.up, efci.low, efci.up)
		colnames(effv.CI) <- c('eOCI.low', 'eOCI.up', 'eFCI.low', 'eFCI.up')
		return(effv.CI)
	}

	#############################################################
	figPlot <- function(crcInfo, xl = "lg[concentration, mol/L]", yl = "Inhibition [%]"){
		# plot the concentration-response curves
		#tiff(file = paste(root_name, "_04-12.tiff", sep = ""), res = 100)
		size <- dim(crcInfo)
		x <- crcInfo[, 1]
		yhat <- crcInfo[, 2]
		expr <- crcInfo[, 3 : (size[2] - 4)]
		if(is.vector(expr)) expr <- as.matrix(expr)
		oci <- crcInfo[, (size[2] - 3) : (size[2] - 2)]
		fci <- crcInfo[,  (size[2] - 1) : size[2]]
		
		par(mar=c(5,5,1,1))
		plot(rep(log10(x), ncol(expr)), expr * 100, , ylim = c(-10, 110), pch = 16, xlab = xl, ylab = yl, cex = 1.8, cex.lab = 1.8, cex.axis = 1.8)
		lines(log10(x), yhat * 100, col = 1, lwd = 1.9)
		lines(log10(x), oci[, 1] * 100, col = 'blue', lwd = 1.9)
		lines(log10(x), oci[, 2] * 100, col = 'blue', lwd = 1.9)
		lines(log10(x), fci[, 1] * 100, col = 'red', lwd = 1.9)
		lines(log10(x), fci[, 2] * 100, col = 'red', lwd = 1.9)
		
		#legend("topleft", inset = 0.01, root_name, box.col = 'white', cex = 1.9) 
		#dev.off()
	}
	#############################################################	
	## source('ECx.R')
	
	## checking experimental data (expr)
	
	if (missing(x) || missing(expr) || missing(eq) || missing(param)) stop('argument missing')
	
	n <- length(x)
	if (is.vector(expr)){
		if (n != length(expr)) stop("x and y should be in the same length")
		y <- expr
		expr <- as.matrix(expr)
		nrep <- 1		
	}else if (is.matrix(expr)){
		size <- dim(expr)
		nrep <- size[2]
		y <- rowMeans(expr)
		if(n != size[1]) stop("x and dim(y)[1] should be in the same length")
	}
	
	## deploying the equation
	fun <- switch(eq,
		Hill = 'y ~ Alpha * x / (Beta + x)',
		Weibull = 'y ~ 1 - exp(-exp(Alpha + Beta * log10(x)))',
		Logit = 'y ~ 1/(1 + exp((-Alpha) - Beta * log10(x)))',
		BCW = 'y ~ 1 - exp(-exp(Alpha + Beta * ((x^Gamma - 1) / Gamma)))',
		BCL = 'y ~ (1 + exp(-Alpha - Beta *((x^Gamma - 1) / Gamma)))^(-1)',
		GL = 'y ~ 1 / (1 + exp(-Alpha - Beta * log10(x)))^Gamma'
	)

	
	## checking nls2 package, use the nls2 or built-in nls for curve fitting
	#if(require(nls2)){
	if(requireNamespace("nls2", quietly = TRUE)){
		print("use the nls2 package")
		if(eq == "Hill" || eq == "Weibull" || eq == "Logit"){
			m <- 2 # the number of parameters
			mode(param) <- "numeric"
			fit <- nls2::nls2(fun, start = list(Alpha = param[1], Beta = param[2]), control = nls.control(maxiter = 1000), algorithm = algo)
		}else if(eq == "BCW" || eq == "BCL" || eq == "GL"){
			m <- 3 # the number of parameters
			mode(param) <- "numeric"
			fit <- nls2::nls2(fun, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3]), control = nls.control(maxiter = 1000), algorithm = algo)
		}
		#detach(package: nls2)
	}else {
		## use the built-in nls
		if(eq == "Hill" || eq == "Weibull" || eq == "Logit"){
			m <- 2 # the number of parameters
			mode(param) <- "numeric"
			fit <- nls(fun, start = list(Alpha = param[1], Beta = param[2]), control = nls.control())
		}else if(eq == "BCW" || eq == "BCL" || eq == "GL"){
			m <- 3 # the number of parameters
			mode(param) <- "numeric"
			fit <- nls(fun, start = list(Alpha = param[1], Beta = param[2], Gamma = param[3]), control = nls.control())
		}
	}
	
	fitInfo <- summary(fit) # fitting information
	
	yhat <- predict(fit, x) # y prediction
	sst <- sum((y - mean(y))^2) # total sum of squares
	sse <- sum((y - yhat)^2) # sum of squared errors
	r2 <- 1 - sse / sst # coefficient of determination
	adjr2 <- 1 - sse * (n - 1) / (sst * (n - m)) # adjusted coefficient of determination
	rmse <- sqrt(sse / (n - m)) # root-mean-square error
	mae <- sum(abs(y - yhat)) / n # mean absolute error
	aic <- n * log10(2 * pi) + n * log10(sse / n) + n + 2 * (m + 1) # Akaike information criterion 
	sta <- t(c(r2, adjr2, mae, rmse, aic))
	colnames(sta) <- c('r2', 'adjr2', 'MAE', 'RMSE', 'AIC')
	
	paramHat <- t(as.matrix(summary(fit)$parameters[, 1]))	
	jac <- jacobian(eq, x, paramHat) # jacobian matrix calculation
	probT <- qt(1 - sigLev / 2, n - m) # the student t distribution
	mse <- rmse^2  # squared residual standard error
	covPara <- mse * solve(t(jac) %*% jac)  # covariance matrix of the parameter estimates
	
	gap.OCI <- sqrt(mse + diag(jac %*% covPara %*% t(jac))) # observation based confidence intervals
	gap.FCI <- sqrt(diag(jac %*% covPara %*% t(jac))) # function based confidence intervals
	
	OCI.up <- yhat + probT * gap.OCI # OCI upper bound
	OCI.low <- yhat - probT * gap.OCI # OCI lower bound
	FCI.up <- yhat + probT * gap.FCI # FCI upper bound
	FCI.low <- yhat - probT * gap.FCI # FCI lower bound
	
	crcInfo <- cbind(x, yhat, expr, OCI.low, OCI.up, FCI.low, FCI.up)
	ciInfo <- cbind(x, yhat, OCI.low, OCI.up, FCI.low, FCI.up)
	
	## confidence intervals for effect concentration and effect
	## checking argument
	if(!missing(effv)) {
		## effect concentration and confidence intervals 
		ecx <- ECx(eq, paramHat, effv)
		ecx.ci <- ecxCI(ciInfo, effv)
		ecx.ci <- cbind(t(ecx), ecx.ci)
		
		## effect confidence intervals 
		effv.ci <- effvCI(ciInfo, effv, ecx)
		effv.vec <- t(t(effv))
		rownames(effv.vec) <- paste('E', effv * 100, sep = '')
		colnames(effv.vec) <- 'Effect'
		effv.ci <- cbind(effv.vec, effv.ci)
	}else{
		ecx.ci = NULL
		effv.ci = NULL
	}
	
	## non-observed effect concentration
	## least observed effect concentration
	## at least 3 repetition
	if(nrep >= 3 && noec == TRUE){
		## source('NOEC.R')
		noecInfo <- NOEC(x, expr, sigLev)	
	}else {
		noecInfo = NULL
	}
	
	## show concentration-response curve
	if(fig == TRUE){
		## source('figPlot.R')
		figPlot(crcInfo)
	}
	list(fitInfo = fitInfo, p = paramHat, sta = sta, crcInfo = crcInfo, eci = ecx.ci, effvci = effv.ci, noecInfo = noecInfo)
	
}
