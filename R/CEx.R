CEx <- function(model, param, conc){
	# calculate response based on concentration
	if (missing(model) || missing (param)) stop('argument missing')
	if (missing(conc)) conc = 0.00005
	if (is.vector(param)) param <- t(param)
	
	effv <- matrix(0, length(model), length(conc))
	
	for (i in seq(model)){
		fun <- model[i]
		p <- param[i, ]
		
		for (j in seq(conc)){
			if (fun == 'Hill')
				ev <- p[1] * conc[j] /(p[2] + conc[j])
			else if(fun == 'Weibull')
				ev <- 1 - exp(-exp(p[1] + p[2] * log10(conc[j])))
			else if (fun == "Logit")
				ev <- 1 / (1 + exp(-p[1] - p[2] * log10(conc[j])))
			else if (fun == "BCW")
				ev <- 1 - exp(-exp(p[1] + p[2] * ((conc[j]^p[3] - 1) / p[3])))
			else if (fun == "BCL")
				ev <- 1 / (1 + exp(-p[1] - p[2]((conc[j]^p[3] - 1) / p[3])))
			else if (fun == "GL")
				ev <- 1 / (1 + exp(-p[1] - p[2] * log10(conc[j])))^p[3]
			effv[i, j] <- ev
		}
	}
	
	colName <- paste('effect_@_', conc, sep = '')
	colnames(effv) <- colName
	rownames(effv) <- rownames(param)
	return(effv)
}
