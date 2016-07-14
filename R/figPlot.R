figPlot <- function(crcInfo, ylimit, xlabel = "log[concentration, mol/L]", ylabel = "Inhibition [%]", lgd = NULL){
	# plot the concentration-response curves
	#tiff(file = paste(root_name, "_04-12.tiff", sep = ""), res = 100)
	size <- dim(crcInfo)
	x <- crcInfo[, 1]
	yhat <- crcInfo[, 2]
	expr <- crcInfo[, 3 : (size[2] - 4)]
	if(is.vector(expr)) expr <- as.matrix(expr)
	oci <- crcInfo[, (size[2] - 3) : (size[2] - 2)]
	fci <- crcInfo[,  (size[2] - 1) : size[2]]
	
	if(missing(ylimit)) ylimit <- c((min(expr) * 100 -20), (max(expr) * 100 + 20))
	par(mar=c(5,5,1,1))
	plot(rep(log10(x), ncol(expr)), expr * 100, ylim = ylimit, pch = 16, xlab = xlabel, ylab = ylabel, cex = 1.8, cex.lab = 1.8, cex.axis = 1.8)
	lines(log10(x), yhat * 100, col = 1, lwd = 1.9)
	lines(log10(x), oci[, 1] * 100, col = 'blue', lwd = 1.9)
	lines(log10(x), oci[, 2] * 100, col = 'blue', lwd = 1.9)
	lines(log10(x), fci[, 1] * 100, col = 'red', lwd = 1.9)
	lines(log10(x), fci[, 2] * 100, col = 'red', lwd = 1.9)
	
	if (is.null(lgd) == FALSE) {
		legend('topleft', lgd, cex = 2)
	}
	#legend("topleft", inset = 0.01, root_name, box.col = 'white', cex = 1.9) 
	#dev.off()
}