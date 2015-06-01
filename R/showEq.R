showEq <- function(){
	# default Non-linear least square fitting algorithm is the Gauss-Newton algorithm
	
	print('######################## function selection ###########################')
	print("Hill: y ~ Alpha * x / (Beta + x)")
	print("Weibull: y ~ 1 - exp(-exp(Alpha + Beta * log10(x)))")
	print("Logit: y ~ 1/(1 + exp((-Alpha)- Beta * log10(x)))")
	print("BCW(Box-Cox-Weibull): y ~ 1 - exp(-exp(Alpha + Beta * ((x^Gamma - 1) / Gamma)))")
	print("BCL(Box-Cox-Logit): y ~ (1 + exp(-Alpha - Beta *((x^Gamma - 1) / Gamma)))^(-1)")
	print("GL(Generalized Logit): y ~ 1 / (1 + exp(-Alpha - Beta * log10(x)))^Gamma")

}