genetic.LCBD <- function(Y, D.opt=2, perm.opt=2, sqrt.D=FALSE, nperm=999, save.D=FALSE, silent=FALSE)
#
# Compute and test significance of LCBD based on genetic distances.
# Y : Individual-based data. Rows = individuals, columns = alleles.
#     class = "genind".
#
# Distances available in dist.genpop() of {adegenet} (author: Thibaut Jombart):
#    method=1 : Nei (1972); not Euclidean.
#    method=2 : Edwards (1971); Cavalli-Sforza & Edwards (1967); Euclidean.
#               Note PL: The genetic chord D is the Hellinger D in ecology.
#    method=3 : Reynolds, Weir and Cockerman (1983); Euclidean.
#    method=4 : Rogers (1972); Euclidean.
#    method=5 : Prevosti (1975); not Euclidean.
#
# perm.opt = 1  # Individual-based: Permute whole rows of Y
# perm.opt = 2  # Individual-based: Permute genes separately
# perm.opt = 3  # Population-based: Permute genes separately across populations
#
# Author:: Pierre Legendre, April 2013, edited Thomas Lamy November 2017
{
eps <- 1e-10
### Internal functions
###  
gen.D <- function(Y, method)
# Compute genetic distance.
# Functions genind2genpop() and dist.genpop() are in {adegenet}
{
	# Convert from GENIND to GENPOP object (for D computed using "dist.genpop")
	# Population-level column sums of frequencies are multiplied by 2
	Dat <- genind2genpop(Y)
	# The population-level allele raw frequencies (6x93) are in Dat@tab

	# Compute genetic distance using dist.genpop() {adegenet}
	# dist.genpop() converts population-level raw freq. to relative frequencies
	# before computing the selected distance function
	Dist <- dist.genpop(Dat, method = method)

	list(Dist=Dist, Dat=Dat)
	}
###
centre <- function(D,n)
	# Centre a square matrix D by matrix algebra
	# mat.cen = (I - 11'/n) D (I - 11'/n)
	{	One <- matrix(1,n,n)
		mat <- diag(n) - One/n
		mat.cen <- mat %*% D %*% mat
	}
### End internal functions

# ===== The main function begins here =====
	require(adegenet)
	require(ade4)
	# epsilon <- .Machine$double.eps

	res <- gen.D(Y, method=D.opt)
	D <- res$Dist
	if(sqrt.D) D <- sqrt(D)
	alleles.per.locus <- nAll(res$Dat)
	nloci <- length(nAll(res$Dat))
	nn <- nrow(Y@tab)        # Number of individual observations
	n <- nrow(res$Dat@tab)   # Number of populations
	if(!silent){
		cat("n.indiv =",nn,"  n.pop =",n,"  nloci =",nloci,"\n")
		cat("Number of individuals per population:")
		print(table(Y@pop))
		}
	#
	### Compute from-to for permutation tests
	to <- cumsum(alleles.per.locus)
	from <- c(1,to[1:(nloci-1)]+1)
	# cat("from =",from,"\n")
	# cat("to   =",to  ,"\n")
	
	### Compute Gower transformation. 
	# Legendre & Legendre (2012), eqs. 9.40 and 9.42
	G <- centre(-0.5*as.matrix(D^2), n)
	SStotal <- sum(D^2)/n      # eq. 8
	BDtotal <- SStotal/(n-1)
	LCBD <- diag(G)/SStotal
	# Check SStotal
	SStotal.alt <- sum(diag(G))
	if(abs(SStotal-SStotal.alt)>eps) stop("SStotal =",SStotal,"  SStotal.alt =",SStotal.alt)  # TL
	# if(round(SStotal,5)!=round(SStotal.alt,5)) stop("SStotal =",SStotal,"  SStotal.alt =",SStotal.alt)
	#
	# Permutation test for LCBD indices
	# perm.opt: There are 3 ways of permuting the data --
	# 1. Individual-based: Permute whole rows of individual data
	# 2. Individual-based: Permute genes separately, individual data
	# 3. Population-based: Permute genes separately, pop. data. Not recommended.
	#    Test not valid if the sampling units are not of the same size (same n).
	#    (OK for ecology where individuals of a species are considered similar,
	#    if the sampling units are of the same size.)
	# => Permutation as well as parametric tests assume that the data are 
	# permutable. For that, they must at least correspond to sampling units of 
	# the same size. In the present tests, that means populations with the same  
	# number of sampled individuals.
	if(nperm > 0) {
		nGE.L <- rep(1,n)
		
		if(perm.opt == 1) { # Individual-based: Permute whole rows of Y
			Y.perm <- Y
			for(iperm in 1:nperm) {
				Y.perm@tab <- Y@tab[sample(nn),]
				res.perm <- gen.D(Y.perm, method=D.opt)
				D.perm <- res.perm$Dist
				if(sqrt.D) D.perm <- sqrt(D.perm)
				G.perm <- centre(-0.5*as.matrix(D.perm^2), n)
				
				SStotal.perm <- sum(diag(G.perm))
				LCBD.perm <- diag(G.perm)/SStotal.perm

				# cat(LCBD.perm,"\n")
				ge <- which(LCBD.perm >= LCBD)
				nGE.L[ge] <- nGE.L[ge] + 1
				}
			p.LCBD <- nGE.L/(nperm+1)
		
		} else if(perm.opt == 2) { # Individual-based: Permute genes separately
			Y.perm <- Y
			for(iperm in 1:nperm) {
				temp <- rep(NA,nn)
				for(k in 1:nloci) {
					order <- sample(nn)
					temp = cbind(temp, Y@tab[order,from[k]:to[k]])
					}
				Y.perm@tab <- temp[,-1]
				# cat(dim(Y.perm@tab),"\n")
				res.perm <- gen.D(Y.perm, method=D.opt)
				D.perm <- res.perm$Dist
				if(sqrt.D) D.perm <- sqrt(D.perm)
				G.perm <- centre(-0.5*as.matrix(D.perm^2), n)
				
				SStotal.perm <- sum(diag(G.perm))
				LCBD.perm <- diag(G.perm)/SStotal.perm

				# cat(LCBD.perm,"\n")
				ge <- which(LCBD.perm >= LCBD)
				nGE.L[ge] <- nGE.L[ge] + 1
				}
			p.LCBD <- nGE.L/(nperm+1)
	
		} else { # Population-based: Permute genes separately across populations
			Dat.perm <- res$Dat
			Y.pop <- Dat.perm@tab
			for(iperm in 1:nperm) {
				temp <- rep(NA,n)
				for(k in 1:nloci) {
					order <- sample(n)
					temp = cbind(temp, Y.pop[order,from[k]:to[k]])
					}
				Dat.perm@tab <- temp[,-1]
				D.perm <- dist.genpop(Dat.perm, method=D.opt)
				if(sqrt.D) D.perm <- sqrt(D.perm)
				G.perm <- centre(-0.5*as.matrix(D.perm^2), n)
				
				SStotal.perm <- sum(diag(G.perm))
				LCBD.perm <- diag(G.perm)/SStotal.perm

				# cat(LCBD.perm,"\n")
				ge <- which(LCBD.perm >= LCBD)
				nGE.L[ge] <- nGE.L[ge] + 1
				}
			p.LCBD <- nGE.L/(nperm+1)
		} 
	} else { p.LCBD <- NA }
#
if(is.euclid(D)) { 
	note <- paste("D =",D.opt,"  Matrix D is Euclidean")
	} else {
	note <- paste("D =",D.opt,"  Matrix D is not Euclidean; sqrt(D) may be Euclidean")
	}
if(!save.D) D <- NA

list(BDtotal=BDtotal, SStotal=SStotal, LCBD=LCBD, p.LCBD=p.LCBD, D=D, note=note)
}


simul.gen.LCBD <- function(Y, nsim=100, nperm=999, D.opt=2, perm.opt=3, sqrt.D=FALSE)
# Simulations to estimate type I error rates 
# using file inputDat25 where 25 indiv. have been pre-selected per population:
#    n.indiv = 150   n.pop = 6   nloci = 8 
#    Number of individuals per population:
#    GP WK FC CA CN FV 
#    25 25 25 25 25 25
#
# Y : data in Fstat format
# nsim  : Number of simulations
# nperm : Number of permutations
# D.opt : Reference number of the distance coefficient in dist.genpop()
# perm.opt : Reference number of the permutation method
# sqrt.D=TRUE : Compute LCBD from the square root of the distances in D
#
# Author:: Pierre Legendre, April 2013
{
out.LCBD = matrix(NA,nsim,6)
out.p = matrix(NA,nsim,6)
YY <- Y

aa <- system.time({
for(i in 1:nsim) {
	YY@tab <- Y@tab[sample(150),]
	tmp <- genetic.LCBD(YY, D.opt=D.opt, perm.opt=perm.opt, sqrt.D=sqrt.D, nperm=nperm, save.D=FALSE, silent=TRUE)
	out.LCBD[i,] = tmp$LCBD
	out.p[i,] = tmp$p.LCBD
	}
})
aa[3] <- sprintf("%2f",aa[3])
cat("Time for simulation =",aa[3]," sec",'\n')

LCBD.means <- apply(out.LCBD, 2, mean)
rejection1 <- function(vec, alpha=0.20) length(which(vec <= alpha))
rejection2 <- function(vec, alpha=0.10) length(which(vec <= alpha))
rejection3 <- function(vec, alpha=0.05) length(which(vec <= alpha))
rejection4 <- function(vec, alpha=0.01) length(which(vec <= alpha))
rej.rates <- matrix(NA,4,length(Y@pop.names))
rownames(rej.rates) <- c("0.20", "0.10", "0.05", "0.01")
colnames(rej.rates) <- Y@pop.names
rej.rates[1,] <- apply(out.p, 2, rejection1)
rej.rates[2,] <- apply(out.p, 2, rejection2)
rej.rates[3,] <- apply(out.p, 2, rejection3)
rej.rates[4,] <- apply(out.p, 2, rejection4)

list(rej.rates=rej.rates/nsim, LCBD.means=LCBD.means, out.LCBD=out.LCBD, out.p=out.p)
}

simul.gen.LCBD2 <- function(Y, n.ind=25, nsim=100, nperm=999, D.opt=2, perm.opt=3, sqrt.D=FALSE)
# Simulations to estimate type I error rates 
# Resample file inputDat
#    n.indiv = 469   n.pop = 6   nloci = 8 
#    Number of individuals per population:
#     GP  WK  FC  CA  CN  FV 
#    157  69  41  97  27  78
#
# Y : data in Fstat format
# n.ind : How many individuals will be selected per population
# nsim  : Number of simulations
# nperm : Number of permutations
# D.opt : Reference number of the distance coefficient in dist.genpop()
# perm.opt : Reference number of the permutation method
# sqrt.D=TRUE : Compute LCBD from the square root of the distances in D
#
# Author:: Pierre Legendre, April 2013
{
### Internal function

sample.Fstat.file <- function(Y, n.ind)
# Sample at random 'n.ind' per population in a file with Fstat format
# Replace list elements @tab, @ind.names, @pop
{
	YY <- Y
	n.alleles <- sum(Y@loc.nall)
	nn <- nrow(Y@tab)        # Number of individual observations
	indiv.per.pop <- table(inputDat@pop)
	npop <- length(indiv.per.pop)
	to <- cumsum(indiv.per.pop)
	from <- c(1,to[1:(npop-1)]+1)
	Y.sel <- rep(NA, n.alleles)
	vec.ind.names <- NA
	vec.pop <- NA
	# vec.pop <- Y@pop[1]
	# print(class(vec.pop))
	# print(vec.pop)
	for(ii in 1:npop) {
		order <- sample(to[ii]-from[ii]+1)
		temp.data <- Y@tab[from[ii]:to[ii],]
		temp.ind.names <- Y@ind.names[from[ii]:to[ii]]
		temp.pop <- Y@pop[from[ii]:to[ii]]
		# print(temp.pop[order[1:n.ind]])
		Y.sel <- rbind(Y.sel, temp.data[order[1:n.ind],])
		vec.ind.names <- c(vec.ind.names, temp.ind.names[order[1:n.ind]])
		vec.pop <- as.factor(c(vec.pop, temp.pop[order[1:n.ind]]))
		# print(class(vec.pop))
		# print(vec.pop)
		}
	YY@tab <- Y.sel[-1,]
	YY@ind.names <- vec.ind.names[-1]
	# print(vec.pop)
	YY@pop <- as.factor(vec.pop[-1])
	YY
	}
#
### End internal function

# Create output matrices
out.LCBD = matrix(NA,nsim,6)
out.p = matrix(NA,nsim,6)

aa <- system.time({
for(i in 1:nsim) {
	YY <- sample.Fstat.file(Y, n.ind)
	tmp <- genetic.LCBD(YY, D.opt=D.opt, perm.opt=perm.opt, sqrt.D=sqrt.D, nperm=nperm, save.D=FALSE, silent=TRUE)
	out.LCBD[i,] = tmp$LCBD
	out.p[i,] = tmp$p.LCBD
	}
})
aa[3] <- sprintf("%2f",aa[3])
cat("Time for simulation =",aa[3]," sec",'\n')

LCBD.means <- apply(out.LCBD, 2, mean)
rejection1 <- function(vec, alpha=0.20) length(which(vec <= alpha))
rejection2 <- function(vec, alpha=0.10) length(which(vec <= alpha))
rejection3 <- function(vec, alpha=0.05) length(which(vec <= alpha))
rejection4 <- function(vec, alpha=0.01) length(which(vec <= alpha))
rej.rates <- matrix(NA,4,length(Y@pop.names))
rownames(rej.rates) <- c("0.20", "0.10", "0.05", "0.01")
colnames(rej.rates) <- Y@pop.names
rej.rates[1,] <- apply(out.p, 2, rejection1)
rej.rates[2,] <- apply(out.p, 2, rejection2)
rej.rates[3,] <- apply(out.p, 2, rejection3)
rej.rates[4,] <- apply(out.p, 2, rejection4)

list(rej.rates=rej.rates/nsim, LCBD.means=LCBD.means, out.LCBD=out.LCBD, out.p=out.p)
}