x <- Frogs.gstudio$A

ploidy: (N - N.missing) * ploidy / N







snp <- read.delim(system.file("extdata", "WWP_SNP_genotypes.txt", 
                              package = "LandGenCourse"), sep = "\t", header = T)
env <- read.delim(system.file("extdata", "WWP_environmental_data.txt", 
                              package = "LandGenCourse"),sep = "\t", header = T)
row.names(snp) <- snp$family 
row.names(env) <- env$family
WWP <- EcoGenetics::ecogen(XY = env[,3:4], G = snp[,-c(1:2)], 
                           E = env[,-c(1:4)], S = env[,1:2], order.G = FALSE)
WWP
WWP.genind <- ecogen2genind(WWP)
WWP.genind
WWP2 <- genind2ecogen(WWP.genind)




data(eco.test)
eco.genind <- ecogen2genind(eco)
eco2 <- genind2ecogen(eco.genind)

test <- EcoGenetics::ecogen(XY = env[,3:4], P = trait, G = snp, 
                            E = env[,-c(1:4)], S = env[,1:2], order.G = FALSE)

test.genind@other$latlong <- test@XY
test.genind@other$site <- test@E
test.genind@other$trait <- test@P
test.genind
test2 <- genind2ecogen(test.genind)

eco@A[1,1:2] <- NA

test1 <- EcoGenetics::ecogen(XY = test.genind@other$latlong, 
                             P = test.genind@other$trait, 
                             G = eco.convert(test.genind@tab, "alleles.matrix", 
                                             "matrix", ploidy = 2), 
                             E = test.genind@other$site, 
                             S = test.genind@strata, order.G = FALSE)

tmp <- eco.convert(test@A, "alleles.matrix", "matrix", ploidy = 2)
loc2al <- eco.convert(tmp, "matrix", "alleles.matrix", ploidy = 2)








# Example with a single character
ex1 <- c("AA","AC","AC","CC")
ex2 <- c("GG","GT","GT","TT")
ex <- c(sample(ex1, 20, rep= T), sample(ex2, 20, rep= T)) 
ex <- matrix(ex, 10, 4)
colnames(ex) <- letters[1:4]
rownames(ex) <- LETTERS[1:10]
# example data
ex  
recoded <- eco.format(ex, ploidy = 2, nout = 1,  recode = "all", show.codes = TRUE) 
# recoded data 
recoded

G.num <- sapply(data.frame(WWP@G), function(ls) eco.format(data.matrix(ls), ploidy = 2, nout = 1,  recode = "all", show.codes = TRUE)$out)




eco.format(data.frame(A=ex[,1]), ploidy = 2, nout = 1,  recode = "all", show.codes = TRUE)
 ex.df <- data.frame(ex)
 eco.format(ex.df[[1]], ploidy = 2, nout = 1,  recode = "all", show.codes = TRUE)
 sapply(data.frame(ex), function(ls) eco.format(ls, ploidy = 2, nout = 1,  recode = "all", show.codes = TRUE)$out)
 
eco.format(WWP0@G, ploidy = 2, nout = 1,  recode = "all", show.codes = TRUE)$out

ecogen2hierfstat(WWP, pop="population")
