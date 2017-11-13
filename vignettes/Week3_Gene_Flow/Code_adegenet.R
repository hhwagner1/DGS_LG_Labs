# From adegenet vignette

data(microbov)
toto <- genind2genpop(microbov)
popNames(toto)

titi <- toto[1:3,]
popNames(titi)
nAll(titi)
tata <- titi[,loc=c(1,3)]
tata
nAll(tata)
locNames(titi)
hel5 <- titi[,loc="HEL5"]
hel5
locNames(hel5)

data(nancycats)
sepCats <- seploc(nancycats)
class(sepCats)
names(sepCats)
sepCats$fca45
identical(tab(sepCats$fca45), tab(nancycats[,loc="fca45"]))

data(microbov)
obj <- seppop(microbov)
class(obj)
names(obj)
obj$Borgou

obj <- lapply(obj,seploc)
names(obj)
class(obj$Borgou)
names(obj$Borgou)
obj$Borgou$INRA63

# Merging data sets
obj <- seppop(microbov)
names(obj)
newObj <- repool(obj$Borgou, obj$Charolais)
newObj
popNames(newObj)

# process other:
data(sim2pop)
sim2pop
nInd(sim2pop)
head(other(sim2pop)$xy)
dim(other(sim2pop)$xy)
other(genind2genpop(sim2pop, process.other=TRUE))

## Using Summaries

toto <- summary(nancycats)
names(toto)
par(mfrow=c(2,2))
plot(toto$pop.eff, toto$pop.nall, xlab="Colonies sample size",
     ylab="Number of alleles",main="Alleles numbers and sample sizes",
     type="n")
text(toto$pop.eff,toto$pop.nall,lab=names(toto$pop.eff))
barplot(toto$loc.n.all, ylab="Number of alleles",
        main="Number of alleles per locus")
barplot(toto$Hexp-toto$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
barplot(toto$pop.eff, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)

# Is mean observed H significantly lower than mean expected H ?
bartlett.test(list(toto$Hexp,toto$Hobs))
t.test(toto$Hexp,toto$Hobs,pair=T,var.equal=TRUE,alter="greater")

# HWE
library(pegas)
data(nancycats)
cats.hwt <- hw.test(nancycats, B=0)
cats.hwt

#  Measuring and testing population structure (a.k.a F statistics)
library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")
fstat(nancycats)

# by locus:
library(pegas)
Fst(as.loci(nancycats))

# testing:
Gtest <- gstat.randtest(nancycats,nsim=99)
Gtest
plot(Gtest)

# Pairwise fst
matFst <- pairwise.fst(nancycats[1:50,])
matFst
#The resulting matrix is Euclidean when there are no missing values:
is.euclid(matFst)

## Estimating inbreeding
data(microbov)
sal <- seppop(microbov)$Salers
sal
temp <- inbreeding(sal, N=100)
class(temp)
head(names(temp))
head(temp[[1]],20)
Fbar <- sapply(temp, mean)
hist(Fbar, col="firebrick", main="Average inbreeding in Salers cattles")
which(Fbar>0.4)
F <- inbreeding(sal, res.type="function")[which(Fbar>0.4)]
F
plot(F$FRBTSAL9266, main=paste("Inbreeding of individual",names(F)),
     xlab="Inbreeding (F)", ylab="Probability density")

## Isolation by distance

data(nancycats)
toto <- genind2genpop(nancycats)
Dgen <- dist.genpop(toto,method=2)
Dgeo <- dist(nancycats$other$xy)
ibd <- mantel.randtest(Dgen,Dgeo)
ibd
plot(ibd)



############################
# Adegenet and pegas vignette:

# Convert:

genind2loci

amova


#############################
# gstudio vignette:




################################
# spca vignette

plot(rupica.smry$Hexp, rupica.smry$Hobs, main = "Observed vs expected heterozygosity")
abline(0, 1, col = "red")


################################
# EcoGenetics 1.2.1
# CONTACT AUTHORS TO FIX IMPORT FUNCTIONS!

require(EcoGenetics)   # ADD TO LandGenCourse
outEco <- genind2ecogen(Frogs.genind)

gstudio2ecogen(Frogs.gstudio, lat = "Latitude", lon = "Longitude", ID = "ID", struct = "pop")
outEco <- gstudio2ecogen(Frogs.gstudio)


#########################
# Adegenet strata

# let's look at the microbov data set:
data(microbov)
microbov

# We see that we have three vectors of different names in the 'other' slot. 
# ?microbov
# These are Country, Breed, and Species
names(other(microbov))

# Let's set the strata
strata(microbov) <- data.frame(other(microbov))
microbov

# And change the names so we know what they are
nameStrata(microbov) <- ~Country/Breed/Species
