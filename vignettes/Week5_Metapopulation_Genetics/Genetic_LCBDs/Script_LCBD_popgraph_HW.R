
require(popgraph)
require(igraph)
require(gstudio)

source('./vignettes/Week5_Habitat_Connectivity/Genetic_LCBDs/popgraph_HW.R')
source('./vignettes/Week5_Habitat_Connectivity/Genetic_LCBDs/dist_amova_HW.R')


graph.arapat <- popgraph::popgraph(x = to_mv(arapat), groups = arapat$Population)
graph.arapat.HW <- popgraph_HW(x = to_mv(arapat), groups = arapat$Population)
graph.arapat.amova.HW <- popgraph_HW(x = dist_amova_HW(arapat)$data, groups = arapat$Population)

V(graph.arapat.HW)$LCBD
sum(V(graph.arapat.HW)$LCBD)
V(graph.arapat.amova.HW)$LCBD

genetic_distance(arapat, mode="amova")
graph.arapat.HW <- popgraph_HW(x = to_mv(arapat), groups = arapat$Population)

graph.arapat.HW.genind <- popgraph_HW(arapat.genind@tab, groups = arapat.genind@pop)
V(graph.arapat.HW.genind)$LCBD
sum(V(graph.arapat.HW.genind)$LCBD)
cbind(V(graph.arapat.HW)$LCBD, V(graph.arapat.HW.genind)$LCBD)

data_ind <- import2genind("./vignettes/Week5_Habitat_Connectivity/Genetic_LCBDs/Data_25pops.gen", ncode=3)
write.table(genind2df(data_ind, sep=":"), "data_df.csv", sep=",")
data_ind.gstudio <- gstudio::read_population("data_df.csv", type="separated", locus.columns=c(2:11))


graph.data_ind <- popgraph::popgraph(gstudio::to_mv(data_ind.gstudio), groups = data_ind.gstudio$pop)
graph.data_ind.HW <- popgraph_HW(gstudio::to_mv(data_ind.gstudio), groups = data_ind.gstudio$pop)
V(graph.data_ind.HW)$LCBD
sum(V(graph.data_ind.HW)$LCBD)

### NEED TO SORT BY NAME!

LCBD.amova <- data.frame(Site = V(graph.data_ind.HW)$name, LCBD.amova = V(graph.data_ind.HW)$LCBD)
LCBD.orig <- data.frame(Site = names(table(data_ind@pop)), LCBD.orig = res.LCBD$LCBD)

LCBD.merged <- merge(LCBD.orig, LCBD.amova)
cor(LCBD.merged$LCBD.amova, LCBD.merged$LCBD.orig)


graph.data_ind.amova.HW <- popgraph_HW(x = dist_amova_HW(data_ind.gstudio)$data, 
                                       groups = data_ind.gstudio$pop)
V(graph.data_ind.amova.HW)$LCBD


# Compare with GESTE

Table4 <- read.csv("./vignettes/Week5_Habitat_Connectivity/Lamy/Lamy_etal_MolEco2012_Table4.csv")

cor(Table4$FST.GESTE, res.LCBD$LCBD)



# Create pop-level distance matrix from amova distances:
x = dist_amova_HW(data_ind.gstudio)$data
groups = data_ind.gstudio$pop
D <- centroid_distance(x, groups)

allLD <- centroid_distance(LDValues, groups)   # Centroid distance...?
allSD <- centroid_variance(LDValues, groups)   # Centroid variance...?
D <- matrix(0, nrow = K, ncol = K)
for (i in seq(1, K)) for (j in seq(i, K)) {
  if (i != j) {
    p1 <- unlist(allLD[i, ])
    p2 <- unlist(allLD[j, ])
    D[i, j] <- D[j, i] <- sqrt(sum((p1 - p2)^2))
  }
}





















## Transform data from "genind" to "genpop"
data_pop <- genind2genpop(data_ind)
popNames(data_pop)
site_names <- c("Pico","Roc","Sen","Vee","Gef","Rej","Ecl","Cou","Fdr","Mah", "Pis","Pev","Blo","Pou","Des","Bam","Baz","Hen","Mam","Tit","Del","Pav","Stj", "Kan","Ptc")

## compute genetic LCBD based on genetic Hellinger distance
res.LCBD <- genetic.LCBD(data_ind, D.opt=2, perm.opt=3, nperm=999) # 1 min
# perm.opt = 3  # Population-based: Permute genes separately across populations
res.LCBD$BDtotal
# [1] 0.1971397
res.LCBD$SStotal
# [1] 4.731352
res.LCBD$LCBD



## Transform data from "genind" to "genpop"
data_pop <- genind2genpop(data_ind)
popNames(data_pop)
site_names <- c("Pico","Roc","Sen","Vee","Gef","Rej","Ecl","Cou","Fdr","Mah", "Pis","Pev","Blo","Pou","Des","Bam","Baz","Hen","Mam","Tit","Del","Pav","Stj", "Kan","Ptc")



## load packages
library(ade4)
library (adegenet)

require(popgraph)
require(igraph)

data(arapat)
summary(arapat)

x <- to_mv(arapat)
nodes <- arapat$Population
graph <- popgraph(x, nodes)
graph

V(graph)$name
V(graph)$size
E(graph)$weight



## conditional genetic distance:
graph <- popgraph::popgraph(x = to_mv(arapat), groups = arapat$Population)
ret <- popgraph::to_matrix(graph, mode = "shortest path")
ret


require(EcoGenetics)
arapat.ecogen <- gstudio2ecogen(arapat, ID="ID", struct = c("Species", "Cluster", "Population"))
arapat.genind <- ecogen2genind(arapat.ecogen)
arapat.genind@pop <- arapat.genind@strata$Population
arapat.genpop <- genind2genpop(arapat.genind)

# LCBD:
# -----
sqrt.D=FALSE

Y <- arapat.genind
Dat <- genind2genpop(Y)
Dist <- dist.genpop(Dat, method = 2)
res <- list(Dist=Dist, Dat=Dat)

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



plot(LCBD, V(graph)$size)
plot(LCBD, V(graph)$LCBD)
cor(LCBD, V(graph)$LCBD)


## Popgraph code:
# ---------------

x <- to_mv(arapat)
groups <- arapat$Population
alpha = 0.05
tol = 1e-04


groups <- factor(as.character(groups))
t <- table(groups)
if (any(t < 4)) {
  popnames <- paste(names(which(t < 4)), collapse = ", ")
  warning(paste("You have strata (", popnames, ") that have fewer than 4 individuals. This anlaysis needs to have a good estimate of within stratum variance."))
}

N <- length(groups)                    # Number of individuals
K <- length(levels(groups))            # Number of groups
critVal <- qchisq(1 - alpha, 1) 
EdgeStr <- matrix(0, K, K)             # K x K matrix of zero's
pcfit <- prcomp(x, retx = T)           # PCA, returning rotated variables
mv <- pcfit$x[, pcfit$sdev > tol]      ## sqrt(diag(cov(pcfit$x))) = pcfit$sdev
P <- ncol(mv)                          # Dropped 2 PCA axes
Pop.priors <- as.numeric(table(groups)/N)     # Proportion of all individuals in each group
K <- length(Pop.priors)
pop.means <- tapply(mv, list(rep(groups, P), col(mv)), mean)
                                       # K x P matrix: mean of each column by group
                                       # Group centroids
sigma.w <- sqrt(diag(var(mv - pop.means[groups, ])))
                                       # For each column in P, the sd within groups. 

scaling <- diag(1/sigma.w, , P)        # Same as diag(1/sigma.w)?? Diagonal matrix with 1/sigma.w
fac <- 1/(N - K)
X <- sqrt(fac) * (mv - pop.means[groups, ]) %*% scaling
                                      # For each individual, the scaled deviation from group means.  
X.s <- svd(X, nu = 0)                 # Eigen analysis (svd) or residuals?
rank <- sum(X.s$d > tol)              # (drop axes with zero variance)
if (rank < P) 
  warning(paste((P - rank), " variables are collinear and being dropped from the discriminant rotation.", 
                sep = ""))            # Get rid of any collinear variables in P
scaling <- scaling %*% X.s$v[, 1:rank] %*% diag(1/X.s$d[1:rank], , rank)
                                      # Modify scaling
mu <- colSums(Pop.priors %*% pop.means)
                                      # For each column in P, the weighted mean (global centroid)
X <- sqrt((N * Pop.priors)/(K - 1)) * scale(pop.means, center = mu, 
                                            scale = FALSE) %*% scaling
                                      # Express population means as deviation from global centroid,
                                      # multiplied by sqrt of group size/(nPops - 1)
X.s <- svd(X, nu = 0)                 # Eigen analysis of weighted group means?
rank <- sum(X.s$d > tol * X.s$d[1L])
scaling <- scaling %*% X.s$v[, 1L:rank]
means <- colMeans(pop.means)          # Means of original group means
LDValues <- scale(mv, center = means, scale = FALSE) %*% 
  scaling                             # For each individual, deviation from each group mean?
allLD <- centroid_distance(LDValues, groups)   # Centroid distance...?
allSD <- centroid_variance(LDValues, groups)   # Centroid variance...?
D <- matrix(0, nrow = K, ncol = K)
for (i in seq(1, K)) for (j in seq(i, K)) {
  if (i != j) {
    p1 <- unlist(allLD[i, ])
    p2 <- unlist(allLD[j, ])
    D[i, j] <- D[j, i] <- sqrt(sum((p1 - p2)^2))
  }
}
rownames(D) <- colnames(D) <- rownames(allLD)        # A symmetric distance matrix between populations
totMean <- mean(D)
colMean <- colMeans(D)
colMeanMatrix <- matrix(colMean, K, K, byrow = T)
rowMeanMatrix <- matrix(colMean, K, K, byrow = F)
C <- -0.5 * (D - colMeanMatrix - rowMeanMatrix + totMean)
                                                     # Why 0.5? Is sum(LCBD) always 0.5 then? 

  # HW: add LCBD code
  SStotal <- sum(D)/K
  BDtotal <- SStotal/(K-1)
  LCBD.2 <- diag(C)/SStotal
  LCBD.2 <- diag(C)/sum(diag(C))
  
  # # Original LCBD code:
  # n <- nrow(D)
  # D <- sqrt(D)
  # G <- centre(-0.5*as.matrix(D^2), n)
  # cor(as.dist(C), as.dist(G))                      # Correlation = 1! The two are actually identical.
  # SStotal <- sum(D^2)/n      # eq. 8
  # BDtotal <- SStotal/(n-1)
  # LCBD <- diag(G)/SStotal
  # cbind(LCBD, LCBD.2)                              # Identical
  # 
  
  CI <- MASS::ginv(C)

R <- matrix(1, K, K)
for (i in 1:K) for (j in 1:K) if (i != j) 
  R[i, j] = C[i, j]/sqrt(C[i, i] * C[j, j])
SRI <- matrix(1, K, K)
RI <- MASS::ginv(R)
EED <- matrix(0, K, K)
for (i in seq(1, K)) for (j in seq(1, K)) if (i != j) 
  SRI[i, j] <- -1 * RI[i, j]/sqrt(RI[i, i] * RI[j, j])
SRI <- 1 - SRI^2
SRI[SRI < 0] <- 0
EED <- -N * log(SRI)
D[EED <= critVal] <- 0
graph <- graph.adjacency(D, mode = "undirected", weighted = TRUE, 
                         diag = FALSE)
igraph::V(graph)$name <- row.names(D)
igraph::V(graph)$LCBD <- LCBD.2
popSD <- scale(allSD, center = min(allSD), scale = TRUE) * 
  5 + 5
igraph::V(graph)$size <- popSD
class(graph) <- c("igraph", "popgraph")
return(graph)














# Popgraph:
# --------

l <- layout.fruchterman.reingold(graph)
plot(graph,layout=l,vertex.label=V(graph)$name)

t <- table(arapat$Species,arapat$Population)
t

colors <- rep("lightblue", length(V(graph)))
mainland_pops <- which( t[2,] > 0 )
colors[ mainland_pops ] <- "red"
cape_pops <- which( t[1,] > 0 )
colors[ cape_pops ] <- "green"
mixed_pops <- which( colSums(t>0) > 1 )
colors[ mixed_pops ] <- "yellow"
plot(graph,layout=l,vertex.label=V(graph)$name,vertex.color=colors)

baja <- arapat[arapat$Species=="Peninsula",]
inds.per.pop <- lapply( partition(baja,"Population"), function(x) dim(x)[1] )
## Examine inds per pop to figure out which have <5 individuals save in smPops
smPops <- c("Const","ESan","157","73","Aqu","Mat","98","75")
baja <- baja[ !(baja$Population %in% smPops) , ]
x <- to_mv(baja)
nodes <- as.character( baja$Population )
graph <- popgraph( x, nodes)

graph

l <- layout.fruchterman.reingold(graph)
plot(graph,layout=l,vertex.label=V(graph)$name)









source('~/Desktop/R_GitHub_Projects/DGS_LG_Labs/vignettes/Week5_Habitat_Connectivity/Genetic_LCBDs/genetic.LCBD.R')


data_ind <- import2genind("./vignettes/Week5_Habitat_Connectivity/Genetic_LCBDs/Data_25pops.gen", ncode=3)

## Transform data from "genind" to "genpop"
data_pop <- genind2genpop(data_ind)
popNames(data_pop)
site_names <- c("Pico","Roc","Sen","Vee","Gef","Rej","Ecl","Cou","Fdr","Mah", "Pis","Pev","Blo","Pou","Des","Bam","Baz","Hen","Mam","Tit","Del","Pav","Stj", "Kan","Ptc")


## compute genetic LCBD based on genetic Hellinger distance
res.LCBD <- genetic.LCBD(data_ind, D.opt=2, perm.opt=3, nperm=999) # 1 min
# perm.opt = 3  # Population-based: Permute genes separately across populations
res.LCBD$BDtotal
# [1] 0.1971397
res.LCBD$SStotal
# [1] 4.731352
res.LCBD$LCBD
#  [1] 0.03841087 0.04135809 0.03337402 0.04042985 0.03642595 0.02364188
#  [7] 0.03712992 0.04280018 0.04247896 0.03528766 0.04308931 0.04870993
# [13] 0.04327150 0.03100420 0.06482796 0.04320822 0.03101483 0.05763325
# [19] 0.03164599 0.03380485 0.04208690 0.02577452 0.03091257 0.03137882
# [25] 0.07029978
res.LCBD$p.LCBD
#  [1] 0.628 0.376 0.930 0.458 0.781 1.000 0.706 0.258 0.280 0.845 0.253
# [12] 0.034 0.223 0.986 0.001 0.215 0.989 0.001 0.983 0.935 0.293 1.000
# [23] 0.988 0.982 0.001
p.adjust(res.LCBD$p.LCBD, "holm")
#  [1] 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000
# [12] 0.748 1.000 1.000 0.025 1.000 1.000 0.025 1.000 1.000 1.000 1.000
# [23] 1.000 1.000 0.025
site_names[res.LCBD$p.LCBD < 0.05]
site_names[p.adjust(res.LCBD$p.LCBD, "holm") < 0.05]
## Significant genetic LCDB indices: sites (12), 15, 18, 25
## Pev, Des, Hen, Ptc
## Porte Enfer Vigie, Desbonnes, L'Henriette, Pointe des Ch?teaux


## Plot the LCBD values on maps of the Guadeloupe snail sites
require(SoDA)
env_25pop <- read.table('./vignettes/Week5_Habitat_Connectivity/Genetic_LCBDs/env_25pops.txt', header=TRUE, row.names=1)
XY.cart <- geoXY(env_25pop[,2], env_25pop[,3], unit=1000)
rownames(XY.cart) <- site_names
plot(XY.cart, asp=1, type="n", xlab="x coordinates (km)", ylab="y coordinates (km)", main="Map of LCBD, Guadeloupe snail", xlim=c(-5,40), ylim=c(-5, 35))
points(XY.cart, pch=21, col="white", bg="brown", cex=80*res.LCBD$LCBD)
text(XY.cart[1:24,], labels=site_names[1:24], pos=4)
text(XY.cart[25,1], XY.cart[25,2],labels=site_names[25], pos=2)


## Analysis against site characteristics
LCDB <- res.LCBD$LCBD
mtot <- lm(LCDB ~ ., data=env_25pop[,c(4:7)]) # only 4 environmental variables 
res.back <-  step(mtot, direction="backward") 
summary(res.back) 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.062217   0.007691   8.089 4.91e-08 ***
# Size         -0.012807   0.005117  -2.503  0.02023 *  
# Connectivity -0.007005   0.002179  -3.215  0.00399 ** 
# ---
# Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1

# Residual standard error: 0.009045 on 22 degrees of freedom
# Multiple R-squared:  0.3917,	Adjusted R-squared:  0.3364 
# F-statistic: 7.083 on 2 and 22 DF,  p-value: 0.00422
