

## load packages
library(ade4)
library (adegenet)


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
