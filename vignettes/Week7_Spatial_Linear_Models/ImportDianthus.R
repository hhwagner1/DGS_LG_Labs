######################
# Import dianthus data
######################

load("Dianthus.Rdata")












######################################
# Import Yessica's data:
# -------------------------------------------------
Dianthus_diversityPatchdata.csv
dist_sheepCte.txt
dist_sheepInt.txt
Patch_XY_Dianthus.csv

Patchdata <- read.csv("~/Desktop/R_GitHub_Projects/DGS_LG_Labs/vignettes/Week7_Spatial_Linear_Models/Dianthus_diversityPatchdata.csv")

PatchXY <- read.csv("~/Desktop/R_GitHub_Projects/DGS_LG_Labs/vignettes/Week7_Spatial_Linear_Models/Patch_XY_Dianthus.csv")

Dist_sheepCte <- read.delim("~/Desktop/R_GitHub_Projects/DGS_LG_Labs/vignettes/Week7_Spatial_Linear_Models/dist_sheepCte.txt")

Dist_sheepInt <- read.delim("~/Desktop/R_GitHub_Projects/DGS_LG_Labs/vignettes/Week7_Spatial_Linear_Models/dist_sheepInt.txt")


# Extract corresponding sites:
# ----------------------------
PatchXY$patch
PatchXY$id

Patchdata$patch

a <- a2 <- a3 <- rep(NA, nrow(Patchdata))
for(i in 1:nrow(Patchdata))
{
  a[i] <- grep(as.character(Patchdata$patch[i]), as.character(PatchXY$patch))
  a2[i] <- grep(as.character(Patchdata$patch[i]), names(Dist_sheepCte))  
  a3[i] <- grep(as.character(Patchdata$patch[i]), names(Dist_sheepInt))  
}
# just checking:
cbind(as.character(Patchdata$patch), as.character(PatchXY$patch)[a])
cbind(names(Dist_sheepCte[a2]), Dianthus$patch)
cbind(names(Dist_sheepInt[a3]), Dianthus$patch)
cbind(a, a2, a3)

require(dplyr)
Dianthus <- dplyr::left_join(Patchdata, PatchXY, by=c("patch"))
Dianthus_DsheepC <- Dist_sheepCte[a2,a2]
Dianthus_DsheepI <- Dist_sheepInt[a3,a3]

#####################

# Regression model:
# -----------------
cor(Dianthus$Ar, Dianthus[,sapply(Dianthus, is.numeric)], 
    use = "pairwise.complete.obs")
mod.lm.IBR <- lm(Dianthus.ecopop.subset@C$RA ~ Elements2008.4 + S.sheep.continuous, 
                 data = Dianthus.ecopop.subset@E)
summary(mod.lm.IBR)

mod.lm.XY <- lm(Dianthus.ecopop.subset@C$RA ~ X + Y +
                  Elements2008.4 + S.sheep.continuous, 
                data = cbind(Dianthus.ecopop.subset@XY, Dianthus.ecopop.subset@E))
summary(mod.lm.XY)

mod.lm.IBD <- lm(Dianthus.ecopop.subset@C$RA ~ Elements2008.4 + S.euclidean, 
                 data = Dianthus.ecopop.subset@E)
summary(mod.lm.IBD)

# Test for spatial autocorrelation (Moran's I):
# ---------------------------------------------
coords <- data.matrix(Dianthus.ecopop.subset@XY[,3:4])
nb.gab <- spdep::graph2nb(spdep::gabrielneigh(coords), sym=TRUE)
listw.gab <- spdep::nb2listw(nb.gab)

spdep::moran.test(Dianthus.ecopop.subset@C$RA, listw.gab)             # Y
spdep::moran.test(Dianthus.ecopop.subset@E$Elements2008.4, listw.gab)
spdep::moran.test(Dianthus.ecopop.subset@E$S.euclidean, listw.gab) 
spdep::moran.test(Dianthus.ecopop.subset@E$S.sheep.continuous, listw.gab) 

spdep::lm.morantest(mod.lm.IBR, listw.gab)          # residuals
spdep::lm.morantest(mod.lm.XY, listw.gab)       # residuals|(XY)
spdep::lm.morantest(mod.lm.IBD, listw.gab)       # residuals|(XY)



# Autocorrelation in environmental data:
t(rbind(sapply(Dianthus.ecopop.subset@E, function(ls) 
  spdep::moran.test(ls, listw.gab)$estimate),
  p.value = sapply(Dianthus.ecopop.subset@E, function(ls) 
    spdep::moran.test(ls, listw.gab)$p.value)))


###################################################
# Import data (version with "-9" replaced by "NA"):
# -------------------------------------------------

GenData <- read.csv("~/Desktop/R_GitHub_Projects/DGS_LG_Labs/vignettes/Week7_Spatial_Linear_Models/Dianthus_carthusianorum_genodata_Dec15_2013_NA.csv")

SiteData <- read.csv("~/Desktop/R_GitHub_Projects/DGS_LG_Labs/vignettes/Week7_Spatial_Linear_Models/SummaryData_12Sept2011.csv")

# Convert from two-column to separated format:
# --------------------------------------------
tmp <- tmp2 <- GenData[,6:27]
tmp2[] <- NA
for(i in seq(1,21,2))
{
  #tmp2[,i] <- paste(tmp[,i], tmp[,i+1], sep=":")
  tmp2[,i] <- paste(sprintf("%03d", tmp[,i]), 
                    sprintf("%03d", tmp[,i+1]), sep=":")
}
tmp2 <- tmp2[,seq(1,21,2)]
tmp2[tmp2 == "0NA:0NA"] <- "NA:NA"

# Import into ecogen object:
# --------------------------
Dianthus.ecogen <- EcoGenetics::ecogen(XY = GenData[,c("LAT", "LONG")], 
                   G = tmp2, S = GenData[,c(2:3)], missing="NA", 
                   set.names = GenData[,1], sep=":", type="codominant")
Dianthus.ecogen

# Convert lat-lon coordinates:
# ----------------------------
Dianthus.ecogen@XY <- cbind(Dianthus.ecogen@XY, 
                            SoDA::geoXY(latitude=Dianthus.ecogen@XY$LAT, 
                                        longitude=Dianthus.ecogen@XY$LONG))


# Convert to ecopop object
Dianthus.ecopop <- EcoGenetics::ecogen2ecopop(Dianthus.ecogen, hier="PatchID")
Dianthus.ecopop@S

a4 <- 




# Convert to genind object and calculate allelic richness (RA):
# -------------------------------------------------------------
Dianthus.genind <- EcoGenetics::ecogen2genind(Dianthus.ecogen)
Dianthus.genind@pop <- Dianthus.ecogen@S$PatchID
RA <- PopGenReport::allel.rich(Dianthus.genind)$mean.richness

# Copy RA to ecogen object:
# -------------------------
Dianthus.ecogen <- eco.fill_ecogen_with_df(Dianthus.ecogen, pop="PatchID",
                                             pop_levels = names(RA), 
                                             C=data.frame(RA=RA))

# Subset individuals from patches with site data:
# -----------------------------------------------

a <- rep(NA, nrow(Dianthus.ecogen@S))
for(i in 1:nrow(Dianthus.ecogen@S))
{
    a[i] <- grep(Dianthus.ecogen@S$PatchID[i], SiteData$Label)  
}

Dianthus.ecogen.subset <- Dianthus.ecogen[!is.na(a),]

# Copy site data into ecogen slot E
# ---------------------------------

tmp3 <- SiteData[a[!is.na(a)],]
row.names(tmp3) <- row.names(Dianthus.ecogen.subset@G)
EcoGenetics::ecoslot.E(Dianthus.ecogen.subset) <- tmp3
Dianthus.ecogen.subset

# Aggregate to ecopop:
# --------------------
Dianthus.ecopop.subset <- EcoGenetics::ecogen2ecopop(Dianthus.ecogen.subset,
                                                     hier="PatchID")
Dianthus.ecopop.subset@E <- aue.aggregated_df(Dianthus.ecogen.subset@E, 
                            hier=Dianthus.ecogen.subset@S[,"PatchID"], 
                            fun=function(x) mean(x,na.rm = TRUE), 
                            factor_to_counts = FALSE)
Dianthus.ecopop.subset@E
Dianthus.ecopop.subset@E <- Dianthus.ecopop.subset@E[,!is.na(colSums(Dianthus.ecopop.subset@E))]
Dianthus.ecopop.subset



# Regression model:
# -----------------
cor(Dianthus.ecopop.subset@C$RA, Dianthus.ecopop.subset@E)
mod.lm.IBR <- lm(Dianthus.ecopop.subset@C$RA ~ Elements2008.4 + S.sheep.continuous, 
             data = Dianthus.ecopop.subset@E)
summary(mod.lm.IBR)

mod.lm.XY <- lm(Dianthus.ecopop.subset@C$RA ~ X + Y +
                Elements2008.4 + S.sheep.continuous, 
                data = cbind(Dianthus.ecopop.subset@XY, Dianthus.ecopop.subset@E))
summary(mod.lm.XY)

mod.lm.IBD <- lm(Dianthus.ecopop.subset@C$RA ~ Elements2008.4 + S.euclidean, 
             data = Dianthus.ecopop.subset@E)
summary(mod.lm.IBD)

# Test for spatial autocorrelation (Moran's I):
# ---------------------------------------------
coords <- data.matrix(Dianthus.ecopop.subset@XY[,3:4])
nb.gab <- spdep::graph2nb(spdep::gabrielneigh(coords), sym=TRUE)
listw.gab <- spdep::nb2listw(nb.gab)

spdep::moran.test(Dianthus.ecopop.subset@C$RA, listw.gab)             # Y
spdep::moran.test(Dianthus.ecopop.subset@E$Elements2008.4, listw.gab)
spdep::moran.test(Dianthus.ecopop.subset@E$S.euclidean, listw.gab) 
spdep::moran.test(Dianthus.ecopop.subset@E$S.sheep.continuous, listw.gab) 

spdep::lm.morantest(mod.lm.IBR, listw.gab)          # residuals
spdep::lm.morantest(mod.lm.XY, listw.gab)       # residuals|(XY)
spdep::lm.morantest(mod.lm.IBD, listw.gab)       # residuals|(XY)



# Autocorrelation in environmental data:
t(rbind(sapply(Dianthus.ecopop.subset@E, function(ls) 
  spdep::moran.test(ls, listw.gab)$estimate),
        p.value = sapply(Dianthus.ecopop.subset@E, function(ls) 
          spdep::moran.test(ls, listw.gab)$p.value)))