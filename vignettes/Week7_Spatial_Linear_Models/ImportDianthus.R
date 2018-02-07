######################
# Import dianthus data
######################


GenData <- read.csv()


SiteData <- read.csv("C:/Users/Helene/Desktop/DGS_LG_Labs/vignettes/Week7_Spatial_Linear_Models/SummaryData_12Sept2011.csv")

coords <- data.matrix(SiteData[,c("x", "y")])
a <- c(1:nrow(SiteData))[is.na(coords[,1])]
a2 <-c(1:nrow(SiteData))[is.na(SiteData$Richness2008)]
a <- c(a, a2)

mod.lm <- lm(Richness2008 ~ Area.ha + S.sheep.intermittent, data = SiteData[-a,])
summary(mod.lm)

mod.lm.XY <- lm(Richness2008 ~ Area.ha + S.sheep.intermittent + x + y, data = SiteData[-a,])
summary(mod.lm.XY)


nb.gab <- spdep::graph2nb(spdep::gabrielneigh(coords[-a,]), sym=TRUE)
listw.gab <- spdep::nb2listw(nb.gab)
spdep::moran.test(SiteData$Richness2008[-a], listw.gab)             # Y
spdep::lm.morantest(mod.lm, listw.gab)          # residuals
spdep::lm.morantest(mod.lm.XY, listw.gab)       # residuals|(XY)

# Strong autocorrelation in environmental data:
t(rbind(sapply(WWP.ecogen@E, function(ls) spdep::moran.test(ls, listw.gab)$estimate),
        p.value = sapply(WWP.ecogen@E, function(ls) spdep::moran.test(ls, listw.gab)$p.value)))