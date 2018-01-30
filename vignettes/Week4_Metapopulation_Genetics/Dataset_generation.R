##############################################################
# Assembly and export of data sets 'dd.ecogen' and 'dd.site' #
##############################################################

# Note: fixed some spatial coordinates in import file 'MEC_5478_sm_Data.csv'.


require(LandGenCourse)
require(EcoGenetics)
require(gstudio)
require(dplyr)

Snail.all.gstudio <- gstudio::read_population(
  "./vignettes/Week5_Metapopulation_Genetics/Lamy/MEC_5478_sm_Data.csv", 
  locus.columns = c(7:26), type="column")

# Import tables
Table2 <- read.csv('./vignettes/Week5_Metapopulation_Genetics/Lamy/Lamy_etal_MolEco2012_Table2.csv', header=TRUE)
Table4 <- read.csv('./vignettes/Week5_Metapopulation_Genetics/Lamy/Lamy_etal_MolEco2012_Table4.csv', header=TRUE)
Table2.strata <- Table2[,c("SiteID", "SITE", "YEAR", "Spatial", "MultiYear", "APE")]
Table2.diversity <- Table2[,c("SiteID", "n", "RA", "He", "f", "s")]
Table4.strata <- Table4[,c("SiteID", "Cluster")]
Table4.env <- Table4[,c("SiteID", "Type", "FST.GESTE", "Size", "V", "C", "Stab", "D", "APA", "NLT", "Fst.temp")]


# Add 'Cluster' variable:
Snail.all.gstudio <- dplyr::left_join(Snail.all.gstudio, Table4.strata)


# Import from gstudio to 'ecogen' object:
Snail.all.ecogen <- EcoGenetics::gstudio2ecogen(Snail.all.gstudio, ID = "INDIVIDUAL", 
                         lat = "North.WGS84", lon = "West.WGS84", 
                         struct =c("SiteID", "SITE", "YEAR", "Cluster"))


# Assemble site data:
tmp <- dplyr::left_join(Table2.strata, Table4.strata, by="SiteID")
tmp <- dplyr::left_join(tmp, Table2.diversity, by="SiteID")
tmp <- dplyr::left_join(tmp, Table4.env, by="SiteID")

# Create objects to export:
dd.ecogen <- Snail.all.ecogen
dd.site <- tmp
row.names(dd.site) <- dd.site$SiteID

require(dplyr)
coords <- data.frame(dd.ecogen@XY, SiteID=dd.ecogen@S$SiteID)
coords <- coords %>% 
  group_by(SiteID)%>%
  unique()
dd.site.xy <- dplyr::left_join(dd.site, coords, by="SiteID")
sp::coordinates(dd.site.xy) <- ~ Longitude + Latitude
sp::proj4string(dd.site.xy) <- sp::CRS("+proj=longlat +datum=WGS84")
class(dd.site.xy)
dd.site <- dd.site.xy


# Export as '.Rdata' files:
save(dd.ecogen, file = "dd.ecogen.RData", compress='xz')
save(dd.site, file = "dd.site.RData", compress='xz')
load("dd.ecogen.RData")