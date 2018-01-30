
library(LandGenCourse)
data(dd.site)

# Extract 'spatial' subset
dd.spatial <- dd.site[dd.site@data$Spatial==TRUE,]
dd.df <- as.data.frame(dd.spatial)
a <- is.element(row.names(dd.spatial), c("32", "42"))

# Plot FST.GESTE against NLP, C
# -----------------------------
png(filename="FSTvsX.png", width = 7, height = 3.5, units = "in", 
    res = 300, pointsize = 12)
par(mfrow=c(1,2), mar=c(4.5,4.5,1,1), bg=NA)

# First plot: NLT
  plot(FST.GESTE ~ NLT, data=dd.spatial, 
       xlab = "Long-term population size", ylab = "Site-level Fst",
       ylim=c(0, 0.3))
  abline(lm(FST.GESTE ~ NLT, data=dd.spatial))
  a <- is.element(row.names(dd.spatial), c("32", "42"))
  points(FST.GESTE ~ NLT, data=dd.spatial[a,], col="red", pch=16)
  with(dd.spatial@data[a,], 
       text(x = NLT, y = FST.GESTE, labels = SITE, pos=4))

  # Second plot: C
  plot(FST.GESTE ~ C, data=dd.spatial, 
       xlab = "Hydrological connectivity", ylab = "Site-level Fst",
       ylim=c(0, 0.3))
  abline(lm(FST.GESTE ~ C, data=dd.spatial))
  a <- is.element(row.names(dd.spatial), c("32", "42"))
  points(FST.GESTE ~ C, data=dd.spatial[a,], col="red", pch=16)
  with(dd.spatial@data[a,], 
       text(x = C, y = FST.GESTE, labels = SITE, pos=4))
dev.off()


# With ggplot2::qplot
# --------------------

library(ggplot2)
dd.df <- as.data.frame(dd.spatial)  # ggplot2 can't handle sp objects

#Outlier <- is.element(row.names(dd.spatial), c("32", "42"))

Plot1 = qplot(x=NLT, y=FST.GESTE, data=dd.df, ylim=c(0,0.3),  
      xlab = "Long-term population size", ylab = "Site-level Fst",
      geom = c("point", "smooth"), method="lm")

Plot2 = qplot(x=C, y=FST.GESTE, data=dd.df, ylim=c(0,0.3),  
          xlab = "Hydrological connectivity", ylab = "Site-level Fst",
          geom = c("point", "smooth"), method="lm")

ggsave(filename = "qplot.png", width = 7, height = 3.5, unit = "in", dpi = 300,
       plot = cowplot::plot_grid(Plot1, Plot2, labels = c("A", "B"))) 


# With ggplot2::ggplot
# --------------------

Plot1 <- ggplot(data = dd.df, aes(x = NLT, y = FST.GESTE)) + 
  xlab("Long-term population size") + ylab("Site-level Fst") +
  geom_point() + 
  geom_smooth(col = "blue", method="lm") +
  geom_point(data=dd.df[a,],color="red",size=3) +
  geom_text(data=dd.df[a,], mapping=aes(x=NLT, y=FST.GESTE, label=SITE), 
          size=4, hjust = 0, nudge_x=0.1) 

Plot2 <- ggplot(data = dd.df, aes(x = C, y = FST.GESTE)) + 
  xlab("Hydrological connectivity") + ylab("Site-level Fst") +
  geom_point() + geom_smooth(col = "blue", method="lm") +
  geom_point(data=dd.df[a,],color="red",size=3) +
  geom_text(data=dd.df[a,], mapping=aes(x=C, y=FST.GESTE, label=SITE), 
            size=4, hjust = 0, nudge_x=0.1) 

ggsave(filename = "ggplot.png", width = 7, height = 3.5, unit = "in", dpi = 300, 
       plot = cowplot::plot_grid(Plot1, Plot2, labels = c("A", "B"))) 




Plot1 <- ggplot(data = df, aes( x, y )) + 
  xlab("Label x") + ylab("Label y") +
  geom_point() + 
  geom_smooth(col = "blue", method="lm") +
  geom_point(data = df[a,],color = "red",size = 3) +
  geom_text(data = df[a,], mapping=aes(x, y, label = SITE), 
            size = 4, hjust = 0, nudge_x = 0.1) 



ggplot(data = dd.df, aes(x = NLT, y = FST.GESTE, size = C)) + 
  xlab("Long-term population size") + ylab("Site-level Fst") +
  geom_point() +

  Plot3 <- ggplot(data = dd.df, aes(x = NLT, y = FST.GESTE)) +
  xlab("Long-term population size") + ylab("Site-level Fst") +
  geom_point() +
  facet_grid(. ~ Cluster) +
  geom_smooth(col = "blue", method = "lm", se = FALSE) +
  geom_point(data = dd.df[a, ],
             color = "red",
             size = 3) +
  geom_text(
    data = dd.df[a, ],
    mapping = aes(x = NLT, y = FST.GESTE, label = SITE),
    size = 4,
    hjust = 0,
    nudge_x = 0.1
  )
ggsave(
  filename = "facets.png",
  width = 7,
  height = 3.5,
  unit = "in",
  dpi = 300,
  plot = Plot3
)



Outlier <- a
ggplot(data = dd.df, aes(x = NLT, y = FST.GESTE, colour = Outlier)) + 
  xlab("Long-term population size") + ylab("Site-level Fst") +
  geom_smooth(col = "blue", method="lm") + 
  geom_point(size=3) 




# ggmap
# -----




# Export residuals

Residuals <- residuals(lm(FST.GESTE ~ NLT + C, data=dd.spatial))
Fitted <- fitted(lm(FST.GESTE ~ NLT + C, data=dd.spatial))

# Plot values in space with sp
# ----------------------------
#dd.df <- data.frame(dd.spatial@coords, dd.spatial@data)
dd.spatial@data$Residuals <- residuals(lm(FST.GESTE ~ NLT + C, data=dd.spatial))
dd.spatial@data$Fitted <- fitted(lm(FST.GESTE ~ NLT + C, data=dd.spatial))

sp::spplot(dd.spatial, zcol="FST.GESTE", colorkey=TRUE)
sp::spplot(dd.spatial, zcol="NLT", colorkey=TRUE)
sp::spplot(dd.spatial, zcol="C", colorkey=TRUE)
sp::spplot(dd.spatial, zcol="Fitted", colorkey=TRUE, col.regions=cm.colors(7))
sp::bubble(dd.spatial, zcol="Residuals", col = c("#FF00FF", "#0000FF"))

#sp::panel.spplot(rows=2, cols=2, sp.layout=?)


# Map

ggmap::qmplot(Longitude, Latitude, data=dd.df, 
       color = sign(Residuals), size = abs(Residuals)) +
ggplot2::geom_text(data = dd.df[a,], 
                   mapping=aes(Longitude, Latitude, label = SITE), 
                   size = 4, vjust = 0, nudge_y = 0.-0.015, color="black")
ggplot2::ggsave(filename = "ResidualMap.png", width = 7, height = 5.5, 
       unit = "in", dpi = 300) 



Terrain.google <- qmplot(Longitude, Latitude, data=dd.df, 
                         maptype="terrain", source="google")
Satellite.google <- qmplot(Longitude, Latitude, data=dd.df, 
                           maptype="satellite", source="google")
Roadmap.google <- qmplot(Longitude, Latitude, data=dd.df, 
                         maptype="roadmap", source="google")
Hybrid.google <- qmplot(Longitude, Latitude, data=dd.df, 
                         maptype="hybrid", source="google")
Toner.lite.stamen <- qmplot(Longitude, Latitude, data=dd.df, 
                            maptype="toner-lite", source="stamen")

Watercolor.stamen <- qmplot(Longitude, Latitude, data=dd.df, 
                       maptype="watercolor", source="stamen")

cowplot::plot_grid(Terrain.google, Satellite.google, 
                   Roadmap.google, Hybrid.google,
                   Toner.lite.stamen, Watercolor.stamen, 
                   nrow=2, ncol=3, label_colour="white",
                   labels = c("Terrain (Google)", "Satellite (Google)",
                              "Roadmap (Google)", "Hybrid (Google)",
                              "Toner-lite (Stamen)", "Watercolor (Stamen)"))

ggsave(filename = "Maps6.png", width = 12*0.7, height = 8*0.6, unit = "in", dpi = 300,
       cowplot::plot_grid(Terrain.google, Satellite.google, 
                          Toner.lite.stamen, 
                          Roadmap.google, Hybrid.google,
                          Watercolor.stamen, hjust=-0.4,
                          nrow=2, ncol=3, label_colour="white", align="hv",
                          labels = c("Terrain (Google)", "Satellite (Google)","Toner-lite (Stamen)",
                                     "Roadmap (Google)", "Hybrid (Google)",
                                     "Watercolor (Stamen)"))) 
