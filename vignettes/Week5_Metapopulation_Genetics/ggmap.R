# Map with graph

png(file="Snails_map.png", width=8, height=7, unit="in", res=300)
data(dd.site, package="LandGenCourse")
Snails.spatial <- dd.site[dd.site@data$Spatial==TRUE,]
coords <- as.data.frame(Snails.spatial@coords)

map <- ggmap::get_map(location = Snails.spatial@bbox, maptype = "terrain", source = "google", zoom = 11)

ggmap::ggmap(map) + 
  ggplot2::geom_point(ggplot2::aes(x = Longitude, y = Latitude, 
                                   color = Snails.spatial@data$RA), 
                      size = 3, data=coords) + 
  ggplot2::scale_colour_gradientn(
    colors = colorRampPalette(c("blue","grey", "red"))(10), 
    name="Allelic Richness\n(rarefied)\n") + 
  ggplot2::annotate("text", x=coords$Longitude + 0.018, 
                    y=coords$Latitude, 
                    label = Snails.spatial@data$SITE, 
                    size=2.5, col="black", fontface="bold") +
  ggplot2::labs(title = substitute(paste("Genetic Diversity of ", italic('D. depressissimum'))), x = "", y = "")  
# + ggplot2::theme(legend.position="none")
dev.off()


