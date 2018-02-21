
meig <- spmoran::meigen(coords=xy)

SF <- data.frame(xy, meig$sf)
sp::coordinates(SF) <- ~ x + y
sp::spplot(SF, "X1", colorkey = TRUE)
sp::spplot.grid(SF)

ggplot(as.data.frame(Result), aes(x, y, size=PatchSize)) +
  geom_point(color="darkblue") + coord_fixed()


SF <- data.frame(xy, meigW$sf)
sp::coordinates(SF) <- ~ x + y
sp::spplot(SF, names(SF), colorkey = TRUE)

# a)
meig <- spmoran::meigen(coords=xy)
e_res <- spmoran::esf( y=Dianthus.df$A, x=Dianthus.df[,c("PatchSize", "IBR")],
                       meig=meig, fn = "r2" )
SF <- data.frame(xy, meig$sf, e_res$sf)
sp::coordinates(SF) <- ~ x + y
sp::spplot(SF, "e_res.sf", colorkey = TRUE)
cor(SF$e_res.sf, Dianthus.df$A)^2

# b)
cmat.d1    <- spdep::listw2mat( listw.d1)
meigW  <- spmoran::meigen( cmat = cmat.d1 )
e_res<- spmoran::esf( y=Dianthus.df$A, x=Dianthus.df[,c("PatchSize", "IBR")],
                       meig=meigW, fn = "r2" )

SF <- data.frame(xy, meigW$sf, e_res$sf)
sp::coordinates(SF) <- ~ x + y
sp::spplot(SF, "e_res.sf", colorkey = TRUE)
cor(SF$e_res.sf, Dianthus.df$A)^2
