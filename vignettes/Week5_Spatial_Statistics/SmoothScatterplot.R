names(GD.pop)
i = 3

png(file="Scatterplot1.png", width=6, height=4, res=600, units="in")
plot(Dgeo, GD.pop[[i]], pch=20, cex=0.5, axes=TRUE, 
     xlab="Geographic Distance", ylab="Genetic Distance")
abline(lm(GD.pop[[i]] ~ Dgeo))
dev.off()


png(file="Scatterplot2.png", width=6, height=4, res=600, units="in")
plot(Dgeo, GD.pop[[i]], pch=20, cex=0.5, axes=TRUE, 
     xlab="Geographic Distance", ylab="Genetic Distance")
dens <- MASS::kde2d(Dgeo, GD.pop[[i]], n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(GD.pop[[i]] ~ Dgeo))
lines(loess.smooth(Dgeo, GD.pop[[i]]), col="red", lwd=2)
#lines(loess.smooth(Dgeo, GD.pop[[i]], span=1/4), col="red", 
#      lwd=1, lty=2)
dev.off()


png(file="Scatterplot3.png", width=6, height=4, res=600, units="in")
plot(Dgeo, GD.pop[[i]], pch=20, cex=0.5, axes=TRUE,log="x", 
     xlab="Geographic Distance", ylab="Genetic Distance")
dens <- MASS::kde2d(Dgeo, GD.pop[[i]], n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
image(dens, col=transp(myPal(300), 0.7), add=TRUE)
abline(lm(GD.pop[[i]] ~ Dgeo))
lines(loess.smooth(Dgeo, GD.pop[[i]]), col="red", lwd=2)
dev.off()



png(file="Scatterplot4.png", width=6, height=4, res=600, units="in")
plot(rank(Dgeo), rank(GD.pop[[i]]), pch=20, cex=0.5, axes=TRUE, 
     xlab="Rank of Geographic Distance", ylab="Rank of Genetic Distance")
dens <- MASS::kde2d(rank(Dgeo), rank(GD.pop[[i]]), n=300)
myPal <- colorRampPalette(c("white","blue","gold","orange","red"))
image(dens, col=transp(myPal(300), 0.7), add=TRUE)

abline(lm(rank(GD.pop[[i]]) ~ rank(Dgeo)))
lines(loess.smooth(rank(Dgeo), rank(GD.pop[[i]])), col="red", lwd=2)
dev.off()