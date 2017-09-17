###
### This script is a modification of the last part of DGS - 2016
### week 6 lab. It is not a stand-alone. The beginning of the lab assignment
### are needed to run this part but where removed from this file for the
### sake of clarity.
###

## The issue in the original lab script was that the y axis was not consistent across the
## sub-plots, (as it was mentioned in the instructions).
## This is due to the fact that we want to plot an object of class "mantel.correlog"
## that has a default plot method function that ignore any unimplemented argument
## such as ylim. I was unable to find this function so I wrote my own.
## Note that it overrides the default function and is called automatically when plot()
## is called on a "mantel.correlog" object and this might not be desirable outside
## the framework of this lab.

## (Re)wrighting a ploting method function
plot.mantel.correlog = function(mc, xlim=NULL, ylim=NULL, alpha=0.05){
  # Extracting the correlogram
  mct = as.data.frame(mc$mantel.res)
  mct <- mct[!is.na(mct$Mantel.cor), ]
  # Computing plot limits if needed
  if(is.null(xlim)) xlim=range(pretty(mct$class.index))
  if(is.null(xlim)) xlim=range(pretty(mct$Mantel.cor))
  # Actual plot, matel.correlog style
  plot(mct$Mantel.cor ~ mct$class.index, pch=22, bg=mct$"Pr(corrected)"<=alpha, ylim=ylim, xlim=xlim, type="o")
  abline(h=0, col="red")
  
}

## The following is modified from the instructions (very last section):

# Table of number of unique pairs per distance class:
nUniquePairs.drop1 <- Reduce(cbind, lapply(Mantel.cg.drop1, 
                                           function(x) x$mantel.res[,2]/2))
dimnames(nUniquePairs.drop1)[[2]] <- names(Mantel.cg.drop1)
par(mfrow=c(3,3), mai=c(0.5,0.5,0.1,0.1))
# computing the range of distance on the full datatset, because
# it might varries when basins holding the closest points are dropped
Xlim = range(pretty(Mantel.cg.drop1$AllBasins$mantel.res[
  !is.na(Mantel.cg.drop1$AllBasins$mantel.res)[, "Mantel.cor"], "class.index"]))
# The same won't work on the correlation value so we have to manually (and arbitrarily) set a constant range
Ylim=c(-0.3, 0.5)
# In the rest of the script, the only thing that changes is the addition of the Xlim and Ylim parameters.
plot(Mantel.cg.drop1$AllBasins, xlim=Xlim, ylim=Ylim)
mtext(text= "All basins ", side=3, line=-2, adj=1)
for(i in 1:(length(Mantel.cg.drop1)-1))
{
  plot(Mantel.cg.drop1[[i]], xlim=Xlim, ylim=Ylim)
  mtext(text=paste(names(Mantel.cg.drop1)[i], "dropped ", sep=" "), 
        side=3, line=-2, adj=1)
}