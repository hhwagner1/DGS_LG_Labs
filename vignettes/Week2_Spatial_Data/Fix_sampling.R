The result will be a nested `tibble` containing the `plot_id`, the metrics and a `RasterLayer` with the clipped buffers (as a `list`). Attention: Because several metrics on class- and landscape-level the clipped buffers will be "repeated" several times.

Here we show results for the first four sampling locations:
  
  ```{r, fig.width=8, fig.height=5.5}
unique_plots <- unique(nlcd_sampled_plots$raster_sample_plots)[1:4]

par(mfrow = c(2,2))
plot(unique_plots[[1]], 
     main = paste(Sites.sp$SiteName[1]), 
     col = rev(rainbow(9)))
plot(unique_plots[[2]],
     main = paste(Sites.sp$SiteName[2]),
     col = rev(rainbow(9)))
plot(unique_plots[[3]],
     main = paste(Sites.sp$SiteName[3]), 
     col = rev(rainbow(9)))
plot(unique_plots[[4]],
     main = paste(Sites.sp$SiteName[4]), 
     col = rev(rainbow(9)))
par(mfrow = c(1,1))
```
