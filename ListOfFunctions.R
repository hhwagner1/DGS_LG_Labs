

Vignettes <- paste0("Week", c(1:14), "_vignette.R")
Vignettes <- c(Vignettes, "Week0_BasicR.R", "Week2_bonus_vignette.R", "Week8_bonus_vignette.R")
Vignettes_labels <- paste0("W", c(1:14))
Vignettes_labels <- c(Vignettes_labels, "W0", "W2b", "W8b")

Functions <- lapply(Vignettes, function(x) unlist(NCmisc::list.functions.in.file(
  system.file("doc", x, package = "LandGenCourse"))))

test <- unlist(Functions)
names(test) <- rep(Vignettes_labels, times=sapply(Functions, length))

test2 <- table(test, names(test))
#rio::export(test2, "./output/Index_of_functions.xlsx")
write.csv(test2, "./output/Index_of_functions.csv")


