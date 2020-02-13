#' @title Install DGS packages 
#' @description Installs require packages for DGS and LandGenCourse package
#'
#' @param    install.cran - (TRUE | FALSE) install default packages from CRAN
#' @param    install.git - (TRUE | FALSE) install default packages from GitHub
#' @param    install.versions - (TRUE | FALSE) install archived versions from CRAN
#' @param    rm.installed - (FALSE | TRUE) If installed, remove previous version of packages  
#' @param    new.packages - A vector of package names to be added to default packages
#' @param    new.gits - A vector of GItHub package names to be added to default packages
#'                        must follow "user/package" convention
#' @param    rm.packages - A vector of packages to remove from CRAN defaults
#' @param    repo - A CRAN mirror, default is "revolutionanalytics" 
#' @param    lib - A path to an alternate local library
#' @param    install.LandGenCourse - (FALSE | TRUE) Validate dependencies and install LandGenCourse
#'
#' @examples
#' #### Install CRAN, specific CRAN versions and GitHib
#'   install.DGS()
#'
#' #### Install CRAN, specific CRAN versions and GitHib,
#' ####   removing existing/installed packages first 
#'   install.DGS(rm.installed = TRUE)
#'
#' #### Install only GitHib packages (excluding LandGenCourse)
#'   install.DGS(install.cran = FALSE, install.versions = FALSE)
#'
#' #### This installs all dependencies, validate installs and
#' ####   then install the LandGenCourse from GitHub
#'   install.DGS(install.LandGenCourse = TRUE)
#'
#' #### This allows you to add a CRAN package to the default packages
#'   install.DGS(install.cran = TRUE, new.packages = c("spData")) 
#' 
#' #### This allows you to remove a CRAN package from the default packages
#'   install.DGS(install.cran = TRUE, rm.packages = c("base64enc")) 
#' 
#' @export
install.DGS <- function(install.cran = TRUE, install.git = TRUE, 
                        install.versions = TRUE, rm.installed = FALSE,
                        new.packages = NULL, new.gits = NULL,  
						rm.packages = NULL,repo = NULL, lib = NULL,
						install.LandGenCourse = FALSE) {
  #### set a CRAN mirror
  if(is.null(repo)){
    repo <- "https://cran.revolutionanalytics.com" 
  }
  local({r <- getOption("repos")
    r["CRAN"] <- repo
      options(repos=r)})
 
  ##### Set local library path to R install / library
  if(is.null(lib)){
    if( Sys.info()['sysname'] != "Windows") {
      lib <- file.path(chartr("\\", "/", R.home()), "library")
    } else {
      lib <- file.path(R.home(), "library")
    }
    .libPaths(lib)
  } else {
    .libPaths(lib)
  }
  
  # Specify CRAN packages and optionally append with newly requested packages or 
  #   remove specified packages that are no longer needed
  packages <- c("base64enc", "BiocManager", "car", "cowplot", "data.table", "doParallel",
                "dplyr", "EcoGenetics", "effsize", "formatR", "gdistance", "GeNetIt",
                "geosphere", "ggmap", "ggplot2", "here", "hierfstat", "httpuv", "igraph",
                "landscapemetrics", "lattice", "lme4", "mapplots", "maps", "maptools",
                "microbenchmark", "mmod", "MuMIn", "mvtnorm", "nlme", "pegas", "poppr",
                "predictmeans", "profvis", "proto", "purrr", "pwr", "RANN", "raster",
                "rasterVis", "RColorBrewer", "readr", "rgdal", "rgeos", "RgoogleMaps", 
                "rio", "rlang", "rmarkdown", "sampling", "seqinr", "sf", "SoDA", "sp",
                "spacetime", "spatialEco", "spatialreg", "spdep", "stringi", "Sunder",
                "testthat", "tibble", "usdm", "vegan") 		
  if(!is.null(rm.packages)){ 
    rm.idx <- which(packages %in% rm.packages)
	if(length(rm.idx) > 0) {
      packages <- packages[-rm.idx]
    }
  }	
  if(!is.null(new.packages)){ packages <- c(packages, new.packages) }
 
  ############################
  #### libraries function ####
  libraries <- function(x, add = FALSE, install = TRUE, check.source = TRUE,
                        repository = NULL, lib = NULL) {
  	if(is.null(repository)){
        if(getOption("repos")["CRAN"] == "@CRAN@")
          options(repos="https://cloud.r-project.org/") 
  	} else {
  	  mirrors <- read.csv(file.path(R.home(), "/doc/CRAN_mirrors.csv"))$URL
  	    if( length(grep(repository, mirrors)) < 1 )
  		  stop("Not a valid repository mirror")
  	  options(repos=repository)  
      }		  
  	if(check.source) {
  	  options(install.packages.check.source = "yes")
      } else {
  	  options(install.packages.check.source = "no")
      }
  	ap <- utils::available.packages()
        pmiss <- x[which(!x %in% rownames(ap))]
      if(length(pmiss) > 0) { 
  	  cat("\n\n", "The following packages are not available on this mirror:", "\n",  
                  paste(pmiss, collapse=", "), "\n\n")    
  	    xprint <- x[-which(x %in% pmiss)]
  	}
      if (is.null(lib)) {
        lib <- .libPaths()[1L]
        if(length(.libPaths()) > 1L) 
          message(sprintf(ngettext(length(pkgs), "Installing package into %s\n(as %s is unspecified)", 
                  "Installing packages into %s\n(as %s is unspecified)"), 
                  sQuote(lib), sQuote("lib")), domain = NA)
      }
      ok <- dir.exists(lib) & (file.access(lib, 2) == 0L)
      if (length(lib) > 1 && any(!ok)) 
          stop(sprintf(ngettext(sum(!ok), "'lib' element %s is not a writable directory", 
              "'lib' elements %s are not writable directories"), 
              paste(sQuote(lib[!ok]), collapse = ", ")), 
              domain = NA)
      if (length(lib) == 1L && .Platform$OS.type == "windows") {
          ok <- dir.exists(lib)
          if (ok) {
              fn <- file.path(lib, paste0("_test_dir_", Sys.getpid()))
              unlink(fn, recursive = TRUE)
              res <- try(dir.create(fn, showWarnings = FALSE))
              if (inherits(res, "try-error") || !res) 
                  ok <- FALSE
              else unlink(fn, recursive = TRUE)
          }
      }
      if (length(lib) == 1L && !ok) {
        warning(gettextf("'lib = \"%s\"' is not writable", lib), 
  		      domain = NA, immediate. = TRUE)
        userdir <- unlist(strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep))[1L]
      }
    lib.path <- normalizePath(lib)
    pkg <- as.data.frame(installed.packages())
    for(i in 1:length(x)) {
      if(!x[i] %in% pkg$Package) { 
        cat(x[i], " is not installed","\n\n")
        if(install) {
          cat("    Installing", x[i], "to", file.path(R.home(), "library"), "\n\n")	  
  	      try( install.packages(x[i], lib = lib.path) )
  		if(add) require(x[i], quietly = TRUE, character.only = TRUE)
        }	  
  	} else {
  	  cat(x[i], as.character(pkg[which(pkg$Package %in% x[i]),]$Version), "is installed", "\n\n")
  	    if(add) require(x[i], quietly = TRUE, character.only = TRUE) 
      }
    }
    if(add) cat("\n\n", "These packages have been added to the current environment: ", "\n",  
                paste(xprint, collapse=", "), "\n\n")   
  }
  
  pkg <- as.data.frame(installed.packages())$Package
  
  # checks and installs devtools/remotes 
  libraries(c("devtools","remotes"), check.source = FALSE, 
            repository = getOption("repos"), lib = .Library)

  # function for removing previous versions 
  #  for rm.installed = TRUE
  rm.pkg <- function(p) {
    ip <- as.data.frame(installed.packages())$Package
      for(i in p) {
        if(i %in% ip)	  
	      remove.packages(i, .Library) 
	  }
  }

  ################################################
  #### Install specific CRAN package versions ####
  if(install.versions == TRUE) {
    vp <- c("ade4", "adegenet",  "rstudioapi", "spmoran", 
            "shiny", "swirl", "swirlify", "formatR", "miniUI")
    r = getOption("repos")["CRAN"]
	  ver <- as.data.frame(installed.packages())$Version
        if(rm.installed) { rm.pkg(vp) }  
	if(!"ade4" %in% pkg || ver[which(pkg %in% "ade4")] != "1.7-8") {
      remotes::install_version("ade4", version = "1.7-8", repos = r, upgrade = "never")}
	if(!"adegenet" %in% pkg || ver[which(pkg %in% "adegenet")] != "2.1.0") {
      remotes::install_version("adegenet", version = "2.1.0", repos = r, upgrade = "never")}
	if(!"rstudioapi" %in% pkg || ver[which(pkg %in% "rstudioapi")] != "0.5") {
      remotes::install_version("rstudioapi", version = "0.5", repos = r, upgrade = "never")}
	if(!"spmoran" %in% pkg || ver[which(pkg %in% "spmoran")] != "0.1.7.2") {
      remotes::install_version("spmoran", version = "0.1.7.2", repos = r, upgrade = "never")}
	if(!"shiny" %in% pkg || ver[which(pkg %in% "shiny")] != "0.13") {
      remotes::install_version("shiny", version = "0.13", repos = r, upgrade = "never")}
	if(!"swirl" %in% pkg || ver[which(pkg %in% "swirl")] != "2.4.3") {
      remotes::install_version("swirl", version = "2.4.3", repos = r, upgrade = "never")}
	if(!"swirlify" %in% pkg || ver[which(pkg %in% "swirlify")] != "0.5.1") {
      remotes::install_version("swirlify", version = "0.5.1", repos = r, upgrade = "never")}
	if(!"formatR" %in% pkg || ver[which(pkg %in% "formatR")] != "1.5") {
      remotes::install_version("formatR", version = "1.5", repos = r, upgrade = "never")}
	if(!"miniUI" %in% pkg || ver[which(pkg %in% "miniUI")] != "0.1.1") {
      remotes::install_version("miniUI", version = "0.1.1", repos = r, upgrade = "never")}   
    packages <- unique(c(packages, "ade4", "adegenet",  "rstudioapi", "spmoran", 
                         "shiny", "swirl", "swirlify", "formatR", "miniUI"))
  }			

  #################################  
  #### Install GitHub packages ####
  if(install.git == TRUE) {
    gp <- c("kjgilbert/QstFstComp", "dyerlab/popgraph", "dyerlab/gstudio")
      if(!is.null(new.gits)){ gp <- c(gp, new.gits) }	
        if(rm.installed) { rm.pkg(basename(gp)) }
	  for(i in gp) {
        remotes::install_github(i, force=FALSE, upgrade = "never")
      }
    packages <- unique(c(packages, "QstFstComp", "popgraph", "gstudio")) 	
    } 
  	
  #############################################	
  #### Install CRAN packages from binaries ####  
  if(install.cran == TRUE) {
    #### Verify that the packages are available in the specified repository
    ap <- utils::available.packages()
      pmiss <- packages[which(!packages %in% rownames(ap))]
        if(length(pmiss) > 0)
    	  cat("\n\n", "The following packages are not available on this mirror:", "\n",  
            paste(pmiss, collapse=", "), "\n\n")  
    
    #### Install available CRAN packages
    if(rm.installed) { rm.pkg(packages) }  
      libraries(packages, check.source = FALSE, 
                repository = getOption("repos"), 
                lib = .Library)
  }
  
  #######################################################
  #### Check installs and then install LandGenCourse #### 
  if(install.LandGenCourse){
    pkg <- as.character(as.data.frame(installed.packages())$Package)
      pkg.idx <- pkg[which(pkg %in% packages)]
        if(length(pkg.idx) != length(packages)) {
          cat("The following packages are missing: ", packages[which(!packages %in% pkg.idx)], "\n\n") 
        } else {
          remotes::install_github("hhwagner1/LandGenCourse", force=TRUE, upgrade = "never")
        }
  }		
} # end function
