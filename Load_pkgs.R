print(R.version.string)
version <- R.Version()
major   <- as.numeric(version$major)
minor   <- as.numeric(version$minor)

stopifnot(major>3 || (major==3 && minor>=2))

pkglist <- c("data.table", "tidyverse", "lubridate", "ModelMetrics", 
             "magrittr", "mgcv", "readr", "ncdf4", "raster", 
             "spdep", "reshape2", "rgdal")

new.packages <- pkglist[!(pkglist %in% installed.packages()[,"Package"])]

if(length(new.packages)){
  install.packages(new.packages,
                   repos="http://cran.ma.imperial.ac.uk/") 
}

lapply(pkglist, require, character.only=TRUE)

if(major>3 || (major==3 && minor>=2)==TRUE){
  message('R packages successfully loaded')
}

rm(pkglist, new.packages, version)


# ----------------
# Eof
# ----------------