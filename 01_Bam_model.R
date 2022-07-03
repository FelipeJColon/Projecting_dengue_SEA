
rm(list = ls())

# ----------------------------------------
# Load packages
# ----------------------------------------
require(data.table)
require(tidyverse)
require(lubridate)
require(ModelMetrics)
require(rgdal)
require(magrittr)
require(mgcv)
require(readr)
require(ncdf4)
require(raster)
require(spdep)

myDir <- "~/Documents/GitLab/dengue/"
output <- "/Volumes/LSHTM/dengueSEA/"

den_data <- fread(file.path(output, "data",
                            "dengue_analysis_data_SEA_lshtm.csv")) %>%
   dplyr::mutate(urban=ifelse(popdens <= 100 , 0,
                              ifelse(popdens > 100 & popdens <=300, 100,
                                     ifelse(popdens > 300 & popdens <=1500 , 
                                            300, 
                                            ifelse(popdens > 1500 , 1500, 
                                                   "none")))),
                 urban=factor(urban, levels=c(0,100,300,1500),
                              labels=c("<= 100",
                                       "101-300",
                                       "301-1500",
                                       "> 1500")))

# Create date variable with monthly steps at day 15 of each month
den_data$date  <- ymd(den_data$date)

# Map
myMap   <- rgdal::readOGR(file.path(myDir, "shapefiles"),
                          "sea_map_simplified")
mapNb        <- spdep::poly2nb(myMap)
names(mapNb) <- myMap$id

# --------------------------------
# Data set 
# --------------------------------

# Load observation data
myData <- den_data %>%
   dplyr::filter(year >=2000) %>%
   dplyr::select(admin1, country, date, month, year, ssp, total,
                 tas02, dry02, population, popdens, gdp, 
                 total_arrivals, flux_rad, flux_grav, cluster, 
                 urban) %>%
   dplyr::mutate(gdp2=gdp/1e6,
                 ssp="historical",
                 admin1=factor(admin1),
                 country=as.factor(country),
                 lpop=log(population),
                 ID.area=factor(as.numeric(as.factor(admin1))),
                 ID.country=factor(as.numeric(as.factor(country))),
                 ID.month=as.numeric(as.character(month)),
                 ID.year=year - min(year) + 1,
                 ID.year=factor(ID.year)) %>%
   drop_na(tas02) %>%
   data.table()

myData$flux_rad[myData$country=="singapore"] <- quantile(myData$flux_rad,
                                                         prob=0.95,
                                                         na.rm=T)

rhs <- function(x, threshold){ ifelse(x >= threshold, x - threshold, 0)}

rmvn <- function(n,mu,sig) { ## MVN random deviates
   L <- mroot(sig)
   m <- ncol(L)
   t(mu + L%*%matrix(rnorm(m*n),m,n)) 
}

# --------------------
# Base model
# --------------------

f0 <- total ~ 
   offset(lpop) +
   s(admin1, bs='mrf', xt=list(nb=mapNb)) +
   s(ID.area, bs="re") +
   s(ID.month, by=ID.country, bs="cr") +
   s(ID.year, bs="re") +
   tas02 + rhs(tas02, 22) + rhs(tas02, 29.5) +
   dry02 + rhs(dry02, 7) + rhs(dry02, 22) +
   urban +
   gdp2 +
   log1p(flux_rad) +
   log1p(total_arrivals)

training <- myData %>% filter(ssp=="historical") %>%
   dplyr::mutate(ID.year=factor(ID.year),
                 ID.country=factor(ID.country),
                 ID.area=factor(ID.area))

model0 <- bam(formula=as.formula(f0), 
              data=training, 
              na.action=na.exclude, 
              family=nb(),
              gamma=2,
              rho=0.8,
              drop.unused.levels=FALSE)

saveRDS(model0, file.path(myDir, "bam_model_sea.rds"))

