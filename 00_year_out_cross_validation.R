
rm(list = ls())

# ----------------------------------------
# Load packages
# ----------------------------------------
require(data.table)
require(tidyverse)
require(lubridate)
require(ModelMetrics) # for mae() and mase()
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

bks <- quantile(den_data$gdp/1e9, probs=seq(0,1, by =0.2))

den_data$qgdp <- cut(den_data$gdp/1e9, breaks=bks, include.lowest=TRUE)

# Create date variable with monthly steps at day 15 of each month
den_data$date  <- ymd(den_data$date)

# Map
mymap   <- rgdal::readOGR(file.path(myDir, "shapefiles"),
                          "sea_map_simplified")

mymap@data$admin1 <- as.character(mymap@data$id)
mapNb             <- spdep::poly2nb(mymap)
names(mapNb)      <- mymap$id

#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
   dataz <- as.data.frame(dataz)
   #vif part
   form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
   dataz   <- data.frame(fooy=1,dataz)
   lm_mod  <- lm(form,dataz)
   
   cat("\n\nVariance inflation factors\n\n")
   print(myvif(lm_mod))
}
#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
   v <- vcov(mod)
   assign <- attributes(model.matrix(mod))$assign
   if (names(coefficients(mod)[1]) == "(Intercept)") {
      v <- v[-1, -1]
      assign <- assign[-1]
   } else warning("No intercept: vifs may not be sensible.")
   terms <- labels(terms(mod))
   n.terms <- length(terms)
   if (n.terms < 2) stop("The model contains fewer than 2 terms")
   if (length(assign) > dim(v)[1] ) {
      diag(tmp_cor)<-0
      if (any(tmp_cor==1.0)){
         return("Sample size is too small, 100% collinearity is present")
      } else {
         return("Sample size is too small")
      }
   }
   R <- cov2cor(v)
   detR <- det(R)
   result <- matrix(0, n.terms, 3)
   rownames(result) <- terms
   colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
   for (term in 1:n.terms) {
      subs <- which(assign == term)
      result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
      result[term, 2] <- length(subs)
   }
   if (all(result[, 2] == 1)) {
      result <- data.frame(GVIF=result[, 1])
   } else {
      result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
   }
   invisible(result)
}
#END VIF FUNCTIONS


# --------------------------------
# Data set 
# --------------------------------

# Load observation data
myData <- den_data %>%
   dplyr::select(admin1, country, date, month, year, total,
                 tas02, pre02, hur02, wnd02, dry02, wet02, 
                 sum02, population, popdens, gdp, total_arrivals, 
                 flux_rad, flux_grav, cluster, urban, qgdp) %>%
   dplyr::mutate(admin1=factor(admin1),
                 country=as.factor(country),
                 cluster=as.factor(cluster),
                 time=as.numeric(as.factor(date)),
                 lpop=log(population),
                 ID.area=factor(as.numeric(as.factor(admin1))),
                 ID.country=factor(as.numeric(as.factor(country))),
                 ID.month=as.numeric(as.character(month)),
                 ID.year=year - min(year) + 1,
                 ID.year=as.factor(ID.year)) %>%
   drop_na(pre02) %>%
   data.table()

myData$flux_rad[is.na(myData$flux_rad)] <- quantile(myData$flux_rad,
                                                    prob=0.95,
                                                    na.rm=T)
corData <- myData %>%
   dplyr::select( tas02, pre02, hur02, 
                  wnd02, dry02, wet02, 
                  sum02, popdens, gdp, 
                  total_arrivals, 
                  flux_rad, flux_grav)
corvif(corData)

# --------------------
# Base model
# --------------------

# Linear spline function
rhs <- function(x, threshold) ifelse(x >= threshold, x - threshold, 0)

years  <- 2000:2017

myData$tas02b <- rhs(myData$tas02, 22)
myData$tas02c <- rhs(myData$tas02, 29.5)

myData$dry02b <- rhs(myData$dry02, 7)
myData$dry02c <- rhs(myData$dry02, 22)

myData$gdp2 <- myData$gdp/1e6

fx <- total ~ 
   offset(lpop) +
   s(admin1, bs='mrf', xt=list(nb=mapNb)) +
   s(ID.area, bs="re") +
   s(ID.month, by=ID.country, bs="cr") +
   s(ID.year, bs="re") +
   tas02 + rhs(tas02, 22) + rhs(tas02, 29.5) +
   dry02 + rhs(dry02,7) + rhs(dry02,22) +
   urban +
   gdp2 +
   log(flux_rad) +
   log(total_arrivals) 


fy <- total ~ 
   offset(lpop) +
   s(admin1, bs='mrf', xt=list(nb=mapNb)) +
   s(ID.area, bs="re") +
   s(ID.month, by=ID.country, bs="cr") +
   s(ID.year, bs="re") +
   tas02 + rhs(tas02, 22) + rhs(tas02, 29.5) +
   dry02 + rhs(dry02,7) + rhs(dry02,22) +
   urban + 
   log(gdp2) +
   log(flux_rad) +
   log(total_arrivals) 

equations <- c(fx, fy)

myRmse1 <- myRmse2 <- vector("list")

for(i in seq_along(years)){
   cat("\n", i, "\n")
   
   for(j in seq_along(equations)){
      cat("\t", j)
      
      try({
         train_set <- filter(myData, year!=years[i]) %>%
            drop_na(total)
         test_set  <- filter(myData, year==years[i])
         
         mod0 <- bam(formula=as.formula(equations[[j]]), 
                     data=train_set, 
                     na.action=na.exclude, 
                     family=nb(),
                     gamma=2,
                     rho=0.8,
                     drop.unused.levels=FALSE,
                     select=TRUE,
                     method="fREML")
         
         test_set$pred0 <- predict.gam(mod0,
                                       newdata=test_set,
                                       type='response')
         
         test_set %<>% drop_na(total) %>% drop_na(pred0)
         
         myRmse1[[j]] <- data.table(equation=j,
                                    year=i,
                                    rmse=ModelMetrics::rmse(test_set$total, 
                                                            test_set$pred0))
         
      })
      
   }
   myRmse2[[i]] <- rbindlist(myRmse1)
}


Rmse1 <- rbindlist(myRmse2) %>%
   group_by(equation) %>%
   dplyr::summarise(rmse=mean(rmse)) %>%
   arrange(rmse)

Rmse1
