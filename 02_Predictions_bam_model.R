
rm(list = ls())

require(data.table)
require(tidyverse)
require(lubridate)
require(reshape2)
require(raster)
require(magrittr)
require(ncdf4)
require(rgdal)
require(matchingR)
require(mgcv)

# Directory tree
myDir <- "~/Documents/GitLab/dengue/"
output <- "/Volumes/LSHTM/dengueSEA/"

# ----------------------------------------
# Extract data for variable from
# netcdf files
# ----------------------------------------
den_data <- fread(file.path(output, "data",
                            "dengue_analysis_data_SEA_lshtm.csv"))

# Create date variable with monthly steps at day 15 of each month
den_data$date  <- ymd(den_data$date)

# Map
mymap <- rgdal::readOGR(file.path(myDir, "shapefiles"),
                        "sea_map_simplified")

# Load model
mod0 <- readRDS(file.path(myDir, "bam_model_sea.rds"))

# Define models
models <- c("gfdl", "ipsl", "mpi", "mri", "uksem")
# models <- c("uksem")

# Define scenarios
rcps <- c('ssp126', 'ssp370', 'ssp585')

"%ni%" <- Negate("%in%")
nsim  <- 1e3
betas <- rmvn(nsim, coef(mod0), mod0$Vp)
phi   <- mod0$family$getTheta(TRUE)         # dispersion parameter

# Linear spline function
rhs <- function(x, threshold) ifelse(x >= threshold, x - threshold, 0)

for(model in models){
   cat("\n", model, "\n")
   
   for(rcp in rcps){
      cat("\n", rcp, "\n")
      
      futureData <- fread(
         file.path(output, "data",
                   paste0("prediction_set_", model, "_",
                          rcp, "_sea.csv"))) %>%
         dplyr::select(-flux_grav) %>%
         dplyr::mutate(gdp2=gdp/1e6,
                       admin1=factor(admin1),
                       country=as.factor(country),
                       cluster=as.factor(cluster),
                       time=as.numeric(as.factor(date))+252,
                       lpop=log(population),
                       ID.area=factor(as.numeric(as.factor(admin1))),
                       ID.country=factor(as.numeric(as.factor(country))),
                       ID.month=as.numeric(as.character(month)),
                       ID.year=year - min(year) + 1,
                       ID.year=as.factor(ID.year+21))  %>%
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
                                             "> 1500"))) %>%
         dplyr::filter(year >= 2018)
      
      futureData$flux_rad[futureData$country=="singapore"] <- 
         quantile(futureData$flux_rad, prob=0.95, na.rm=T) 
      
      X <- predict(mod0, 
                   exclude = "s(ID.year)",
                   newdata = futureData,
                   type="lpmatrix") 
      
      futureData$pred_freq <- predict(mod0,
                                      exclude="s(ID.year)",
                                      type="response",
                                      newdata=futureData)
      
      mean_nb <- exp( X %*% t(betas) + futureData$lpop) #  mean of Neg Bin (need to manually add offset)
      preds   <- apply(mean_nb,1,function(x){rnbinom(length(x), mu=x, size=phi)})
      
      futureData$meanPreds <- apply(preds,2,mean)
      futureData$upperConf <- apply(preds,2,quantile,probs=0.975, na.rm=T)
      futureData$lowerConf <- apply(preds,2,quantile,probs=0.025, na.rm=T)
      
      futureData %<>%
         dplyr::filter(year >= 2018) %>%
         dplyr::mutate(cir=meanPreds/population*1e5,
                       cirlb=lowerConf/population*1e5,
                       cirub=upperConf/population*1e5,
                       lts=ifelse(cir > 10, 1, 0),
                       rcp=rcp,
                       model=model)
      
      annualData <- futureData %>%
         group_by(admin1, country, year, ssp) %>%
         dplyr::summarise(meanPreds=sum(meanPreds),
                          lowerConf=sum(lowerConf),
                          upperConf=sum(upperConf),
                          lts=sum(lts),
                          population=max(population),
                          popdens=max(popdens),
                          tas02=mean(tas02),
                          dry02=mean(dry02),
                          pre02=mean(pre02)) %>%
         dplyr::mutate(rcp=rcp,
                       model=model)
      
      fwrite(futureData, file.path(output, "output",
                                   paste0("predictions_", model, "_",
                                          rcp, "_sea.csv")))
      
      fwrite(annualData, file.path(output, "output",
                                   paste0("annual_estimates_", model, "_",
                                          rcp, "_sea.csv")))
   }
}


# End of file


