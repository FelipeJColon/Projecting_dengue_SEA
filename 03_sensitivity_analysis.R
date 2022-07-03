
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
require(spdep)

myDir <- "~/Documents/GitLab/dengue/"
output <- "/Volumes/LSHTM/dengueSEA/"

# Map
mymap   <- rgdal::readOGR(file.path(myDir, "shapefiles"),
                          "sea_map_simplified")

mymap@data$admin1 <- as.character(mymap@data$id)
mapNb             <- spdep::poly2nb(mymap)
names(mapNb)      <- mymap$id

rhs <- function(x, threshold){ ifelse(x >= threshold, x - threshold, 0)}

mod0 <- readRDS(file.path(myDir, "bam_model_sea.rds"))

# historical data
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
                                       "> 1500")),
                 date=ymd(date),
                 gdp2=gdp/1e6)


# Provinces
provs <- den_data %>%
   dplyr::select(admin1) %>%
   distinct()

# Define models
models <- c("gfdl", "ipsl", "mpi", "mri", "uksem")

# Define scenarios
rcps <- c('ssp126', 'ssp370', 'ssp585')

for(model in models){
   cat("\n", model, "\n")
   
   for(rcp in rcps){
      cat("\n", rcp, "\n")
      
      futureData <- fread(
         file.path(output, "data",
                   paste0("prediction_set_", model, "_",
                          rcp, "_sea.csv"))) %>%
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
                                             "> 1500")))
      
      futureData$flux_rad[futureData$country=="singapore"] <- 
         quantile(futureData$flux_rad, prob=0.95, na.rm=T)
      
      # Climate constant at their mean for baseline period
      # All other variables change as per SSPs
      fd_clim <- futureData %>%
         dplyr::filter(year >= 2020) %>%
         dplyr::select(admin1, country, date, month,
                       popdens, lpop, gdp2, flux_rad,
                       total_arrivals, urban,
                       ID.area, ID.month, ID.country, ID.year) %>%
         dplyr::mutate(admin1=as.character(admin1))
      
      h_clim <- den_data %>%
         dplyr::select(admin1, month, tas02, dry02) %>%
         group_by(admin1, month) %>%
         dplyr::summarise(tas02=mean(tas02),
                          dry02=mean(dry02)) %>%
         data.table() 
      
      fd_clim %<>% left_join(h_clim)
      
      # Pop density constant at their mean for baseline period
      # All other variables change as per SSPs
      fd_pdn <- futureData %>%
         dplyr::filter(year >= 2020) %>%
         dplyr::select(admin1, country, date, month,
                       tas02, dry02, gdp2, lpop,
                       flux_rad, total_arrivals,
                       ID.area, ID.month, ID.country, 
                       ID.year) %>%
         dplyr::mutate(admin1=as.character(admin1))
      
      # Set popdens to mean values
      h_pden <- den_data %>%
         dplyr::select(admin1, month, popdens, lpop) %>%
         group_by(admin1, month) %>%
         dplyr::summarise(popdens=mean(popdens)) %>%
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
         data.table() 
      
      fd_pdn %<>% left_join(h_pden)
      
      # Total Pop constant at their mean for baseline period
      # All other variables change as per SSPs
      fd_pop <- futureData %>%
         dplyr::filter(year >= 2020) %>%
         dplyr::select(admin1, country, date, month,
                       tas02, dry02, gdp2, popdens,
                       urban, flux_rad, total_arrivals,
                       ID.area, ID.month, ID.country, 
                       ID.year) %>%
         dplyr::mutate(admin1=as.character(admin1))
      
      # Set total pop to mean values
      h_pop <- den_data %>%
         dplyr::select(admin1, month, popdens, lpop) %>%
         group_by(admin1, month) %>%
         dplyr::summarise(lpop=mean(lpop)) %>%
         data.table() 
      
      fd_pop %<>% left_join(h_pop)
      
      # GDP constant at their mean for baseline period
      # All other variables change as per SSPs
      fd_gdp <- futureData %>%
         dplyr::filter(year >= 2020) %>%
         dplyr::select(admin1, country, date, month,
                       tas02, dry02, lpop, popdens,
                       urban, flux_rad, total_arrivals,
                       ID.area, ID.month, ID.country, ID.year) %>%
         dplyr::mutate(admin1=as.character(admin1))
      
      h_gdp <- den_data %>%
         dplyr::select(admin1, month, gdp2) %>%
         group_by(admin1, month) %>%
         dplyr::summarise(gdp2=mean(gdp2)) %>%
         data.table() 
      
      fd_gdp %<>% left_join(h_gdp)
      
      # Human mobility constant at their mean for baseline period
      # All other variables change as per SSPs
      fd_mob <- futureData %>%
         dplyr::filter(year >= 2020) %>%
         dplyr::select(admin1, country, date, month,
                       tas02, dry02, lpop, popdens,
                       urban, gdp2, ID.area, ID.month,
                       ID.country, ID.year) %>%
         dplyr::mutate(admin1=as.character(admin1))
      
      h_mob <- den_data %>%
         dplyr::select(admin1, month, flux_rad, total_arrivals) %>%
         group_by(admin1, month) %>%
         dplyr::summarise(flux_rad=mean(flux_rad),
                          total_arrivals=mean(total_arrivals)) %>%
         data.table() 
      
      fd_mob %<>% left_join(h_mob)
      
      p_clim <- predict(mod0, type='response', newdata=fd_clim, se.fit=TRUE)
      p_pden <- predict(mod0, type='response', newdata=fd_pdn, se.fit=TRUE)
      p_pop  <- predict(mod0, type='response', newdata=fd_pop, se.fit=TRUE)
      p_gdp  <- predict(mod0, type='response', newdata=fd_gdp, se.fit=TRUE)
      p_mob  <- predict(mod0, type='response', newdata=fd_mob, se.fit=TRUE)
      
      fd_clim$meanPreds <- p_clim$fit
      fd_pdn$meanPreds  <- p_pden$fit
      fd_pop$meanPreds  <- p_pop$fit
      fd_gdp$meanPreds  <- p_gdp$fit
      fd_mob$meanPreds  <- p_mob$fit
      
      fd_clim$upperConf <- p_clim$fit + (1.96 * p_clim$se.fit)
      fd_pdn$upperConf  <- p_pden$fit + (1.96 * p_pden$se.fit)
      fd_pop$upperConf  <- p_pop$fit + (1.96 * p_pop$se.fit)
      fd_gdp$upperConf  <- p_gdp$fit + (1.96 * p_gdp$se.fit)
      fd_mob$upperConf  <- p_mob$fit + (1.96 * p_mob$se.fit)
      
      fd_clim$lowerConf <- p_clim$fit - (1.96 * p_clim$se.fit)
      fd_pdn$lowerConf  <- p_pden$fit - (1.96 * p_pden$se.fit)
      fd_pop$lowerConf  <- p_pop$fit - (1.96 * p_pop$se.fit)
      fd_gdp$lowerConf  <- p_gdp$fit - (1.96 * p_gdp$se.fit)
      fd_mob$lowerConf  <- p_mob$fit - (1.96 * p_mob$se.fit)
      
      fd_clim %<>% dplyr::mutate(model=model, rcp=rcp)
      fd_pdn %<>% dplyr::mutate(model=model, rcp=rcp)
      fd_pop %<>% dplyr::mutate(model=model, rcp=rcp)
      fd_gdp %<>% dplyr::mutate(model=model, rcp=rcp)
      fd_mob %<>% dplyr::mutate(model=model, rcp=rcp)
      
      fwrite(fd_clim, file.path(output, "output", "mean_vals",
                                paste0("predictions_climate_", model, "_",
                                       rcp, "_sea.csv")))
      fwrite(fd_pdn, file.path(output, "output", "mean_vals",
                               paste0("predictions_popdens_", model, "_",
                                      rcp, "_sea.csv")))
      fwrite(fd_pop, file.path(output, "output", "mean_vals",
                               paste0("predictions_population_", model, "_",
                                      rcp, "_sea.csv")))
      fwrite(fd_gdp, file.path(output, "output", "mean_vals",
                               paste0("predictions_gdp_", model, "_",
                                      rcp, "_sea.csv")))
      fwrite(fd_mob, file.path(output, "output", "mean_vals",
                               paste0("predictions_human_mobility_", model, "_",
                                      rcp, "_sea.csv")))
   }
}


# EoF


