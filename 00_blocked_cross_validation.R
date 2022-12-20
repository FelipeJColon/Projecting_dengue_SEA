
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
                flux_rad, flux_grav, cluster) %>%
  drop_na(pre02) %>%
  dplyr::mutate(gdp2=gdp/1e6,
                ssp="historical",
                admin1=factor(admin1),
                country=as.factor(country),
                time=as.numeric(as.factor(date)),
                lpop=log(population),
                ID.area=factor(as.numeric(as.factor(admin1))),
                ID.country=factor(as.numeric(as.factor(country))),
                ID.month=as.numeric(as.character(month)),
                ID.year=year - min(year) + 1,
                ID.year=factor(ID.year)) %>%
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

# Add variables to mix
myVars2 <- c("tas02 + rhs(tas02,22) + rhs(tas02,29.5)", 
             "dry02 + rhs(dry02,7) + rhs(dry02,22)", 
             "urban", "gdp2", 
             "log(total_arrivals)", "log(flux_rad)")

# Add possible combinations
nvar <- length(myVars2)

myCombs <- unlist(lapply(1:nvar,
                         function(x){combn(1:nvar, x, 
                                           simplify=FALSE)}),
                  recursive=FALSE)

# Baseline model
fn <- total ~ 
   offset(lpop) +
   s(admin1, bs='mrf', xt=list(nb=mapNb)) +
   s(ID.area, bs="re") +
   s(ID.month, by=ID.country, bs="cr") +
   s(ID.year, bs="re") 


# Generate all possible equations
equations2 <- sapply(myCombs, function(x){
  paste("total", "~", "offset(lpop) +",
        "s(admin1, bs='mrf', xt=list(nb=mapNb)) +",
        "s(ID.area, bs='re') +",
        "s(ID.month, by=ID.country, bs='cr') +",
        "s(ID.year, bs='re') +",
        paste(myVars2[x], collapse=" + "))
})

# Join all equations
formulae <- c(fn, equations2)



no_cores <- detectCores()

# Define cluster
cl <- makeCluster(no_cores)

registerDoParallel(detectCores())  

years <- 1:24
"%ni%" <- Negate("%in%")

# Export vars

x <- foreach(i=0:215, .combine='rbind', .inorder=FALSE) %:%
  # cat("\n", i, "\n")
  
  foreach(j=1:length(formulae), .combine='rbind', .inorder=FALSE) %dopar% {
    cat("\t", j)
    
    try({
      
      my_years  <- years + i
      train_set <- filter(myData, time%ni%my_years)
      test_set  <- filter(myData, time%in%my_years)
      
      basemod <- bam(formula=as.formula(fn),
                     data=train_set,
                     na.action=na.exclude,
                     family=nb(),
                     gamma=2,
                     rho=0.8,
                     drop.unused.levels=FALSE,
                     select=TRUE,
                     method="fREML")
      
      mod0 <- bam(formula=as.formula(formulae[[j]]),
                  data=train_set,
                  na.action=na.exclude,
                  family=nb(),
                  gamma=2,
                  rho=0.8,
                  drop.unused.levels=FALSE,
                  select=TRUE,
                  method="fREML")
      
      
      train_set$preds <- predict.gam(mod0, 
                                     type="response")
      
      test_set$pred0 <- predict.gam(mod0, newdata=test_set,
                                    type='response')
      
      train_set %<>% drop_na(total) %>% drop_na(preds)
      
      test_set %<>% drop_na(total) %>% drop_na(pred0)
      
      myMae1 <- ModelMetrics::mae(test_set$total,
                                    test_set$pred0)
      aic1 <- AIC(mod0)
      
      bic1 <- BIC(mod0)
      
      cor_in <- cor(train_set$total,
                    train_set$preds,
                    method="spearman")
      
      cor_out <- cor(test_set$total,
                     test_set$pred0,
                     method="spearman")
      
      
      write.csv(test_set,
                file.path(myDir, "output_new",
                          paste0("test_set_equation_", j, "_block_",
                                 i, ".csv")),
                row.names=FALSE)
      
      write.csv(myMae1,
                file.path(myDir, "output_new",
                          paste0("mae_equation_", j, "_block_",
                                 i, ".csv")),
                row.names=FALSE)
      
      write.csv(aic1,
                file.path(myDir, "output_new",
                          paste0("aic_equation_", j, "_block_",
                                 i, ".csv")),
                row.names=FALSE)
      
      write.csv(bic1,
                file.path(myDir, "output_new",
                          paste0("bic_equation_", j, "_block_",
                                 i, ".csv")),
                row.names=FALSE)
      
      write.csv(cor_in,
                file.path(myDir, "output_new",
                          paste0("cor_in_equation_", j, "_block_",
                                 i, ".csv")),
                row.names=FALSE)
      
      
      write.csv(cor_out,
                file.path(myDir, "output_new",
                          paste0("cor_out_equation_", j, "_block_",
                                 i, ".csv")),
                row.names=FALSE)
      
      
    })
  }



