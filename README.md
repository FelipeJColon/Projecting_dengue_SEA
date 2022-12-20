# Projecting the future incidence and burden of dengue in Southeast Asia

Repository for initial submission


## Software dependencies

For the development of the code in this repository we made use of the following `R` packages (in a macOS environment):

- `data.table` version 1.14.0
- `tidyverse` version 1.3.1
- `lubridate` version 1.7.10
- `ModelMetrics` version 1.2.2.2
- `rgdal` version 1.5-23
- `magrittr` version 2.0.1
- `mgcv` version 1.8-35
- `readr` version 1.4.0
- `spdep` version 1.1-11
- `reshape2` version 1.4.4

There is no need for non-standard hardware. 

## Installation

To run the code presented here simply run the `Load_pkgs.R` to install all required packages onto your computer. 
Then, source each file one at the time. The typicall installation time of all packages is about 7 minutes, but
times may vary depending on the specification of your computer and the OS.

## Demo

Save data in any directory in your desktop. Then make sure **shapefiles** are in a folder called `shapefiles` and 
all `.csv` files are stored in a `data` folder.

### Expected output

After running all files you should get the following output

- `00_blocked_cross_validation.R`: Generates  a suite of `.csv` files named `metric_equation_x_block_y.csv` 
where 'metric' indicates the name of the metric used (i.e., MAE, BIC, AIC, Spearman correlation),  `x` 
indicates the equation being investigated, and `y` the cross-validation block. Estimated running time: 3 days
using a High Performance Computer if you want to run all possible combinations of covariates.

- `01_Bam_model.R`: an `RDS` file named `bam_model_sea.rds` which includes the model specification. Estimated running 
time: 35 minutes.

- `02_Predictions_bam_model.R`: a suite of `.csv` files named `annual_estimates_model_rcp_sea.csv` where `model` indicates
the name of the GCM model used and `rcp` indicates the SSP scenario. Estimated running time: 2h.

- `03_sensitivity_analysis.R`: a suite of `.csv` files named `predictions_variable_model_rcp_sea.csv` where `variable` is
the name of the group of predictors kept constant, `model` indicates the name of the GCM model used and `rcp` indicates 
the SSP scenario. Estimated running time: 1h.


## Instructions for use

To run the code on our data set make sure you replace the lines corresponding to `myDir` and `output` on the preamble of
each of the files with the directory where you (i) have saved the data files and (ii) where you would like to save 
output files.

Make sure you pay attention to any **subfoldder** in the code which will be found after the `file.path()` command as in 
- `den_data <- fread(file.path(output, "data",
                            "dengue_analysis_data_SEA_lshtm.csv"))`
                            
where you will find that the file `dengue_analysis_data_SEA_lshtm.csv` is in subfolder `data` included in the path to `output`.


