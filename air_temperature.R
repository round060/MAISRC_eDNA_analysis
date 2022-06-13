require(raster)
require(rgdal)
require(ncdf4)

# I'd install the new geoknife - it will give you a progress bar for the job when you are using `wait` = TRUE
#devtools::install_github('usgs-r/geoknife')

library(geoknife)
library(reshape2)
library(tidyverse)
library(readr)
library(lubridate)


sites <- read.csv("./data/2021_lakes_lat_long.csv", header = T)

sites = select(sites, -c(Lake, site))
colnames(sites) = c("id", "lat", "lon")

for (i in 1:length(sites$id)) {
  df <- data.frame(c(sites$lon[i], sites$lat[i]))
  names(df) <- sites$id[i]
  if (i == 1)
    geom <- df
  else 
    geom <- cbind(geom, df)
}


stencil <- simplegeom(geom)
fabric <- webdata(url = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/NARR/Dailies/monolevel/air.2m.2021.nc', variables = 'air')
job <- geoknife(geom, fabric, wait = T)
airtemp_data <- result(job, with.units = TRUE)
airtemp_data <- airtemp_data %>% mutate_if(is.numeric, funs(. - 273.15)) 
#convert from K to C should give depreciation error
airlong = melt(select(airtemp_data, -variable, -statistic, -units), id.vars = "DateTime", variable.name = "DOW")

all.temps = rename(airlong, mean_daily_temp_c = value)

monthly_temps <- all.temps %>%
  group_by(DOW, month = lubridate::floor_date(DateTime, "month")) %>%
  summarize(summary_variable = mean(mean_daily_temp_c)) 


#write file
write.csv(all.temps,"./data/2021_temps.csv", row.names = F)
