## Filter the OCO-2 Lite, OCO-2 FRK, and TCCON data to produce CSV files at each TCCON site for analysis.

library(ncdf4)
library(ggplot2)
library(dplyr)
library(readr)

# Calculate the local solar time from UTC time and longitude.
solar_time_function <- function(inputDate, longitude) {
  inputTime <- as.numeric(format(inputDate, format = "%H")) + as.numeric(format(inputDate, format = "%M"))/60
  offset <- longitude / 180 * 12
  (inputTime + offset) %% 24
}

# Read in the OCO-2 V8r Lite data. This is only used here to find out what time of day
# the satellite passes over each TCCON site for filtering TCCON readings.
oco2lite <- read.csv('oco2v8lite.csv')
oco2lite$day <- as.POSIXct(oco2lite$day, tz="UTC")
oco2lite$solar_time <- solar_time_function(oco2lite$day,oco2lite$lon)

# Read in the OCO-2 FRK V7r 'Level 3' data
frkfiles7 <- list.files(path="oco2v7level3", pattern="*.csv$", full.names=TRUE, recursive=FALSE)
oco2frk7 = lapply(frkfiles7, read_csv) %>% bind_rows()
oco2frk7$date<-as.POSIXct(oco2frk7$date,tz="UTC")

# Read in the OCO-2 FRK V8r 'Level 3' data
frkfiles8 <- list.files(path="oco2v8level3", pattern="*.csv$", full.names=TRUE, recursive=FALSE)
oco2frk8 = lapply(frkfiles8, read_csv) %>% bind_rows()
oco2frk8$date<-as.POSIXct(oco2frk8$date,tz="UTC")

tcconfiles <- list.files(path="tccon", pattern="*.nc$", full.names=TRUE, recursive=FALSE)

for (i in 1:length(tcconfiles)) {
  tccon_file <- tcconfiles[i]
  nc <- nc_open(tccon_file)
  tccon_attributes <- ncatt_get(nc,0)
  print(tccon_attributes$Location)
  tccon_time <- as.POSIXct(ncvar_get(nc,"time")*24*3600, origin = "1970-01-01", tz="UTC")
  long_deg <- ncvar_get(nc,"long_deg")
  lat_deg <- ncvar_get(nc,"lat_deg")
  xco2_ppm <- ncvar_get(nc,"xco2_ppm")
  xco2_ppm_error <- ncvar_get(nc,"xco2_ppm_error")

  tccon_solar_time <- solar_time_function(tccon_time,long_deg[1])

  tccon_data <- data.frame(xco2_ppm,xco2_ppm_error,tccon_time, tccon_solar_time)
  sitename <- tccon_attributes$longName

  # Select OCO-2 Lite data at the latitude of the TCCON site to calculate time of day when the satellite passes over it.
  oco2lite_subset <- oco2lite[oco2lite$lat > lat_deg[1]-.5 & oco2lite$lat < lat_deg[1]+.5 ,]
  if (nrow(oco2lite_subset) == 0) {
    print("There is no OCO-2 Lite data at this latitude.")
    next
  }

  oco2lite_mean_solar_time <- summary(oco2lite_subset$solar_time)[[4]]

  # Select TCCON data within +/- 30 minutes of the mean OCO-2 Lite local solar time.
  tccon_subset <- tccon_data[tccon_data$tccon_solar_time>(oco2lite_mean_solar_time-0.5) &
                                 tccon_data$tccon_solar_time<(oco2lite_mean_solar_time+0.5) &
                                 tccon_data$xco2_ppm_error<1,]

  # Find the nearest 1x1 grid centre for comparison with OCO-2 FRK
  tccon_rounded_lat <- floor(lat_deg[1])+.5
  tccon_rounded_lon <- floor(long_deg[1])+.5

  names(tccon_subset) <- c("tccon_xco2","tccon_error","tccon_datetime","tccon_solar_time")
  tccon_subset$tccon_latitude <- lat_deg[1]
  tccon_subset$tccon_longitude <- long_deg[1]
  tccon_subset$date <- as.Date(tccon_subset$tccon_datetime,tz="UTC")
  tccon_subset$latitude <- tccon_rounded_lat
  tccon_subset$longitude <- tccon_rounded_lon
  write.csv(tccon_subset,file=paste0(tccon_attributes$longName,"_tccon.csv"), row.names=FALSE)

  oco2frk7_subset <- oco2frk7[oco2frk7$lat==tccon_rounded_lat &
                         oco2frk7$lon==tccon_rounded_lon ,]

  names(oco2frk7_subset) <- c('longitude','latitude','oco2_frk_xco2','oco2_frk_sd','date')
  write.csv(oco2frk7_subset,file=paste0(tccon_attributes$longName,"_oco2_v7r_frk.csv"), row.names=FALSE)

  oco2frk8_subset <- oco2frk8[oco2frk8$lat==tccon_rounded_lat &
                         oco2frk8$lon==tccon_rounded_lon ,]

  names(oco2frk8_subset) <- c('longitude','latitude','oco2_frk_xco2','oco2_frk_sd','date')
  write.csv(oco2frk8_subset,file=paste0(tccon_attributes$longName,"_oco2_v8r_frk.csv"), row.names=FALSE)

  if (nrow(tccon_subset) != 0 & nrow(oco2frk7_subset) != 0 & nrow(oco2frk8_subset) != 0) {

  # Plot data points near each TCCON site with +/- 2 s.e. error bars.
  ggsave(
    ggplot() +
      theme_bw(base_size = 20) +
      theme(legend.position="bottom") +
      geom_pointrange(data = tccon_subset,
                      aes(x = tccon_datetime, y=tccon_xco2, ymin = tccon_xco2 - 2*tccon_error, ymax = tccon_xco2 + 2*tccon_error, colour="TCCON"),
                      size=.2, alpha=.3) +
      geom_pointrange(data = oco2frk7_subset,
                      aes(x = date, y = oco2_frk_xco2, ymin = oco2_frk_xco2 - 2*oco2_frk_sd, ymax = oco2_frk_xco2 + 2*oco2_frk_sd, colour="OCO2 V7r FRK"), size=.2, alpha=.3) +
      geom_pointrange(data = oco2frk8_subset,
                      aes(x = date, y = oco2_frk_xco2, ymin = oco2_frk_xco2 - 2*oco2_frk_sd, ymax = oco2_frk_xco2 + 2*oco2_frk_sd, colour="OCO2 V8r FRK"), size=.2, alpha=.3) +
      scale_colour_manual(values=c("OCO2 V7r FRK"="Blue","OCO2 V8r FRK"="Green","TCCON"="Red"), name="")+
      labs(x="Date (Year-Month)", y="XCO2 (ppm)", title=tccon_attributes$Location)+
      coord_cartesian(ylim=c(390,410),xlim=c(as.POSIXct('2015-01-01'), as.POSIXct('2017-01-01')))
    ,filename = paste0(tccon_attributes$longName,"_oco2_tccon.png"), width=12, height=8)
  }
}
