library(ncdf4)
library(ggplot2)
library(dplyr)
library(readr)

if (!dir.exists("tccon")) {
  dir.create("tccon")
}

# Calculate the local solar time from UTC time and longitude.
solar_time_function <- function(inputDate, longitude) {
  inputTime <- as.numeric(format(inputDate, format = "%H")) + as.numeric(format(inputDate, format = "%M"))/60
  offset <- longitude / 180 * 12
  (inputTime + offset) %% 24
}

# Read in the OCO-2 Lite data
oco2lite <- read.csv('oco2lite.csv')
oco2lite$day <- as.POSIXct(oco2lite$day, tz="UTC")
oco2lite$solar_time <- solar_time_function(oco2lite$day,oco2lite$lon)

# Read in the OCO-2 FRK 'Level 3' data
frkfiles <- list.files(path="level3", pattern="*.csv$", full.names=TRUE, recursive=FALSE)
oco2frk = lapply(frkfiles, read_csv) %>% bind_rows()
oco2frk$date<-as.POSIXct(oco2frk$date,tz="UTC")

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

  # Find the nearest 1x1 grid centre for comparison with OCO-2 FRK
  tccon_rounded_lat <- floor(lat_deg[1])+.5
  tccon_rounded_lon <- floor(long_deg[1])+.5
  
  oco2frk_subset <- oco2frk[oco2frk$lat==tccon_rounded_lat &
                         oco2frk$lon==tccon_rounded_lon ,]
  if (nrow(oco2frk_subset) == 0) {
    print("There is no OCO-2 FRK data at this TCCON site.")
    next
  }
  
  # Select OCO-2 Lite data within the FRK 1x1 grid cell.
  target_range <- 0.5
  oco2lite_subset <- oco2lite[oco2lite$lat > tccon_rounded_lat-target_range & oco2lite$lat < tccon_rounded_lat+target_range &
                                oco2lite$lon > tccon_rounded_lon-target_range & oco2lite$lon < tccon_rounded_lon+target_range,]
  if (nrow(oco2lite_subset) == 0) {
    print("There is no OCO-2 Lite data at this TCCON site.")
    next
  }
  
  oco2lite_mean_solar_time <- summary(oco2lite_subset$solar_time)[[4]]
  
  # Select TCCON data within +/- 30 minutes of the mean OCO-2 Lite local solar time.
  tccon_filtered <- tccon_data[tccon_data$tccon_solar_time>(oco2lite_mean_solar_time-0.5) &
                                 tccon_data$tccon_solar_time<(oco2lite_mean_solar_time+0.5) &
                                 tccon_data$xco2_ppm_error<1,]
  
  if (nrow(tccon_filtered) == 0) {
    print("There is no TCCON data within the OCO-2 time window at this TCCON site.")
    next
  }

  oco2_export <- oco2frk_subset
  names(oco2_export) <- c('longitude','latitude','oco2_frk_xco2','oco2_frk_sd','date')
  write.csv(oco2_export,file=file.path("tccon",paste0(tccon_attributes$longName,"_oco2_frk.csv")), row.names=FALSE)
  
  tccon_export <- tccon_filtered
  names(tccon_export) <- c("tccon_xco2","tccon_error","tccon_datetime","tccon_solar_time")
  tccon_export$tccon_latitude <- lat_deg[1]
  tccon_export$tccon_longitude <- long_deg[1]
  tccon_export$date <- as.Date(tccon_export$tccon_datetime,tz="UTC")
  tccon_export$latitude <- tccon_rounded_lat
  tccon_export$longitude <- tccon_rounded_lon
  write.csv(tccon_export,file=file.path("tccon",paste0(tccon_attributes$longName,"_tccon.csv")), row.names=FALSE)

  oco2lite_subset_export <- oco2lite_subset
  names(oco2lite_subset_export) <- c("date", "longitude", "latitude", "oco2_lite_xco2", "oco2_lite_xco2_uncertainty", "oco2_lite_solar_time")
  write.csv(oco2lite_subset_export,file=file.path("tccon",paste0(tccon_attributes$longName,"_oco2_lite.csv")), row.names=FALSE)

  # Plot data points near each TCCON site with +/- 2 s.e. error bars.
  ggsave(
    ggplot() +
      theme_bw(base_size = 20) +
      theme(legend.position="bottom") +
      geom_pointrange(data = oco2lite_subset,
                      aes(x = day, y = xco2, ymin = xco2 - 2*std, ymax = xco2 + 2*std, colour="OCO2 L2 Lite"), size=.2, alpha=.3) +
      geom_pointrange(data = tccon_filtered,
                      aes(x = tccon_time, y=xco2_ppm, ymin = xco2_ppm - 2*xco2_ppm_error, ymax = xco2_ppm + 2*xco2_ppm_error, colour="TCCON"),
                      size=.2, alpha=.3) +
      geom_pointrange(data = oco2frk_subset,
                      aes(x = date, y = mu, ymin = mu - 2*sd, ymax = mu + 2*sd, colour="OCO2 FRK"), size=.2, alpha=.3) +
      scale_colour_manual(values=c("OCO2 FRK"="Blue","OCO2 L2 Lite"="Green","TCCON"="Red"), name="")+
      labs(x="Date (Year-Month)", y="XCO2 (ppm)", title=tccon_attributes$Location)+
      coord_cartesian(ylim=c(390,410),xlim=c(as.POSIXct('2015-01-01'), as.POSIXct('2017-01-01')))
    ,filename = file.path("tccon",paste0(tccon_attributes$longName,"_oco2_tccon.png")), width=12, height=8)
}
