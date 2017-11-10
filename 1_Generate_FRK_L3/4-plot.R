## Produces daily map plots from the Fixed Rank Kriging results.

# Change this to "oco2v7" or "oco2v8"
data_version <- "oco2v7"

library(FRK)
library(dplyr)
library(ggplot2)

theme_set(theme_grey(base_size = 20))

my_colours <- c("#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7","#0f81f3",
                "#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e","#aefc2a","#e9fc0d",
                "#f6da0c","#f5a009","#f6780a","#f34a09","#f2210a","#f50008","#d90009","#a80109","#730005")

my_theme <- theme(panel.background = element_rect(fill = "white",colour = "white"), panel.grid = element_blank(), axis.ticks = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
		  plot.title = element_text(hjust = 0.5))

plotOneDay <- function(selecteddate,sixteendays) {
  print("Plotting one day")

  oneday <- sixteendays[as.Date(sixteendays$day,tz="UTC")==as.Date(selecteddate,tz="UTC"),]

	ggsave(
    (ggplot(oneday) +
       my_theme +
       geom_point(aes(lon,lat,colour=pmin(pmax(xco2,390),410))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_colour_gradientn(colours=my_colours, limits=c(390,410))  +
       labs(x="lon (deg)", y="lat (deg)", colour="XCO2\n(ppm)\n", title=paste(selecteddate,data_version,"DAILY DATA"))+
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = file.path(paste0(data_version,"plots"),paste0(selecteddate,"_data.png")), width=16, height=9, dpi=120)
}

plotSixteenDays <- function(selecteddate,sixteendays) {
  print("Plotting 16 days")

	ggsave(
    (ggplot(sixteendays) +
       my_theme +
       geom_point(aes(lon,lat,colour=pmin(pmax(xco2,390),410))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_colour_gradientn(colours=my_colours, limits=c(390,410))  +
       labs(x="lon (deg)", y="lat (deg)", colour="XCO2\n(ppm)\n", title=paste(selecteddate,data_version,"16-DAY MOVING WINDOW"))+
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = file.path(paste0(data_version,"plots"),paste0(selecteddate,"_16days.png")), width=16, height=9, dpi=120)
}

plotPredictions <- function(selecteddate, level3) {
  print("Plotting FRK Predictions")

  ggsave(
    (ggplot(level3) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(mu,390),410))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradientn(colours=my_colours, limits=c(390,410)) +
       labs(x="lon (deg)", y="lat (deg)", fill="pred\n(ppm)\n", title=paste(selecteddate,data_version," FIXED RANK KRIGING (FRK)")) +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = file.path(paste0(data_version,"plots"),paste0(selecteddate,"_prediction.png")), width=16, height=9, dpi=120)
}

plotUncertainty <- function(selecteddate, level3){
  print("Plotting FRK Uncertainty")

  ggsave(
    (ggplot(level3) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(sd,0.00),2.00))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradient(low="Green",high="Brown", limits=c(0.00,2.00)) +
       labs(x="lon (deg)", y="lat (deg)", fill="s.e.\n(ppm)\n", title=paste(selecteddate,data_version," FRK STANDARD ERROR")) +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = file.path(paste0(data_version,"plots"),paste0(selecteddate,"_uncertainty.png")),width=16,height=9,dpi=120)
}

plotAnomaly <- function(selecteddate, level3) {
  print("Plotting Anomaly")

  mu_mean <- mean(level3$mu)
  level3$anomaly <- level3$mu - mu_mean

  ggsave(
    (ggplot(level3) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(anomaly,-5),5))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradientn(colours=my_colours, limits=c(-5,5)) +
       labs(x="lon (deg)", y="lat (deg)", fill="anomaly\n(ppm)\n", title=paste0(selecteddate,data_version," Anomaly (pred - pred mean ",round(mu_mean,2),"ppm)")) +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = file.path(paste0(data_version,"plots"),paste0(selecteddate,"_anomaly.png")), width=16,height=9,dpi=120)
}

oco2lite <- read.csv(paste0(data_version,'lite.csv'))
oco2lite$day <- as.Date(oco2lite$day, tz="UTC")
if (!dir.exists(paste0(data_version,"plots"))) {
  dir.create(paste0(data_version,"plots"))
}

inputfiles <- list.files(path=paste0(data_version,"level3"), pattern="*.csv$", full.names=FALSE, recursive=FALSE)

for (i in 1:length(inputfiles)) {
  selecteddate <- as.Date(strsplit(inputfiles[i],"[.]")[[1]][1], tz="UTC")

  if ( file.exists(file.path(paste0(data_version,"plots"),paste0(selecteddate,"_anomaly.png"))) ) {
    # This date has already been plotted.
    next
  }
  file.create(file.path(paste0(data_version,"plots"),paste0(selecteddate,"_anomaly.png")))

  print(selecteddate)

  startdate <- as.Date(selecteddate,tz="UTC")-7
  enddate <- as.Date(selecteddate,tz="UTC")+8

  sixteendays <- oco2lite[oco2lite$day >= startdate & oco2lite$day <= enddate,]

  # Create a dummy data frame if there is no data
  if (is.null(sixteendays)) {
    sixteendays <- data.frame("date"=selecteddate,"lat"=0,"lon"=0,"xco2"=0,"std"=0)
  }

  plotOneDay(selecteddate,sixteendays)
  plotSixteenDays(selecteddate,sixteendays)

  if (!file.exists(file.path(paste0(data_version,"level3"),paste0(selecteddate,".csv"))) | file.size(file.path(paste0(data_version,"level3"),paste0(selecteddate,".csv"))) == 0) {
    # Input data does not exist for this date, create an empty data frame instead.
    level3 <- data.frame("date"=selecteddate,"lat"=0,"lon"=0,"mu"=0,"sd"=0)
  } else {
    level3 <- read.csv(file.path(paste0(data_version,"level3"),paste0(selecteddate,".csv")))
    level3$date <- as.Date(level3$date, tz="UTC")
  }

  plotPredictions(selecteddate, level3)
  plotUncertainty(selecteddate, level3)
  plotAnomaly(selecteddate, level3)
}
