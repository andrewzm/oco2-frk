## Produce plots showing the difference between V7r FRK and V8r FRK on 13th May 2016.

library(sp)
library(spacetime)
library(ggplot2)
library(dplyr)
library(FRK)

oco2frk7 <- read.csv("oco2v7level3/2016-05-13.csv")
names(oco2frk7) <- c("lon","lat","mu7","sd7","date")
oco2frk8 <- read.csv("oco2v8level3/2016-05-13.csv")
names(oco2frk8) <- c("lon","lat","mu8","sd8","date")
combined <- inner_join(oco2frk7,oco2frk8)

theme_set(theme_grey(base_size = 20))

my_colours <- c("#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7","#0f81f3",
                "#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e","#aefc2a","#e9fc0d",
                "#f6da0c","#f5a009","#f6780a","#f34a09","#f2210a","#f50008","#d90009","#a80109","#730005")

my_theme <- theme(panel.background = element_rect(fill = "white",colour = "white"), panel.grid = element_blank(), axis.ticks = element_blank(),
                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(),
		  plot.title = element_text(hjust = 0.5))

print("Plotting V8r FRK Predictions")
ggsave(
    (ggplot(oco2frk8) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(mu8,390),410))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradientn(colours=my_colours, limits=c(390,410)) +
       labs(x="lon (deg)", y="lat (deg)", fill="pred\n(ppm)\n", title="2016-05-13 OCO-2 V8r Lite Fixed Rank Kriging (FRK)") +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = "2016-05-13_v8_prediction.png", width=16, height=9, dpi=120)

print("Plotting V8r FRK Uncertainty")
ggsave(
    (ggplot(oco2frk8) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(sd8,0.00),2.00))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradient(low="Green",high="Brown", limits=c(0.00,2.00)) +
       labs(x="lon (deg)", y="lat (deg)", fill="s.e.\n(ppm)\n", title="2016-05-13 OCO-2 V8r Lite FRK Standard Error") +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = "2016-05-13_v8_uncertainty.png",width=16,height=9,dpi=120)

print("Plotting V7r FRK Predictions")
ggsave(
    (ggplot(oco2frk7) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(mu7,390),410))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradientn(colours=my_colours, limits=c(390,410)) +
       labs(x="lon (deg)", y="lat (deg)", fill="pred\n(ppm)\n", title="2016-05-13 OCO-2 V7r Lite Fixed Rank Kriging (FRK)") +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = "2016-05-13_v7_prediction.png", width=16, height=9, dpi=120)

print("Plotting V7r FRK Uncertainty")
ggsave(
    (ggplot(oco2frk7) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(sd7,0.00),2.00))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradient(low="Green",high="Brown", limits=c(0.00,2.00)) +
       labs(x="lon (deg)", y="lat (deg)", fill="s.e.\n(ppm)\n", title="2016-05-13 OCO-2 V7r Lite FRK Standard Error") +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = "2016-05-13_v7_uncertainty.png",width=16,height=9,dpi=120)

combined$diffmu <- combined$mu8-combined$mu7
combined$diffsd <- combined$sd8-combined$sd7

print("Plotting FRK Prediction difference (V8r - V7r)")
ggsave(
    (ggplot(combined) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(diffmu,-5),5))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradientn(colours=my_colours, limits=c(-3,3)) +
       labs(x="lon (deg)", y="lat (deg)", fill="diff\n(ppm)\n", title="2016-05-13 FRK Difference (V8r - V7r)") +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = "2016-05-13_diff_prediction.png", width=16, height=9, dpi=120)

print("Plotting FRK Uncertainty difference (V8r - V7r)")
ggsave(
    (ggplot(combined) +
       my_theme +
       geom_tile(aes(lon,lat,fill=pmin(pmax(diffsd,-1.00),1.00))) +
       lims(x = c(-180, 180), y = c(-90, 90)) +
       scale_fill_gradientn(colours=my_colours, limits=c(-.5,.5)) +
       labs(x="lon (deg)", y="lat (deg)", fill="diff\n(ppm)\n", title="2016-05-13 FRK Standard Error Difference (V8r - V7r)") +
       coord_map("mollweide")) %>%
      draw_world(inc_border=TRUE),
    filename = "2016-05-13_diff_uncertainty.png",width=16,height=9,dpi=120)
