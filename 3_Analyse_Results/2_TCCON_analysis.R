library(ggplot2)
library(dplyr)
library(lubridate)
library(xtable)
rm(list=ls())

source("TCCON_analysis_utils.R")
compare_from <- "FRKv7_no_target" # Can be FRKv7_no_target OR FRKv8_no_target
compare_to <- "TCCON" # Leave as TCCON
save_images <- 1

#####################
### Part 1: LOAD DATA
#####################
if(compare_from == "FRKv7_no_target") {
  dname <- "../data/v7r_frk/"
} else if(compare_from == "FRKv8_no_target") {
  dname <- "../data/v8r_frk/"
}

allcsv <- dir(dname,pattern = "*.csv")
oco2csv <- paste0(dname,allcsv[which(grepl("oco2",allcsv) & !grepl("lite",allcsv))])
#tcconcsv <- paste0("../data/tccon/",allcsv[which(grepl("tccon",allcsv))])
tcconcsv <- paste0("../data/tccon/", dir("../data/tccon",pattern = "*.csv"))

##################################
### Part 2: Summarise TCCON by day
##################################

## For each TCCON station find mean and total variability
comp_df <- NULL
for(i in 1:length(oco2csv)) {
  TCCON <- read.csv(tcconcsv[i])
  
  ## Fix names because column names changed between folders...
  names(TCCON) <- c(names(TCCON)[1:7],"oco2_latitude","oco2_longitude")
  
  TCCON$date <- date(TCCON$date)
  TCCON_summ <- group_by(TCCON,date) %>% summarise(tccon_xco2_av = mean(tccon_xco2),
                                                   tccon_error_av = sqrt(mean(tccon_error^2) + 
                                                                              safe.var(tccon_xco2)),
                                                   #tccon_error_av = pmax(mean(tccon_error),0.4),
                                                   tccon_latitude = tccon_latitude[1],
                                                   tccon_longitude = tccon_longitude[1],
                                                   filename = tcconcsv[i],
                                                   name = tools::file_path_sans_ext(basename(tcconcsv[i])))
  
  TCCON_summ$name <- sapply(strsplit(TCCON_summ$name,"0"),function(x) x[1])
  TCCON_summ$name <- as.factor(firstup(TCCON_summ$name))
  
  ## Now load OCO2 and join up with TCCON after changing names...
  OCO2 <- read.csv(oco2csv[i])
  names(OCO2) <- c("oco2_longitude","oco2_latitude","oco2_xco2","oco2_sd","date")
  OCO2$date <- date(OCO2$date)
  comp_df <- rbind(comp_df,inner_join(OCO2,TCCON_summ))
  if(any(is.na(comp_df$name))) stop()
}

##################################
### Part 3: Preprocess data frame
##################################

## Remove days with NAs in either TCCON or OCO2
#comp_df_sub <- filter(comp_df,!is.na(oco2_xco2) & !is.na(tccon_xco2_av)) 

## Fix time range
#comp_df_sub <- filter(comp_df_sub, !(month(date) >= 3 & year(date) >= 2017)) 

comp_df_sub <- filter(comp_df, !(month(date) >= 3 & year(date) >= 2017)) 

## Sort TCCON by latitude
TCCON_latitudes  <- unique(comp_df_sub[,c("name","tccon_latitude")])
TCCON_latitudes$order <- order(TCCON_latitudes$tccon_latitude)
comp_df_sub <- left_join(comp_df_sub,TCCON_latitudes)
comp_df_sub <- arrange(comp_df_sub,order)
comp_df_sub$name <- factor(comp_df_sub$name,levels=levels(TCCON_latitudes$name)[TCCON_latitudes$order],ordered=TRUE)

## Add season information
comp_df_sub <- mutate(comp_df_sub,
                      season = case_when(month(date)==1 ~ "1",
                                         month(date)==2 ~ "1",
                                         month(date)==3 ~ "2",
                                         month(date)==4 ~ "2",
                                         month(date)==5 ~ "2",
                                         month(date)==6 ~ "3",
                                         month(date)==7 ~ "3",
                                         month(date)==8 ~ "3",
                                         month(date)==9 ~ "4",
                                         month(date)==10 ~ "4",
                                         month(date)==11 ~ "4",
                                         month(date)==12 ~ "1"))
                                         

########################################
### Part 4: Compute diagnostic summaries
########################################

station_summary <- function(df) {
  data.frame(N = nrow(df),
             Bias = with(df,mean((oco2_xco2 - tccon_xco2_av))),
             MAE = with(df,mean(abs(oco2_xco2 - tccon_xco2_av))),
             RMSE = with(df,sqrt(mean((oco2_xco2 - tccon_xco2_av)^2))),
             R2 = with(df,cor(tccon_xco2_av, oco2_xco2)^2),
             Slope = with(df,lm(oco2_xco2~tccon_xco2_av + 0)) %>% coefficients(),
             Cov95 = with(df,cov95(tccon_xco2_av,oco2_xco2 - mean(oco2_xco2 - tccon_xco2_av),
                                            sqrt(oco2_sd^2 + tccon_error_av^2))))
}

format_rows <- function(df) {
  df %>% xtable(digits=c(0,0,0,2,2,2,2,3,2)) %>% 
    print.xtable(include.rownames=FALSE)
}

Rows1 <- group_by(comp_df_sub,name) %>% do(station_summary(.)) %>% 
         arrange(desc(name)) %>%
         as.data.frame() %>% 
         format_rows() %>%
         strsplit("\\n")
fileConn <- file(paste0(compare_from,".tex"))
writeLines(Rows1[[1]][9:33], fileConn)
close(fileConn)

if(grepl("v7",compare_from)) {
  rowtitle1 <- "Total v7r"
  rowtitle2 <- "Total v7r (w/o Pas.)"
} else {
  rowtitle1 <- "Total v8r"
  rowtitle2 <- "Total v8r (w/o Pas.)"
}


Rows2 <- cbind(data.frame(name = rowtitle1),
               station_summary(comp_df_sub)) %>%
         format_rows() %>%
         strsplit("\\n")

Rows3 <- filter(comp_df_sub,!(name == "Pasadena")) %>%
  station_summary()  %>% 
  cbind(data.frame(name = rowtitle2), .) %>% 
  format_rows() %>%
  strsplit("\\n")


fileConn <- file(paste0(compare_from,"_TOTAL.tex"))
writeLines(c(Rows2[[1]][9], Rows3[[1]][9]), fileConn)
close(fileConn)



# Plot scatter plot of differences B
gScatter <- ggplot(comp_df_sub %>% group_by(name,month(date),year(date)) %>% 
                     summarise(tccon_xco2_av = median(tccon_xco2_av),
                               oco2_xco2 = median(oco2_xco2))) + 
  geom_point(aes(tccon_xco2_av,oco2_xco2,colour=name,pch=name),size=2) +  
  theme_bw() + geom_abline(slope=1) + scale_shape_manual(values=seq(0,24)) + 
  xlab("TCCON retrievals averaged by station and month (ppm)") + 
  ylab("Level 3 XCO2 predictions averaged by station and month (ppm)") +
  theme(plot.title = element_text(size=20))


if(compare_from == "FRKv8_no_target") {
  ggsave(gScatter  +  ggtitle("(b)") + theme(plot.title = element_text(size=20)),
         filename=paste0("../img/TCCON_scatterv8.png"),width=8,height=6)
} else {
  ggsave(gScatter + theme(legend.position = "none") +  ggtitle("(a)") + theme(plot.title = element_text(size=20)),
         filename=paste0("../img/TCCON_scatterv7.png"),width=6,height=6)
}

## Do some exploratory plots
if(save_images) {
  
  ## All plots of predictions
  gpred <- ggplot(comp_df_sub) + geom_point(aes(date,tccon_xco2_av)) + 
    geom_point(aes(date,oco2_xco2),col="red") + facet_wrap(~name) + theme_bw()
  
  # All plots of differences
  gdiffs <- ggplot(comp_df_sub) + geom_point(aes(date,tccon_xco2_av - oco2_xco2))  + 
    facet_wrap(~name) + theme_bw()
  
  # Box plots of ALL differences
  gboxes <- ggplot(comp_df_sub) + geom_boxplot(aes(y = oco2_xco2 - tccon_xco2_av,x=name)) + 
            theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  
  # Box plots of differences by season
  gseasons <- ggplot(comp_df_sub) + 
    geom_boxplot(aes(y = oco2_xco2 - tccon_xco2_av,x=name,fill=season)) + 
    theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
    scale_fill_manual(values=c("blue","orange","yellow","purple"))
  
  
  # Plot scatter plot of differences A
  gScatter_ALL <- ggplot(comp_df_sub) + 
    geom_point(aes(tccon_xco2_av,oco2_xco2,colour=name,pch=name),size=2) +
    theme_bw() + geom_abline(slope=1) + scale_shape_manual(values=seq(0,22))
  lm(comp_df_sub$oco2_xco2 ~ comp_df_sub$tccon_xco2_av + 0) %>% summary()
  
  # Plot Time series, differences, and histograms, of Lamont and Wollongong
  options(digits=10)
  part_label = list(Lamont = c("(a)","(b)","(c)"),
                    Wollongong = c("(d)","(e)","(f)"))
  for(name_station in c("Lamont","Wollongong")) {
    X <- filter(comp_df_sub,name==name_station)
    Summary_results <- data.frame(x = c(2.5,2.5,2.5),
                                  y=c(0.5,0.45,0.4),
                                  txt=c(paste0("Mean = ",mean(X$oco2_xco2 - X$tccon_xco2_av) %>% round(2) %>% format(nsmall=2)),
                                        paste0("Median = ",median(X$oco2_xco2 - X$tccon_xco2_av) %>% round(2) %>% format(nsmall=2)),
                                        paste0("SD = ",sd(X$oco2_xco2 - X$tccon_xco2_av) %>% round(2) %>% format(nsmall=2))))
    
    g1 <- ggplot(X) + geom_point(aes(date,oco2_xco2)) + 
      geom_point(aes(date,tccon_xco2_av),col="red",alpha=0.36) + 
      geom_errorbar(aes(date,
                        ymin = tccon_xco2_av - 2*tccon_error_av,
                        ymax = tccon_xco2_av + 2*tccon_error_av),alpha=0.4,col="red") +
      geom_errorbar(aes(date,
                        ymin = oco2_xco2 - oco2_sd,
                        ymax = oco2_xco2 + oco2_sd),alpha=0.4,col="black") +
      theme_bw() + ggtitle(part_label[[name_station]][1]) + 
      theme(text = element_text(size=20)) +
      ylab("XCO2 (ppm)") + geom_text(data = data.frame(x=date("2015-01-01"),
                                                       y=408,
                                                       txt=paste0("r = ",cor(X$oco2_xco2,X$tccon_xco2_av) %>% round(2) %>% format(nsmall=2))),aes(x,y,label=txt),size=7) +
      coord_cartesian(ylim = c(392,409))
    
    g2 <- ggplot(X) + geom_point(aes(date,oco2_xco2 - tccon_xco2_av))   +
      geom_errorbar(aes(date,
                        ymin = oco2_xco2 - tccon_xco2_av  - 2*sqrt(oco2_sd^2 + tccon_error_av^2),
                        ymax = oco2_xco2 - tccon_xco2_av  + 2*sqrt(oco2_sd^2 + tccon_error_av^2)),
                    alpha=0.4,col="black") +
      ylab(expression(paste(Delta," XCO2 (ppm)"))) +
      theme_bw() + ggtitle(part_label[[name_station]][2]) + 
      theme(text = element_text(size=20)) +
      coord_cartesian(ylim = c(-4,4.5))
    
    g3 <- ggplot(X) + geom_histogram(aes(oco2_xco2 - tccon_xco2_av,y=..density..),
                                     col="black",bins=20) + 
      theme_bw() + xlab(expression(paste(Delta," XCO2 (ppm)"))) + 
      geom_text(data = Summary_results,aes(x,y,label=txt),hjust="right") +
      theme_bw() + ggtitle(part_label[[name_station]][3]) + 
      theme(text = element_text(size=20)) + ylim(c(0,0.7)) + xlim(c(-3,3))
    
    ggsave(g1,filename=paste0("../img/",name_station,"_ts.png"),width=8,height=4)
    ggsave(g2,filename=paste0("../img/",name_station,"_diff.png"),width=8,height=4)
    ggsave(g3,filename=paste0("../img/",name_station,"_hist.png"),width=8,height=4)
  }
  
}

## UNCERTAINTY DIAGNOSTICS taking into account sample size (on hold)
if(0) {
  # QQplot
  my_qqplot <- function(q,x,y,y_sd) {
    cum <- q*0
    for(i in seq_along(q)){
      cum[i] <- mean(x > y - abs(qnorm(p = q[i],mean = 0,sd = y_sd)) &
                       x < y + abs(qnorm(p = q[i],mean = 0,sd = y_sd)))
    }
    cum
  }
  
  quants <- seq(0.5,1,length=100)
  cum <- quants*0
  cum_df <- NULL
  reg_df <- data.frame(date = seq(date("2014-10-01"),date("2017-02-28"),by="day"))
  
  for(name_station in levels(comp_df_sub$name)) {
    x11()
    sub_df <- subset(comp_df_sub,name == name_station)
    cum <- with(sub_df,
                my_qqplot(quants,tccon_xco2_av,oco2_xco2,
                          sqrt(oco2_sd^2 + tccon_error_av^2 )))
    cum2 <- with(sub_df,
                 my_qqplot(quants,tccon_xco2_av,oco2_xco2 - mean(oco2_xco2 - tccon_xco2_av),
                           sqrt(oco2_sd^2 + tccon_error_av^2)))
    plot((quants - 0.5)*2,cum,asp=1,xlim=c(0,1),ylim=c(0,1),col="green",main=name_station)
    lines((quants - 0.5)*2,cum2,asp=1,xlim=c(0,1),ylim=c(0,1),col="blue",main=name_station)
    
    ## Effective sample size
    if(0){
      #if(nrow(sub_df) > 30) {
      sub_df$diff <- sub_df$oco2_xco2 - sub_df$tccon_xco2_av
      reg_df2 <- left_join(reg_df,sub_df)
      tau <- 1 + 2*sum(acf(reg_df2$diff,na.action = na.pass)$acf[2:5],na.rm=TRUE)
      neff <-  nrow(sub_df) / tau
      print(tau)
    } else {
      neff <- nrow(sub_df)
    }
    
    cum_mat <- matrix(0,length(quants),500)
    warning("NOT Using only neff samples")
    for(i in 1:ncol(cum_mat)) {
      sim <- rnorm(mean=sub_df$oco2_xco2,
                   sd = sqrt(sub_df$oco2_sd^2 + sub_df$tccon_error_av^2),
                   n = nrow(sub_df))
      idx <- sample(1:nrow(sub_df),neff)
      cum_mat[,i] <-  my_qqplot(quants,sim[idx],sub_df$oco2_xco2[idx],
                                sqrt(sub_df$oco2_sd^2 + sub_df$tccon_error_av^2)[idx])
    }
    envup <- apply(cum_mat,1,function(x) quantile(x,0.995))
    envlo <- apply(cum_mat,1,function(x) quantile(x,0.005))
    lines((quants - 0.5)*2,envup,ylim=c(0,1),asp=1,xlim=c(0,1))
    lines((quants - 0.5)*2,envlo,asp=1,xlim=c(0,1),ylim=c(0,1))
    lines(c(0,1),c(0,1),col="red")
    
    cat(cum[90],cum2[90],cum2[90] > min(cum_mat[90,]) & cum2[90] < max(cum_mat[90,]),"\n")
    
    cum_df <- rbind(cum_df,
                    data.frame(station = name_station, 
                               quants = quants, 
                               cum = cum, cum2 = cum2,
                               envlo = envlo,envup = envup))
    
  }
  
  ggplot(cum_df %>% filter(!station %in% c("Pasadena","Sodankyla"))) + geom_line(aes(quants,cum),col="red") + 
    geom_line(aes(quants,cum2),col="blue",linetype=2) + 
    geom_line(aes(quants,envlo)) + geom_line(aes(quants,envup)) +
    facet_wrap(~station) + theme_bw()
  
  
}