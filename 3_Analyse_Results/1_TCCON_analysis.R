#######################################################################
### Compare the FRK Version 7r and Version 8r products to TCCON
#######################################################################

## Set working directory to this directory. Then load packages
library("ggplot2")
library("dplyr")
library("lubridate")
library("xtable")

## Clear workspace
rm(list=ls())

## Load utility functions
source("TCCON_analysis_utils.R")

compare_to <- "TCCON"             # Leave as TCCON
save_images <- 0                  # Set to 0 to not save images

#####################
### Part 1: LOAD DATA
#####################

dname_V7 <- "../data/v7r_frk/"
dname_V8 <- "../data/v8r_frk/"


## Get all CSVs in the product directory.
## If there are the "coincident" lite retrievals don't load those.
if(!identical(dir(dname_V7,pattern = "*.csv"),dir(dname_V8,pattern = "*.csv")))
    stop("CSV files in V7 and V8 folders do not have identical names")
allcsv <- dir(dname_V7,pattern = "*.csv")
oco2csv_V7 <- paste0(dname_V7,allcsv[which(grepl("oco2",allcsv) & !grepl("lite",allcsv))])
oco2csv_V8 <- paste0(dname_V8,allcsv[which(grepl("oco2",allcsv) & !grepl("lite",allcsv))])

## Get all CSVs in the TCCON directory
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

  ## Cast "date" to Date object
  TCCON$date <- date(TCCON$date)

  ## summarise TCCON. When computing error also consider the variability within the 60min. time window
  TCCON_summ <- group_by(TCCON,date) %>% summarise(tccon_xco2_av = mean(tccon_xco2),
                                                   tccon_error_av = sqrt(mean(tccon_error^2) +
                                                                              safe.var(tccon_xco2)),
                                                   tccon_latitude = tccon_latitude[1],
                                                   tccon_longitude = tccon_longitude[1],
                                                   filename = tcconcsv[i],
                                                   name = tools::file_path_sans_ext(basename(tcconcsv[i])))

  ## Add name of station to data frame
  TCCON_summ$name <- sapply(strsplit(TCCON_summ$name,"0"),function(x) x[1])
  TCCON_summ$name <- as.factor(firstup(TCCON_summ$name))

  ## Now load OCO2 and join up with TCCON after changing names...
  OCO2_V7 <- read.csv(oco2csv_V7[i])
  OCO2_V8 <- read.csv(oco2csv_V8[i])
  names(OCO2_V7) <- c("oco2_longitudeV7","oco2_latitudeV7",
                      "oco2_xco2V7","oco2_sdV7","dateV7")
  names(OCO2_V8) <- c("oco2_longitudeV8","oco2_latitudeV8",
                      "oco2_xco2V8","oco2_sdV8","dateV8")
  OCO2_V7$date <- date(OCO2_V7$date)
  OCO2_V8$date <- date(OCO2_V8$date)
  comp_df <- comp_df %>%
                rbind(inner_join(OCO2_V7, TCCON_summ) %>%
                      inner_join(OCO2_V8))
  if(any(is.na(comp_df$name))) stop("Something is wrong. Please check code in Part 2")
}

##################################
### Part 3: Preprocess data frame
##################################

## Fix time range
comp_df_sub <- filter(comp_df, !(month(date) >= 3 & year(date) >= 2017))

## Sort TCCON by latitude
TCCON_latitudes  <- unique(comp_df_sub[,c("name","tccon_latitude")])
TCCON_latitudes$order <- order(TCCON_latitudes$tccon_latitude)
comp_df_sub <- left_join(comp_df_sub,TCCON_latitudes)
comp_df_sub <- arrange(comp_df_sub,order)
comp_df_sub$name <- factor(comp_df_sub$name,
                           levels = levels(TCCON_latitudes$name)[TCCON_latitudes$order],ordered=TRUE)

## Add season information (not used in the analysis)
comp_df_sub <- mutate(comp_df_sub,
                      season = case_when(month(date) == 1 ~ "1",
                                         month(date) == 2 ~ "1",
                                         month(date) == 3 ~ "2",
                                         month(date) == 4 ~ "2",
                                         month(date) == 5 ~ "2",
                                         month(date) == 6 ~ "3",
                                         month(date) == 7 ~ "3",
                                         month(date) == 8 ~ "3",
                                         month(date) == 9 ~ "4",
                                         month(date) == 10 ~ "4",
                                         month(date) == 11 ~ "4",
                                         month(date) == 12 ~ "1"))

########################################
### Part 4: Compute diagnostic summaries
########################################

## Function that takes the data frame and computes N, Bias, MAE etc.
station_summaryV7 <- function(df) {
  data.frame(N = nrow(df),
             Bias = with(df,mean((oco2_xco2V7 - tccon_xco2_av))),
             MAE = with(df,mean(abs(oco2_xco2V7 - tccon_xco2_av))),
             RMSE = with(df,sqrt(mean((oco2_xco2V7 - tccon_xco2_av)^2))),
             R2 = with(df,cor(tccon_xco2_av, oco2_xco2V7)^2),
             Slope = with(df,lm(oco2_xco2V7~tccon_xco2_av + 0)) %>% coefficients(),
             Cov95 = with(df,cov95(tccon_xco2_av,oco2_xco2V7 - mean(oco2_xco2V7 - tccon_xco2_av),
                                            sqrt(oco2_sdV7^2 + tccon_error_av^2))))
}

station_summaryV8 <- function(df) {
    data.frame(N = nrow(df),
               Bias = with(df,mean((oco2_xco2V8 - tccon_xco2_av))),
               MAE = with(df,mean(abs(oco2_xco2V8 - tccon_xco2_av))),
               RMSE = with(df,sqrt(mean((oco2_xco2V8 - tccon_xco2_av)^2))),
               R2 = with(df,cor(tccon_xco2_av, oco2_xco2V8)^2),
               Slope = with(df,lm(oco2_xco2V8~tccon_xco2_av + 0)) %>% coefficients(),
               Cov95 = with(df,cov95(tccon_xco2_av,oco2_xco2V8 - mean(oco2_xco2V8 - tccon_xco2_av),
                                     sqrt(oco2_sdV8^2 + tccon_error_av^2))))
}


## Formats rows for printing tables
format_rows <- function(df) {
  df %>% xtable(digits=c(0,0,0,2,2,2,2,3,2)) %>%
    print.xtable(include.rownames=FALSE)
}

## Generate the big table (diagnostics for each station)
Rows1V7df <- group_by(comp_df_sub,name) %>% do(station_summaryV7(.)) %>%
         arrange(desc(name)) %>%
         as.data.frame()

Rows1V7 <- Rows1V7df %>%
         format_rows() %>%
         strsplit("\\n")

Rows1V8df <- group_by(comp_df_sub,name) %>% do(station_summaryV8(.)) %>%
    arrange(desc(name)) %>%
    as.data.frame()

Rows1V8 <- Rows1V8df %>%
    format_rows() %>%
    strsplit("\\n")

## Write to file
fileConn <- file("FRKv7_no_target.tex")
writeLines(Rows1V7[[1]][9:33], fileConn)
close(fileConn)


fileConn <- file("FRKv8_no_target.tex")
writeLines(Rows1V8[[1]][9:33], fileConn)
close(fileConn)

## Generate the summary table (diagnostics across all stations) and save to file
rowtitle1V7 <- "Total v7r"
rowtitle2V7 <- "Total v7r (w/o Pas.)"
rowtitle1V8 <- "Total v8r"
rowtitle2V8 <- "Total v8r (w/o Pas.)"


## With Pasadena
Rows2V7 <- cbind(data.frame(name = rowtitle1V7),
               station_summaryV7(comp_df_sub)) %>%
         format_rows() %>%
         strsplit("\\n")

Rows2V8 <- cbind(data.frame(name = rowtitle1V8),
                 station_summaryV8(comp_df_sub)) %>%
    format_rows() %>%
    strsplit("\\n")

## Without Pasadena
Rows3V7 <- filter(comp_df_sub,!(name == "Pasadena")) %>%
  station_summaryV7()  %>%
  cbind(data.frame(name = rowtitle2V7), .) %>%
  format_rows() %>%
  strsplit("\\n")

Rows3V8 <- filter(comp_df_sub,!(name == "Pasadena")) %>%
    station_summaryV8()  %>%
    cbind(data.frame(name = rowtitle2V8), .) %>%
    format_rows() %>%
    strsplit("\\n")

## Write to file
fileConn <- file("FRKv7_no_target_TOTAL.tex")
writeLines(c(Rows2V7[[1]][9], Rows3V7[[1]][9]), fileConn)
close(fileConn)

fileConn <- file("FRKv8_no_target_TOTAL.tex")
writeLines(c(Rows2V8[[1]][9], Rows3V8[[1]][9]), fileConn)
close(fileConn)


# Plot scatter plot of differences B
gScatter <- function(x) {
  comp_df_sub$dep <- comp_df_sub[,x]
  ggplot(comp_df_sub %>% group_by(name,month(date),year(date)) %>%
                     summarise(tccon_xco2_av = median(tccon_xco2_av),
                               oco2_xco2 = median(dep))) +
  geom_point(aes(tccon_xco2_av,oco2_xco2,colour=name,pch=name),size=2) +
  theme_bw() + geom_abline(slope=1) + scale_shape_manual(values=seq(0,24)) +
  xlab("TCCON retrievals averaged by station and month (ppm)") +
  ylab("Level 3 XCO2 predictions averaged by station and month (ppm)") +
  theme(text = element_text(size=13)) +
  theme(legend.text=element_text(size=15), legend.title = element_blank())
    }

########################
### Part 5: Plot figures
########################

if(save_images) {

    ## Figure 10: Scatter plots
    ggsave(gScatter("oco2_xco2V8")  +  ggtitle("(b)") + theme(plot.title = element_text(size=20)),
               filename = paste0("../img/TCCON_scatterv8.png"), width = 9, height = 6)
    ggsave(gScatter("oco2_xco2V7") + theme(legend.position = "none") +  ggtitle("(a)") +
                   theme(plot.title = element_text(size=20)),
               filename = paste0("../img/TCCON_scatterv7.png"), width = 6, height = 6)

    # Scatter plot ALL comparisons V7 (not used in papers)
    gScatter_ALL <- ggplot(comp_df_sub) +
        geom_point(aes(tccon_xco2_av,oco2_xco2V7,colour=name,pch=name),size=2) +
        theme_bw() + geom_abline(slope=1) + scale_shape_manual(values=seq(0,23))

    ## All plots of predictions V7 (not used in paper)
    gpred <- ggplot(comp_df_sub) + geom_point(aes(date,tccon_xco2_av)) +
        geom_point(aes(date,oco2_xco2V7),col="red") + facet_wrap(~name) + theme_bw()

    # All plots of differences V7 (not used in paper)
    gdiffs <- ggplot(comp_df_sub) + geom_point(aes(date,tccon_xco2_av - oco2_xco2V7))  +
        facet_wrap(~name) + theme_bw()

    # Box plots of ALL differences V7 (not used in paper)
    gboxes <- ggplot(comp_df_sub) + geom_boxplot(aes(y = oco2_xco2V7 - tccon_xco2_av,x=name)) +
        theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

    # Box plots of differences by season V7 (not used in paper)
    gseasons <- ggplot(comp_df_sub) +
        geom_boxplot(aes(y = oco2_xco2V7 - tccon_xco2_av,x=name,fill=season)) +
        theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
        scale_fill_manual(values=c("blue","orange","yellow","purple"))

    ## Plot change in R2 vs latitude
    x <- (TCCON_latitudes %>% arrange(-tccon_latitude))$tccon_latitude
    g1 <- ggplot(Rows1V7df) + geom_point(aes(x,Rows1V8df$Bias - Bias, size = N)) +
        geom_hline(yintercept = 0) +
        theme_bw() + xlab("lat (deg)") + ylab("Change in Bias (V8r - V7r)")
    g2 <- ggplot(Rows1V7df) + geom_point(aes(x,Rows1V8df$MAE - MAE, size = N)) +
        geom_hline(yintercept = 0) +
        theme_bw() + xlab("lat (deg)") + ylab("Change in MAPE (V8r - V7r)")
    g3 <- ggplot(Rows1V7df) + geom_point(aes(x,Rows1V8df$RMSE - RMSE, size = N)) +
        geom_hline(yintercept = 0) +
        theme_bw() + xlab("lat (deg)") + ylab("Change in RMSPE (V8r - V7r")
    g4 <- ggplot(Rows1V7df) + geom_point(aes(x,R2 - Rows1V8df$R2, size = N)) +
        geom_hline(yintercept = 0) +
        theme_bw() + xlab("lat (deg)") + ylab("Change in R2 (V7r - V8r)")
    #gridExtra::grid.arrange(g1,g2,g3,g4,ncol = 2)



    ## Find the slope constrained to pass through the origin
    lm(comp_df_sub$oco2_xco2V7 ~ comp_df_sub$tccon_xco2_av + 0) %>% summary()
    lm(comp_df_sub$oco2_xco2V8 ~ comp_df_sub$tccon_xco2_av + 0) %>% summary()

    # Figure 9: Plot Time series, differences, and histograms, of Lamont and Wollongong
    options(digits=10)
    part_label = list(Lamont = c("(a)","(b)","(c)"),
                      Wollongong = c("(d)","(e)","(f)"))
    for(name_station in c("Lamont","Wollongong")) {
        X <- filter(comp_df_sub,name==name_station)
        Summary_results <- data.frame(x = c(2.5,2.5,2.5),
                                      y = c(0.5,0.45,0.4),
                                      txt = c(paste0("Mean = ",mean(X$oco2_xco2V7 - X$tccon_xco2_av) %>% round(2) %>% format(nsmall=2)),
                                            paste0("Median = ",median(X$oco2_xco2V7 - X$tccon_xco2_av) %>% round(2) %>% format(nsmall=2)),
                                            paste0("SD = ",sd(X$oco2_xco2V7 - X$tccon_xco2_av) %>% round(2) %>% format(nsmall=2))))

        g1 <- ggplot(X) + geom_point(aes(date,oco2_xco2V7)) +
            geom_point(aes(date,tccon_xco2_av),col="red",alpha=0.36) +
            geom_errorbar(aes(date,
                              ymin = tccon_xco2_av - 2*tccon_error_av,
                              ymax = tccon_xco2_av + 2*tccon_error_av),alpha=0.4,col="red") +
            geom_errorbar(aes(date,
                              ymin = oco2_xco2V7 - oco2_sdV7,
                              ymax = oco2_xco2V7 + oco2_sdV7),alpha=0.7,col="black") +
            theme_bw() + ggtitle(part_label[[name_station]][1]) +
            theme(text = element_text(size=20)) +
            ylab("XCO2 (ppm)") + geom_text(data = data.frame(x=date("2015-01-01"),
                                                             y=408,
                                                             txt=paste0("r = ",cor(X$oco2_xco2V7,X$tccon_xco2_av) %>% round(2) %>% format(nsmall=2))),aes(x,y,label=txt),size=7) +
            coord_cartesian(ylim = c(392,409))

        g2 <- ggplot(X) + geom_point(aes(date,oco2_xco2V7 - tccon_xco2_av))   +
            geom_errorbar(aes(date,
                              ymin = oco2_xco2V7 - tccon_xco2_av  - 2*sqrt(oco2_sdV7^2 + tccon_error_av^2),
                              ymax = oco2_xco2V7 - tccon_xco2_av  + 2*sqrt(oco2_sdV7^2 + tccon_error_av^2)),
                          alpha=0.4,col="black") +
            ylab(expression(paste(Delta," XCO2 (ppm)"))) +
            theme_bw() + ggtitle(part_label[[name_station]][2]) +
            theme(text = element_text(size=20)) +
            coord_cartesian(ylim = c(-4,4.5))

        g3 <- ggplot(X) + geom_histogram(aes(oco2_xco2V7 - tccon_xco2_av,y=..density..),
                                         col="black",bins=20) +
            theme_bw() + xlab(expression(paste(Delta," XCO2 (ppm)"))) +
            geom_text(data = Summary_results,aes(x,y,label=txt),hjust="right") +
            theme_bw() + ggtitle(part_label[[name_station]][3]) +
            theme(text = element_text(size=20)) + ylim(c(0,0.7)) + xlim(c(-3.5,3.5))

        ggsave(g1,filename=paste0("../img/",name_station,"_ts.png"),width=8,height=4)
        ggsave(g2,filename=paste0("../img/",name_station,"_diff.png"),width=8,height=4)
        ggsave(g3,filename=paste0("../img/",name_station,"_hist.png"),width=8,height=4)
    }

}

########################################################################
## UNCERTAINTY DIAGNOSTICS taking into account sample size (deprecated)
########################################################################
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
                    my_qqplot(quants,tccon_xco2_av,oco2_xco2V7,
                              sqrt(oco2_sdV7^2 + tccon_error_av^2 )))
        cum2 <- with(sub_df,
                     my_qqplot(quants,tccon_xco2_av,oco2_xco2V7 - mean(oco2_xco2V7 - tccon_xco2_av),
                               sqrt(oco2_sdV7^2 + tccon_error_av^2)))
        plot((quants - 0.5)*2,cum,asp=1,xlim=c(0,1),ylim=c(0,1),
             col="green", main=name_station, xlab = "Nominal coverage",
        ylab = "Empirical coverage")
        lines((quants - 0.5)*2,cum2,asp=1,xlim=c(0,1),ylim=c(0,1),col="blue" ,main=name_station)

        ## Effective sample size
        if(0){
            #if(nrow(sub_df) > 30) {
            sub_df$diff <- sub_df$oco2_xco2V7 - sub_df$tccon_xco2_av
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
            sim <- rnorm(mean=sub_df$oco2_xco2V7,
                         sd = sqrt(sub_df$oco2_sdV7^2 + sub_df$tccon_error_av^2),
                         n = nrow(sub_df))
            idx <- sample(1:nrow(sub_df),neff)
            cum_mat[,i] <-  my_qqplot(quants,sim[idx],sub_df$oco2_xco2V7[idx],
                                      sqrt(sub_df$oco2_sdV7^2 + sub_df$tccon_error_av^2)[idx])
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
