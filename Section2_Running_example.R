#######################################################################
### 2D spatial example (Section 2)
#######################################################################

## Load libraries
library("dplyr")
library("FRK")
library("gstat")
library("ggplot2")
library("grid")
library("gridExtra")
library("Matrix")
library("RandomFields")      
library("sp")

save_images = 0  # set to 1 to save images


###################### UTILITY FUNCTIONS  #################################

## Squared prediction error
## z:         data
## pred:      prediction
## root_mean: if true, the RMSPE is returned instead of the SPE
##
## Returns: RMSPE if root_mean == TRUE, otherwise vector of squared prediction errors
spe <- function(z,pred,root_mean = TRUE) {
    Y <- (z - pred)^2
    if(root_mean)
        Y <- sqrt(mean(Y))
    else
        Y
}

## Prediction error
## z:    data
## pred: prediction
## mean: if true, the mean prediction error is returned
##
## Returns: MPE if root_mean == TRUE, otherwise vector of prediction errors
pe <- function(z,pred,mean = TRUE) {
    Y <- (z - pred)
    if(mean)
        Y <- mean(Y)
    else
        Y
}

## 90% coverage
## z:       data
## mu:      prediction
## se:      prediction error
## average: if true, returns the average coverage
##
## Returns: 90% coverage if average == TRUE, otherwise vector of TRUE/FALSE indicators (sample within
##          90% prediction interval or not)
coverage90 <- function(z,mu,se,average=TRUE) {
    lower <- mu - 1.64*se
    upper <- mu + 1.64*se
    sum((z < upper) & (z > lower)) / length(z)
    if(average)
        sum((z < upper) & (z > lower)) / length(z)
    else
        (z < upper) & (z > lower)
}

N <-                1000  # number of data points
l <-                0.15  # length scale in exp. covariance function
SNR <-              1     # signal-to-noise ratio

process_var <-      1  # process marginal variance = 1
measurement_var <-  process_var/SNR  # find meas. variance based on SNR and proc. variance

## Construct a 100 x 100 grid
xaxis <- seq(0, 1, length=100)
yaxis <- seq(0, 1, length=100)

## Simulate a field on this grid (at this stage just to generate indices)
z <- RFsimulate(RMexp(var = process_var, scale = l), 
                x = xaxis, y = yaxis)

## Extact the 1e6 coordinates
zcoords <- coordinates(z)

## Form a long data frame based on these data
z_df <- data.frame(x = zcoords[, 1],
                   y = zcoords[, 2],
                   grid_box = 1:nrow(zcoords))

## Make sure the locations are fixed for reproducibility purposes by fixing seed
set.seed(1)

## Sample N observation locations
d <- d_df <- sample_n(z_df, N, replace = FALSE) %>% dplyr::select(x,y,grid_box)

## Cast d to an sp object
coordinates(d) = ~x + y

## Save all observation indices (w.r.t. the 100 x 100 grid) in OBS
OBS <- d$grid_box

## Now generate indices for th unobserved locations
UNOBS <- setdiff(z_df$grid_box,OBS) 

## Now generage some random indices for diagnostics
DIAG <- sample(UNOBS, 200)

## The prediction locations are the union of the observed and unobserved
pred_idx <- union(OBS, UNOBS)

## Simulate an exponential field with sim parameters
z <- RFsimulate(RMexp(var = process_var, scale=l),
                x = xaxis, y = yaxis)

## Relabel to z
z_df$z <- z$variable1

## Make an sp version of the data for gstat and FRK
z.sp <- as(z, "SpatialPixelsDataFrame")

## Re-assign coordinate names since the above as() function does not keep coordnames
coordnames(z.sp) <- c("x", "y")

## Add a fine-scale variation field (prop. to the identity) for FRK
z.sp$fs <- 1

## Find which simulated locations correspond with the data locations
d_df$z_true <- d$z_true <- over(d,z)[, 1]

## The data is just the process + measuremnet error. Recall that d is of object sp
d_df$z <- d$z <- d$z_true + sqrt(measurement_var) * rnorm(N)

## Do the same but this time for the noisy case
d_df$z_10SNR <- d$z_10SNR <- d$z_true + sqrt(10*measurement_var) * rnorm(N)

################
#### EXAMPLE 1
################

## Create variogram model (also for low SNR case)
z.model = vgm(process_var, model = "Exp", 
              nugget=0, range = l, Err = measurement_var)
z.model_10SNR = vgm(process_var, model = "Exp",
                    nugget=0, range = l, Err = 10*measurement_var)

## Convert z.sp2 to SpatialPixelsDataFrame
z.sp2 <- as(z,"SpatialPointsDataFrame")
gridded(z.sp2) = TRUE

## Do kriging with the true parameters (gold standard)
z.kriged = krige(z~1, d, z.sp2[pred_idx,], model = z.model, beta=0)

## Now extract the predictions and prediction errors at all prediction indices
z_df$gstat_pred <- z_df$gstat_se <- z_df$gstat_se_obs <- NA
z_df$gstat_pred[pred_idx] <- z.kriged$var1.pred
z_df$gstat_se[pred_idx] <- sqrt(z.kriged$var1.var)
z_df$gstat_se_obs[pred_idx] <- sqrt(z.kriged$var1.var + measurement_var)

## Figure 1
my_colours <- c("#03006d","#02008f","#0000b6","#0001ef",
                "#0000f6","#0428f6","#0b53f7","#0f81f3",
                "#18b1f5","#1ff0f7","#27fada","#3efaa3",
                "#5dfc7b","#85fd4e","#aefc2a","#e9fc0d",
                "#f6da0c","#f5a009","#f6780a","#f34a09",
                "#f2210a","#f50008","#d90009","#a80109")

gbase <- ggplot() +  coord_fixed(xlim = c(0,1),ylim = c(0,1)) + 
    xlab(expression(s[1])) + ylab(expression(s[2]))

gexp <- (gbase + geom_point(data = z_df,aes(x,y),pch=46) + theme_bw() +
             geom_point(data=z_df[DIAG,],aes(x,y),pch=3,col="red") + 
             geom_point(data=z_df[OBS,],aes(x,y),pch=1,col="blue") +
             ggtitle("(a)") + theme(text = element_text(size=20)))
if(save_images) 
  ggsave(plot = gexp,filename="img/E1_exp_setup.png",width=5,height=5)

gfield <- gbase + geom_tile(data=z_df,aes(x,y,fill=z)) + 
    scale_fill_gradientn(colours=my_colours, name=expression(Y),limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(b)") + theme(text = element_text(size=20))
if(save_images) 
  ggsave(plot = gfield,filename="./img/E1_field.png",width=5,height=5)

gpred <- gbase + geom_tile(data=z_df,aes(x,y,fill=gstat_pred)) + 
    scale_fill_gradientn(colours=my_colours, name="pred",limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(d)") + theme(text = element_text(size=20))
if(save_images) 
  ggsave(plot = gpred,filename="./img/E1_pred.png",width=5,height=5)

gse <- gbase + geom_tile(data = z_df, aes(x,y,fill=gstat_se)) + 
    scale_fill_distiller(palette="BrBG",name="s.e.") +
    theme_bw() + ggtitle("(c)") + theme(text = element_text(size=20))
if(save_images) 
  ggsave(plot = gse,filename="./img/E1_se.png",width=5,height=5)

MSPE <- spe(z = z_df[DIAG,]$z,pred = z_df[DIAG,]$gstat_pred,
            root_mean = TRUE)
pred_obs <- z_df[DIAG,]$z + rnorm(length(DIAG),mean = 0,
                                  sd=sqrt(measurement_var))
cov90 <- coverage90(z = z_df[DIAG,]$z,mu = z_df[DIAG,]$gstat_pred,
                    se = z_df[DIAG,]$gstat_se)
cov90_obs <- coverage90(z = pred_obs, mu = z_df[DIAG,]$gstat_pred,
                        se = z_df[DIAG,]$gstat_se_obs)
cov90_wrong <- coverage90(z = pred_obs,mu = z_df[DIAG,]$gstat_pred,
                          se = z_df[DIAG,]$gstat_se)
cov90_obs_wrong <- coverage90(z =  z_df[DIAG,]$z, mu = z_df[DIAG,]$gstat_pred,
                              se = z_df[DIAG,]$gstat_se_obs)
rbind(cbind(cov90,cov90_obs_wrong),cbind(cov90_wrong,cov90_obs))

################
#### EXAMPLE 1
################

## Create matrices for simple kriging
Dpp <- fields::rdist(z_df[,c("x","y")])
Doo <- fields::rdist(z_df[OBS,c("x","y")])
Dop <- fields::rdist(z_df[OBS,c("x","y")],z_df[,c("x","y")])
Dpo <- t(Dop)
Cxpxp <- process_var*exp(-Dpp/l)
Cxoxo <- process_var*exp(-Doo/l)
Cxoxp <- process_var*exp(-Dop/l)
Cxpxo <- process_var*exp(-Dpo/l)
B <- chol2inv(chol(Cxoxo + measurement_var*diag(N)))
MU <- Cxpxo %*% B %*% z_df[OBS,"z"]
COV <- Cxpxp - Cxpxo %*% B %*% Cxoxp
L <- t(chol(COV))

## Do the conditional simulations
z_df$Sim1 <- as.numeric(MU + L %*% rnorm(nrow(z_df)))
z_df$Sim2 <- as.numeric(MU + as.numeric(L %*% rnorm(nrow(z_df))))

## Figure 2
gSim1 <- gbase + geom_tile(data=z_df,aes(x,y,fill=pmax(pmin(Sim1,3),-3.2))) + 
    scale_fill_gradientn(colours=my_colours, name=expression(Y[sim1]),limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(a)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gSim1,filename="./img/E1_Sim1.png",width=5,height=5)

gSim2 <- gbase + geom_tile(data=z_df,aes(x,y,fill=pmax(pmin(Sim2,3),-3.2))) + 
    scale_fill_gradientn(colours=my_colours, name=expression(Y[sim2]),limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(c)") + theme(text = element_text(size=20))

if(save_images)
  ggsave(plot = gSim2,filename="./img/E1_Sim2.png",width=5,height=5)

################
#### EXAMPLE 3
################

f <- z ~ 1                      # FRK formula
d$std <- sqrt(measurement_var)  # Add measurement error to spatial data

## Now run FRK with three resolutions
S <- FRK(f = f,
         data = list(d),
         BAUs = z.sp[pred_idx, ],
         nres = 3)

## And predict using the returned SRE model at all BAUs
Pred <- SRE.predict(SRE_model = S,
                    obs_fs = FALSE)

## Now extract the predictions and prediction errors at all locations
z_df$FRK_pred <- z_df$FRK_se <- z_df$FRK_se_obs <- NA
z_df$FRK_pred[pred_idx] <- Pred@data["mu"][,1]
z_df$FRK_se[pred_idx] <- Pred@data["sd"][,1]
z_df$FRK_se_obs[pred_idx] <- sqrt(Pred@data["sd"][,1]^2 + S@Ve[1,1])

## Figure 3
gpredFRK <- gbase + geom_tile(data=z_df,aes(x,y,fill=FRK_pred)) + 
    scale_fill_gradientn(colours=my_colours, name="pred",limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(b)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gpredFRK,filename="./img/E2_pred.png",width=5,height=5)

gseFRK <- gbase + geom_tile(data = z_df, aes(x,y,fill=pmax(pmin(FRK_se,0.7),0.4))) + 
    scale_fill_distiller(palette="BrBG",name="s.e.",limits=c(0.4,0.7)) +
    theme_bw() + ggtitle("(c)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gseFRK,filename="./img/E2_se.png",width=5,height=5)

MSPE_FRK <- spe(z = z_df[DIAG,]$z,pred = z_df[DIAG,]$FRK_pred,
                root_mean = TRUE)
pred_obs_FRK <- z_df[DIAG,]$z + rnorm(length(DIAG),mean = 0,
                                      sd=sqrt(measurement_var))
cov90_FRK <- coverage90(z = z_df[DIAG,]$z,mu = z_df[DIAG,]$FRK_pred,
                        se = z_df[DIAG,]$FRK_se)
cov90_obs_FRK <- coverage90(z = pred_obs, mu = z_df[DIAG,]$FRK_pred,
                            se = z_df[DIAG,]$FRK_se_obs)
cov90_wrong <- coverage90(z = pred_obs,mu = z_df[DIAG,]$gstat_pred,
                          se = z_df[DIAG,]$gstat_se)
cov90_obs_wrong <- coverage90(z =  z_df[DIAG,]$z, mu = z_df[DIAG,]$gstat_pred,
                              se = z_df[DIAG,]$gstat_se_obs)

## CLear memory
rm(L); rm(Dpp); rm(Dpo); rm(Dop); rm(Doo)
rm(COV); gc()

## Conditional simulations from FRK
KL <- as(t(chol(S@Khat)),"matrix")
S_BAUs <- eval_basis(S@basis,z.sp) %>% as("matrix")
S_obs <- eval_basis(S@basis,d) %>% as("matrix")

xx <- sqrt(rowSums((S_BAUs) * S_BAUs))  
xx <- xx + 1*(xx == 0)          
S_BAUs <- S_BAUs / (as.numeric(xx))
S_obs <- S_BAUs[OBS,]

Cxpxp <- tcrossprod(S_BAUs %*% KL) + S@sigma2fshat*diag(nrow(z_df))
Cxoxo <- tcrossprod(S_obs %*% KL) + S@sigma2fshat*diag(N)
Cxoxp <- (S_obs %*% S@Khat %*% t(S_BAUs))
Cxpxo <- t(Cxoxp)
B <- chol2inv(chol(Cxoxo + measurement_var*diag(N)))
BL <- t(chol(B))

MU <- Cxpxo %*% (B %*% z_df[OBS,"z"])
COV <- Cxpxp - tcrossprod(Cxpxo %*% BL)
L <- t(chol(COV))

## Generate the conditional simulations
z_df$Sim1FRK <- as.numeric(MU + (L %*% rnorm(nrow(z_df))))
z_df$Sim2FRK <- as.numeric(MU + (L %*% rnorm(nrow(z_df))))

## Figure 4 (a) and (b)
gSim1FRK <- gbase + geom_tile(data=z_df,aes(x,y,fill=pmax(pmin(Sim1FRK,3),-3.2))) + 
    scale_fill_gradientn(colours=my_colours, 
                         name=expression(Y[sim1]),limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(a)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gSim1FRK,filename="./img/E2_Sim1FRK.png",width=5,height=5)

gSim2FRK <- gbase + geom_tile(data=z_df,aes(x,y,fill=pmax(pmin(Sim2FRK,3),-3.2))) + 
    scale_fill_gradientn(colours=my_colours, 
                         name=expression(Y[sim2]),limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(b)") + theme(text = element_text(size=20))

if(save_image)
  ggsave(plot = gSim2FRK,filename="./img/E2_Sim2FRK.png",width=5,height=5)

## Figure 4 (c) -- the covariance function
Dpp <- fields::rdist(z_df[,c("x","y")])
xvals <- Dpp[5050,]                         # Consider the middle of the domain
xvals[nrow(Dpp)+1] <- 0.0001
orderx <- order(xvals)
xvals <- xvals[orderx]

yvals <- (Cxpxp + S@sigma2fshat*diag(nrow(z_df)))[5050,]
yvals[nrow(Dpp)+1] <- yvals[orderx][3]
yvals <- yvals[orderx]

cov_df1 <- data.frame(h = seq(0,1,by=0.01), 
                      C = process_var*exp(-seq(0,1,by=0.01)/l))
cov_df2 <- data.frame(h = xvals, C = yvals)

gcov <- ggplot() + geom_line(data=cov_df1,aes(h,C),col="red",linetype=2) +
    geom_line(data=cov_df2,aes(h,C)) + theme_bw() + xlim(c(0,0.5)) + 
    ggtitle("(c)") + theme(text = element_text(size=20)) + ylab("C(h)")
if(save_images)
  ggsave(plot = gcov,filename="./img/E2_cov_fns.png",width=5,height=5)

################
#### EXAMPLE 4
################

## FIXED WINDOW
bins <- seq(0,1,by=0.1)
bins[11] <- 1.01
z.sp2$yc <- z.sp2@coords[,2]
z.sp2$xc <- z.sp2@coords[,1]
unique_ys <- unique(z.sp2$yc)
krig_df <- NULL
data_df <- NULL
z_df$bin <- 0
for(i in 10:1) { # for each window do spatial-only kriging
    coords <- data.frame(coordinates(d))
    bc <-  (bins[i+1] + bins[i])/2
    idx <- which(coords$y >= bins[i] & coords$y < bins[i+1])
    coords1D <- cbind(coords[idx,1],bc + runif(length(idx),min=-0.001,max=0.001))
    colnames(coords1D) <- c("x","y")
    dsub <- SpatialPointsDataFrame(coords = coords1D,data = d@data[idx,])
    
    idx2 <- which(z_df$y >= bins[i] & z_df$y < bins[i+1])
    z_df$bin[idx2] <- i
    
    closest_y <- unique_ys[which.min(abs(unique_ys - bc))]
    z.sp2_sub <- subset(z.sp2,yc == closest_y)
    z.kriged = krige(z~1, dsub, z.sp2_sub, model = z.model,beta=0) 
    krig_df <- rbind(krig_df,cbind(data.frame(bin=i,x = z.sp2_sub$xc,
                                              z_true = z.sp2_sub$variable1),
                                   z.kriged@data))
    data_df <- rbind(data_df,cbind(data.frame(bin = i, x = dsub$x),dsub@data))
}
krig_df_full <- krig_df
krig_df_full$z_true <- NULL
krig_df_full <- rename(krig_df_full, Spat_pred = var1.pred, Spat_var = var1.var) 
z_df <- left_join(z_df,krig_df_full,by = c("bin","x"))

## Figure 5 (a)
g_Spatial <- ggplot() + geom_point(data = data_df,aes(x,z),alpha=0.3,col="blue") +
    geom_line(data = krig_df,aes(x,z_true)) + facet_wrap(~bin,nrow=2) + 
    geom_line(data = krig_df,aes(x,var1.pred),col="red",linetype=1,size=1) + 
    xlab(expression(s[1])) + ylab("Y") + ggtitle("(a)                                   Temporal bin") + theme_bw() +
    theme(text = element_text(size=20)) + 
    scale_x_continuous(breaks=c(0,0.5))
if(save_images)
  ggsave(plot = g_Spatial,filename="./img/E3_pred.png",width=10,height=5)

MSPE_Spatial <- spe(z = z_df[DIAG,]$z,pred = z_df[DIAG,]$Spat_pred,
                    root_mean = TRUE)
cov90_Spatial <- coverage90(z = z_df[DIAG,]$z,mu = z_df[DIAG,]$Spat_pred,
                            se = sqrt(z_df[DIAG,]$Spat_var))

## MOVING WINDOW
krig_df <- NULL
data_df <- NULL
for(this_y in unique(z_df$y)) { # For each point in s_2
    coords <- data.frame(coordinates(d))
    idx <- which(coords$y < (this_y + 0.1) & coords$y >= (this_y - 0.1))
    coords1D <- cbind(coords[idx,1],coords[idx,2])
    colnames(coords1D) <- c("x","y")
    dsub <- SpatialPointsDataFrame(coords = coords1D,data = d@data[idx,])
    
    z.sp2_sub <- subset(z.sp2,yc == this_y)
    z.kriged = krige(z~1, dsub, z.sp2_sub, model = z.model,beta=0) 
    krig_df <- rbind(krig_df,cbind(data.frame(y = this_y,x = z.sp2_sub$xc),
                                   z.kriged@data))
    data_df <- rbind(data_df,cbind(data.frame(y = this_y, x = dsub$x),dsub@data))
}
krig_df_full <- krig_df
krig_df_full <- rename(krig_df_full, Spat2_pred = var1.pred, Spat2_var = var1.var) 
z_df <- left_join(z_df,krig_df_full,by = c("y","x"))

MSPE_Spatial2 <- spe(z = z_df[DIAG,]$z,pred = z_df[DIAG,]$Spat2_pred,
                     root_mean = TRUE)
cov90_Spatial <- coverage90(z = z_df[DIAG,]$z,mu = z_df[DIAG,]$Spat2_pred,
                            se = sqrt(z_df[DIAG,]$Spat2_var))

## Figure 5 (b)
gpred_Spat <- gbase + geom_tile(data=z_df,aes(x,y,fill=Spat_pred)) + 
    scale_fill_gradientn(colours=my_colours, name="pred",limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(b)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gpred_Spat,filename="./img/E3_pred_spat.png",width=5,height=5)

## Figure 5 (c)
gpred_Spat2 <- gbase + geom_tile(data=z_df,aes(x,y,fill=Spat2_pred)) + 
    scale_fill_gradientn(colours=my_colours, name="pred",limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(c)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gpred_Spat2,filename="./img/E3_pred_spat2.png",width=5,height=5)

################
#### EXAMPLE 5
################

## MOVING WINDOW (LOW SNR)
krig_df <- NULL
data_df <- NULL
for(this_y in unique(z_df$y)) {
    coords <- data.frame(coordinates(d))
    idx <- which(coords$y < (this_y + 0.1) & coords$y >= (this_y - 0.1))
    coords1D <- cbind(coords[idx,1],coords[idx,2])
    colnames(coords1D) <- c("x","y")
    dsub <- SpatialPointsDataFrame(coords = coords1D,data = d@data[idx,])
    
    z.sp2_sub <- subset(z.sp2,yc == this_y)
    z.kriged = krige(z_10SNR~1, dsub, z.sp2_sub, model = z.model_10SNR,beta=0) 
    krig_df <- rbind(krig_df,cbind(data.frame(y = this_y,x = z.sp2_sub$xc),
                                   z.kriged@data))
    data_df <- rbind(data_df,cbind(data.frame(y = this_y, x = dsub$x),dsub@data))
}
krig_df_full <- krig_df
krig_df_full <- rename(krig_df_full, Spat2_10SNR_pred = var1.pred, 
                       Spat2_10SNR_var = var1.var) 
z_df <- left_join(z_df,krig_df_full,by = c("y","x"))

MSPE_Spatial3 <- spe(z = z_df[DIAG,]$z,pred = z_df[DIAG,]$Spat2_10SNR_pred,
                     root_mean = TRUE)
cov90_Spatial <- coverage90(z = z_df[DIAG,]$z,mu = z_df[DIAG,]$Spat2_10SNR_pred,
                            se = sqrt(z_df[DIAG,]$Spat2_10SNR_var))

## Figure 6 (b)
gpred_10SNR_Spat <- gbase + geom_tile(data=z_df,aes(x,y,fill=Spat2_10SNR_pred)) + 
    scale_fill_gradientn(colours=my_colours, name="pred",limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(b)") + theme(text = element_text(size=20))
ggsave(plot = gpred_10SNR_Spat,filename="./img/E4_pred_10SNR_MW.png",width=5,height=5,dpi=200)

## EXACT KRIGING (LOW SNR)
z.kriged_10SNR = krige(z_10SNR~1, d, z.sp2[pred_idx,], model = z.model_10SNR,beta=0)
z_df$gstat_pred_10SNR_exact[pred_idx] <- z.kriged_10SNR$var1.pred

## Figure 6 (a)
gpred_10SNR_Spat_exact <- gbase + geom_tile(data=z_df,aes(x,y,fill=gstat_pred_10SNR_exact)) + 
    scale_fill_gradientn(colours=my_colours, name="pred",limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(a)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gpred_10SNR_Spat_exact,filename="./img/E4_pred_10SNR_exact.png",
         width=5,height=5)
MSPE_Spatial4 <- spe(z = z_df[DIAG,]$z,pred = z_df[DIAG,]$gstat_pred_10SNR_exact,
                     root_mean = TRUE)

## FRK (LOW SNR)
S_10_SNR <- S
idx <- match(d$z,as.numeric(S@Z))
S_10_SNR@Z[idx] <- d$z_10SNR
S_10_SNR@Ve <- S_10_SNR@Ve*10
S_10_SNR <- FRK:::.SRE.Estep.ind(S_10_SNR)

Pred <- SRE.predict(SRE_model = S_10_SNR,
                    obs_fs = FALSE)

z_df$FRK_pred_10SNR <- z_df$FRK_se_10SNR <- z_df$FRK_se_obs_10SNR <- NA
z_df$FRK_pred_10SNR[pred_idx] <- Pred@data["mu"][,1]
z_df$FRK_se_10SNR[pred_idx] <- Pred@data["sd"][,1]
z_df$FRK_se_obs_10SNR[pred_idx] <- sqrt(Pred@data["sd"][,1]^2 + S@Ve[1,1])

## Figure 6 (c)
gpredFRK_10SNR <- gbase + geom_tile(data=z_df,aes(x,y,fill=FRK_pred_10SNR)) + 
    scale_fill_gradientn(colours=my_colours, name="pred",limits=c(-3.2,3)) + 
    theme_bw() + ggtitle("(c)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gpredFRK_10SNR,filename="./img/E4_FRK_pred_10SNR.png",width=5,height=5)

gseFRK_10SNR <- gbase + geom_tile(data = z_df, aes(x,y,fill=pmax(pmin(FRK_se_10SNR,0.7),0.4))) + 
    scale_fill_distiller(palette="BrBG",name="s.e.",limits=c(0.4,0.7)) +
    theme_bw() + ggtitle("(c)") + theme(text = element_text(size=20))
if(save_images)
  ggsave(plot = gseFRK_10SNR,filename="./img/E4_FRK_se_10SNR.png",width=5,height=5)

MSPE_FRK2 <- spe(z = z_df[DIAG,]$z,pred = z_df[DIAG,]$FRK_pred_10SNR,
                 root_mean = TRUE)
cov90_FRK <- coverage90(z = z_df[DIAG,]$z,mu = z_df[DIAG,]$FRK_pred_10SNR,
                        se = z_df[DIAG,]$FRK_se_10SNR)