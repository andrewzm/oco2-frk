library(FRK)
library(dplyr)
library(spacetime)
library(sp)

BAUs_to_df <- function (selecteddate, BAUs)
{
  print("Converting BAUs to data frame")

  sp_polys <- BAUs[,which(time(BAUs) == as.POSIXct(selecteddate))]

  vars = c("mu","sd")

  polynames <- as.character(row.names(sp_polys))
  list_polys <- lapply(1:length(sp_polys), function(i) {
    coords <- sp_polys@polygons[[i]]@Polygons[[1]]@labpt
    poldf <- data.frame(id = polynames[i], lon = coords[1], lat = coords[2], stringsAsFactors = FALSE)
    poldf
  })
  df_polys <- bind_rows(list_polys)
  df_polys$id <- as.character(df_polys$id)
  sp_polys$id <- row.names(sp_polys)
  cnames <- coordnames(sp_polys)
  vars_no_coords <- vars[which(!vars %in% cnames)]
  if (length(vars_no_coords) > 0)
    df_polys <- left_join(df_polys, sp_polys@data[c("id", vars_no_coords)], by = "id")
  df_polys$id <- NULL
  df_polys$date <- selecteddate
  df_polys
}

computeBAUs <- function(selecteddate, sixteendays) {
  print("Computing BAUs")

  opts_FRK$set("progress",TRUE)
  opts_FRK$set("parallel",8L)
  lo_res <- 0

  space <- sixteendays[,c("lon","lat")]
  coordinates(space) = ~lon+lat
  proj4string(space)=CRS("+proj=longlat +ellps=sphere")
  STobj <- STIDF(space,sixteendays$day,data=sixteendays)

  BAUs <- auto_BAUs(manifold=STsphere(), data=STobj, type="grid", tunit="days", cellsize=c(1,1,1), xlim=c(-180,180), ylim=c(round(min(sixteendays$lat),0),round(max(sixteendays$lat),0)))

  BAUs$fs = 1
  G_spatial <- auto_basis(manifold = sphere(), data=as(STobj,"Spatial"), nres = ifelse(lo_res,1,3), prune=0, type = "bisquare", subsamp = 20000, isea3h_lo = 1)
  G_temporal <- local_basis(manifold=real_line(), loc = matrix(seq(1,16,by=2)), scale = rep(4,8),type="bisquare")
  G_spacetime <- TensorP(G_spatial, G_temporal)

  print("Binning")
  S <- SRE(xco2 ~ lat +1, list(STobj), G_spacetime, BAUs = BAUs, est_error = FALSE)

  print("Fitting")
  S <- SRE.fit(S, n_EM = 10, tol = 0.01, print_lik=FALSE)

  print("Predicting")
  BAUs <- SRE.predict(SRE_model = S, obs_fs=FALSE)

  return(BAUs)
}

oco2lite <- read.csv('oco2lite.csv')
oco2lite$day <- as.Date(oco2lite$day, tz="UTC")

if (!dir.exists("baus")) {
  dir.create("baus")
}
if (!dir.exists("level3")) {
  dir.create("level3")
}

for (selecteddate in as.character(seq(min(oco2lite$day), max(oco2lite$day), "days"))) {

  if ( file.exists(file.path("level3",paste0(selecteddate,".csv"))) ) {
    # This date has already been processed
    next
  }
  file.create(file.path("level3",paste0(selecteddate,".csv")))

  print(selecteddate)

  startdate <- as.Date(selecteddate,tz="UTC")-7
  enddate <- as.Date(selecteddate,tz="UTC")+8

  sixteendays <- oco2lite[oco2lite$day >= startdate & oco2lite$day <= enddate,]

  # Don't compute BAUs if there is not data on both sides of the selected date within the date window
  if ( is.null(sixteendays) | selecteddate > max(sixteendays$day) | selecteddate < min(sixteendays$day) ) {
    print("Skip this date")
  } else {
    if ( file.exists(file.path("baus",paste0(selecteddate,".RData"))) ) {
      print("Loading BAUs")
      load(file.path("baus",paste0(selecteddate,".RData")))
    } else {
      BAUs <- computeBAUs(selecteddate, sixteendays)
      print("Saving BAUs")
      save(BAUs, file=file.path("baus",paste0(selecteddate,".RData")))
    }

      Level3 <- BAUs_to_df(selecteddate, BAUs)
      print("Saving Level 3 data")
      write.csv(Level3, file=file.path("level3",paste0(selecteddate,".csv")), row.names=FALSE)
  }
}
