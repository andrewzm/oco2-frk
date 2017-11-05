library(ncdf4)

inputfiles <- list.files(path="oco2v8", pattern="*.nc4$", full.names=TRUE, recursive=FALSE)

for (i in 1:length(inputfiles)) {
  print(paste("reading",inputfiles[i]))
  nc <- nc_open(inputfiles[i])
  oneday <- data.frame(
    datetime=strptime(format(ncvar_get(nc,'sounding_id'), scientific=FALSE), "%Y%m%d%H%M%S", tz="UTC"),
    longitude=ncvar_get(nc,'longitude'),
    latitude=ncvar_get(nc,'latitude'),
    xco2=ncvar_get(nc,'xco2'),
    xco2_uncertainty=ncvar_get(nc,'xco2_uncertainty'),
    xco2_quality_flag=ncvar_get(nc,'xco2_quality_flag'),
    warn_level=ncvar_get(nc,'warn_level'),
    operation_mode=ncvar_get(nc,'Sounding/operation_mode'),
    stringsAsFactors = FALSE
  )

  # Filter the data before saving
  oneday <- oneday[complete.cases(oneday),]
  oneday <- oneday[oneday$xco2_quality_flag==0,]
  oneday <- oneday[oneday$warn_level<15,]
  oneday <- oneday[oneday$operation_mode<2,]
  oneday <- oneday[(oneday$xco2_uncertainty < 3),]
  oneday$xco2_uncertainty <- pmax(oneday$xco2_uncertainty,2)

  if (i==1) {
    # If this is the first line, write column names to a new file.
    output_names <- TRUE
    output_append <- FALSE
  } else {
    output_names <- FALSE
    output_append <- TRUE
  }

  write.table(data.frame("day"=oneday$datetime, "lon"=oneday$longitude, "lat"=oneday$latitude, "xco2"=oneday$xco2, "std"=oneday$xco2_uncertainty),
              file='oco2lite.csv', row.names=FALSE, col.names=output_names, sep=',', append=output_append)
}
