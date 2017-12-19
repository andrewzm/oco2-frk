#/bin/sh

## Replication instructions for generating the OCO-2 data file OCO2data.Rdata
## ==========================================================================
##
## The following is ideally done on a Linux machine since UNIX commands are needed to download the data.
##
## 1. Please visit https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Download%20Data%20Files%20from%20HTTP%20Service%20with%20wget
## and carefully follow Steps 1 and 2.
##
## 2. In a terminal go to the project folder and edit this script below by replacing USER and PASS with your Earthdata login details from this link: https://disc.gsfc.nasa.gov/data-access
##
## 3. Run this script.
##
## If there is a certificate error, please add the flag --no-check-certificate. The "wget" will download the data and put them into folders named
## "oco2v7" and "oco2v8". When the download is complete the folder will contain several .nc4 data files.

wget --user=USER --password=PASS -r -c -nH -nd -np -A nc4 -P oco2v7 http://oco2.gesdisc.eosdis.nasa.gov/data/s4pa/OCO2_DATA/OCO2_L2_Lite_FP.7r/
wget --user=USER --password=PASS -r -c -nH -nd -np -A nc4 -P oco2v8 http://oco2.gesdisc.eosdis.nasa.gov/data/s4pa/OCO2_DATA/OCO2_L2_Lite_FP.8r/
