First open 1-download.sh, read the instructions carefully and add your Earthdata login details, then begin downloading the input data:

$ sh 1-download.sh

Read the data into R, filter out unwanted readings, and save to a CSV file:

$ Rscript 2-filter.R

Run Fixed Rank Kriging on the filtered data to produce daily level 3 files. The amount of RAM required here is approximately 20 GB with these options and input data. If you have 64 GB RAM you can run three instances of this script in parallel, etc.

$ Rscript 3-frk.R

Produce plots from the level 3 files. You can run multiple instances of this in parallel too.

$ Rscript 4-plot.R