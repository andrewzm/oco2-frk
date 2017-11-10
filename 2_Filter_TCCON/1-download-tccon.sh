#/bin/sh

## Downloads TCCON data files into a folder called "tccon".

wget -r -c -nH -nd -np -A nc -P tccon ftp://tccon.ornl.gov/2014Public/
