#!/bin/sh -f
infile=/spare11/atm350/el381212/finalproject/2020121500_nam211.gem
outfile=/spare11/atm350/el381212/finalproject/2020121500_nam211.nc
java -Xmx1g -classpath /kt11/ktyle/software/netcdf-java/netcdfAll-5.5.2.jar ucar.nc2.dataset.NetcdfDataset -in $infile -out $outfile --netcdf4 -isLargeFile
 
exit 0
