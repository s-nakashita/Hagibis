#!/bin/sh
if [ $# -lt 2 ]; then
    echo "Usage :: $0 yyyymmddhh center"
    exit
fi

yyyymmddhh=${1}
yyyy=${yyyymmddhh:0:4}
dcolat=${dcolat:-6.0}

center=${2}
case $center in
    ecmwf ) orig=ecmf ;;
    jma ) orig=rjtd ;;
    ncep ) orig=kwbc ;;
    ukmo ) orig=egrr ;;
esac
#orig=${center}
#trackf=../pytrack/$orig/track${yyyymmddhh}.txt
trackf=../pytrack/$center/gtrack${yyyymmddhh}_mean.txt
#trackf=../pytrack/$orig/track_anl.txt
#datadir=../../netcdf/tigge/${yyyy}/$center
datadir=./
outdir=../../netcdf/tigge/${yyyy}/$center/rotate
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

#echo $yyyymmddhh $datadir $trackf $nlon $nlat $latmax | python rotate_scalar.py

python inv_rotate_1d.py $yyyymmddhh $datadir $trackf $dcolat

#echo $yyyymmddhh | python create_netcdf.py

#mv np_${yyyymmddhh}_mean.nc ${outdir}/
#rm np_*_${yyyymmddhh}_mean.nc
