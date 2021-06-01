#!/bin/sh
#if [ $# -lt 2 ]; then
#    echo "Usage :: $0 yyyymmddhh center"
#    exit
#fi

yyyymmddhh=2019100912
yyyy=${yyyymmddhh:0:4}
mem=${1}
nlon=${nlon:-100}
nlat=${nlat:-50}
latmax=${latmax:-8}

center=${2}
case $center in
    ecmwf ) orig=ecmf ;;
    jma ) orig=rjtd ;;
    ncep ) orig=kwbc ;;
    ukmo ) orig=egrr ;;
esac

#trackf=../pytrack/track2019100900_gsm.txt
#datadir=../../netcdf/gsm/gl/2019/10
#outdir=../../netcdf/gsm/gl/2019/10/rotate
trackf=../pytrack/$center/track${yyyymmddhh}_${mem}.txt
datadir=../../netcdf/tigge/${yyyy}/$center
outdir=../../netcdf/tigge/${yyyy}/$center/rotate
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

echo $datadir $trackf $nlon $nlat $latmax $yyyymmddhh $mem | python rotate_vort.py

#echo $yyyymmddhh $datadir $trackf $nlon $nlat $latmax | python rotate_uv.py

#python create_netcdf_v.py

mv np_*.nc ${outdir}/
#rm np_*_${yyyymmddhh}_mean.nc
