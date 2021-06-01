#!/bin/sh
if [ $# -lt 2 ]; then
    echo "Usage :: $0 yyyymmddhh center"
    exit
fi

yyyymmddhh=${1}
yyyy=${yyyymmddhh:0:4}
nlon=${nlon:-360}
nlat=${nlat:-181}
latmax=${latmax:-90}

center=${2}
case $center in
    ecmwf ) orig=ecmf ;;
    jma ) orig=rjtd ;;
    ncep ) orig=kwbc ;;
    ukmo ) orig=egrr ;;
esac
#orig=${center}
#trackf=../pytrack/$orig/track${yyyymmddhh}.txt
#trackf=../pytrack/$center/track${yyyymmddhh}_mean.txt
trackf=../pytrack/$center/gtrack${yyyymmddhh}_mean.txt
#trackf=../pytrack/$orig/track_anl.txt
datadir=../../netcdf/tigge/${yyyy}/$center
outdir=../../netcdf/tigge/${yyyy}/$center/rotate
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

echo $yyyymmddhh $datadir $trackf $nlon $nlat $latmax | python rotate_scalar_glb.py

echo $yyyymmddhh $datadir $trackf $nlon $nlat $latmax | python rotate_uv_glb.py

echo $yyyymmddhh | python create_netcdf_glb.py

#mv np_${yyyymmddhh}_mean.nc ${outdir}/
#rm np_*_${yyyymmddhh}_mean.nc
