#!/bin/sh
orig=jma
orig=${1}
case $orig in
    ecmwf ) mem=50 ;;
    jma ) mem=26 ;;
    ncep ) mem=20 ;;
    ukmo ) mem=17 ;;
esac
idate=2019100912
edate=2019101212
smode=2
emode=${smode}
cat > sens.nml << EOF
&sens_nml
orig = "${orig}",
mem = ${mem},
idate = ${idate},
edate = ${edate},
smode = ${smode},
emode = ${emode}
/
EOF
cat sens.nml
ntype=TE
ntype=${2}
ofile=ensvsa-${ntype}-m${emode}-${orig}-${idate}_n-gr
if [ -e $ofile ]; then
    rm ${ofile}
fi
if [ read_netcdf.f90 -nt read_netcdf.mod ]; then
    make clean
fi
make ${ntype}
./ensvsa-${ntype}
./grads-ensvsa-${ntype}

cfile=ensvsa-${ntype}-m${emode}-${orig}-${idate}_n.ctl
ncfile=ensvsa-${ntype}-m${emode}-${orig}-${idate}_n.nc
yyyy=${idate:0:4}
mm=${idate:4:2}
dd=${idate:6:2}
hh=${idate:8:2}
echo $yyyy $mm $dd $hh
case $mm in
    01 ) mon=Jan ;;
    02 ) mon=Feb ;;
    03 ) mon=Mar ;;
    04 ) mon=Apr ;;
    05 ) mon=May ;;
    06 ) mon=Jun ;;
    07 ) mon=Jul ;;
    08 ) mon=Aug ;;
    09 ) mon=Sep ;;
    10 ) mon=Oct ;;
    11 ) mon=Nov ;;
    12 ) mon=Dec ;;
esac
if [ ${ntype} = TE ]; then
cat > $cfile << EOF
DSET ^${ofile}
UNDEF 9.999E+20
OPTIONS big_endian template
ydef 361 linear -90.0 0.5
xdef 720 linear 0.0 0.5
zdef 3 levels 300 500 850
tdef 8 linear ${hh}Z${dd}${mon}${yyyy} 12hr
vars 6
ugrd 3 33,100,0 U component of wind [m/s]
vgrd 3 33,100,0 V component of wind [m/s]
t 3 33,100,0 Temperature [K]
q 3 33,100,0 Specific Humidity [kg/kg]
pres_meansealev 0 33,100,0  Pressure at mean sea level [hPa]
te 0 33,100,0 Total energy [J/kg]
endvars
EOF
elif [ ${ntype} = dTE ]; then
cat > $cfile << EOF
DSET ^${ofile}
UNDEF 9.999E+20
OPTIONS big_endian template
ydef 361 linear -90.0 0.5
xdef 720 linear 0.0 0.5
zdef 3 levels 300 500 850
tdef 8 linear ${hh}Z${dd}${mon}${yyyy} 12hr
vars 5
ugrd 3 33,100,0 U component of wind [m/s]
vgrd 3 33,100,0 V component of wind [m/s]
t 3 33,100,0 Temperature [K]
pres_meansealev 0 33,100,0  Pressure at mean sea level [hPa]
te 0 33,100,0 Total energy [J/kg]
endvars
EOF
fi
cat $cfile
cdo -f nc import_binary ${cfile} ${ncfile}
ncl -nQ d=0 mode=2 orig=\"${orig}\" yyyymmddhh=\"${idate}\" ../ncl/ensvsa-${ntype}.ncl
ncl -nQ d=6 mode=2 orig=\"${orig}\" yyyymmddhh=\"${idate}\" ../ncl/ensvsa-${ntype}.ncl
ls -ltr