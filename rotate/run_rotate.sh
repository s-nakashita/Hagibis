#!/bin/bash
model=${model:-40r1v2}
init=${init:-2018070112}
ntrunc=${ntrunc:-1023}
nlev=${nlev:-60}
nlon=${nlon:-300}
nlat=${nlat:-51}
dlat=${dlat:-0.2}

datadir=/Volumes/Pegasus/OpenIFS/runs/${model}/TL${ntrunc}L${nlev}/${init}
echo ${datadir}

if [ ! -d ${datadir}/rotate ]; then
  mkdir ${datadir}/rotate
fi

nfrpos=$(grep NFRPOS ${datadir}/fort.4 | awk -F = '{print $2}' | sed -e 's/,//')
#nstop=$(grep NSTOP ${datadir}/fort.4 | awk -F = '{print $2}' | sed -e 's/,//')
eid=$(basename ${datadir}/ICMSH*INIT | sed -e 's/[A-Z]//g')

trackfile=${datadir}/track.txt
nl=$(wc -l ${trackfile} | awk '{print $1}')
nstop=$(((nl - 1) * nfrpos))

ncprefix=${datadir}/nc/ICMGG${eid}+
ncsuffix="_isobaricInhPa.nc"

varname=q
outfile=${datadir}/rotate/${varname}.grd
echo ${outfile}
ncl -nQ nlon=${nlon} nlat=${nlat} dlat=${dlat} nfrpos=${nfrpos} nstop=${nstop} \
        varname=\"${varname}\" trackfile=\"${trackfile}\" ncprefix=\"${ncprefix}\" ncsuffix=\"${ncsuffix}\" \
        outfile=\"${outfile}\" rotate_scalar.ncl

ncprefix=${datadir}/nc/ICMSH${eid}+

varname=T
outfile=${datadir}/rotate/${varname}.grd
echo ${outfile}
ncl -nQ nlon=${nlon} nlat=${nlat} dlat=${dlat} nfrpos=${nfrpos} nstop=${nstop} \
        varname=\"${varname}\" trackfile=\"${trackfile}\" ncprefix=\"${ncprefix}\" ncsuffix=\"${ncsuffix}\" \
        outfile=\"${outfile}\" rotate_scalar.ncl

varname=W
outfile=${datadir}/rotate/${varname}.grd
echo ${outfile}
ncl -nQ nlon=${nlon} nlat=${nlat} dlat=${dlat} nfrpos=${nfrpos} nstop=${nstop} \
        varname=\"${varname}\" trackfile=\"${trackfile}\" ncprefix=\"${ncprefix}\" ncsuffix=\"${ncsuffix}\" \
        outfile=\"${outfile}\" rotate_scalar.ncl

varname=VO
outfile=${datadir}/rotate/${varname}.grd
echo ${outfile}
ncl -nQ nlon=${nlon} nlat=${nlat} dlat=${dlat} nfrpos=${nfrpos} nstop=${nstop} \
        varname=\"${varname}\" trackfile=\"${trackfile}\" ncprefix=\"${ncprefix}\" ncsuffix=\"${ncsuffix}\" \
        outfile=\"${outfile}\" rotate_scalar.ncl

varname=D
outfile=${datadir}/rotate/${varname}.grd
echo ${outfile}
ncl -nQ nlon=${nlon} nlat=${nlat} dlat=${dlat} nfrpos=${nfrpos} nstop=${nstop} \
        varname=\"${varname}\" trackfile=\"${trackfile}\" ncprefix=\"${ncprefix}\" ncsuffix=\"${ncsuffix}\" \
        outfile=\"${outfile}\" rotate_scalar.ncl

outfile=${datadir}/rotate/UV.grd
echo ${outfile}
ncl -nQ nlon=${nlon} nlat=${nlat} dlat=${dlat} nfrpos=${nfrpos} nstop=${nstop} \
        trackfile=\"${trackfile}\" ncprefix=\"${ncprefix}\" ncsuffix=\"${ncsuffix}\" \
        outfile=\"${outfile}\" rotate_uv.ncl

ncprefix=${datadir}/nc/ICMGG${eid}+
ncsuffix="_surface.nc"
varname=msl
outfile=${datadir}/rotate/slp.grd
echo ${outfile}
ncl -nQ nlon=${nlon} nlat=${nlat} dlat=${dlat} nfrpos=${nfrpos} nstop=${nstop} \
        varname=\"${varname}\" trackfile=\"${trackfile}\" ncprefix=\"${ncprefix}\" ncsuffix=\"${ncsuffix}\" \
        outfile=\"${outfile}\" rotate_scalar.ncl
