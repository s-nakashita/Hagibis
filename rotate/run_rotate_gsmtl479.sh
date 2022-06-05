#!/bin/sh
init=2019100900
for exp in cntl; do
    echo $exp
    mkdir -p outline_${exp}/${init}
    python rotate_scalar_gsmtl479.py $exp $init
    python rotate_uv_gsmtl479.py $exp $init
    python create_netcdf_gsmtl479.py $exp $init
    rm np_fcst*_${init}.nc
done