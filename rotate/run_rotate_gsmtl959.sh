#!/bin/sh
init="2019100912"
for exp in clim est mgdsst; do
    echo $exp $init
    mkdir -p outline_${exp}/${init}
    python rotate_scalar_gsmtl959.py $exp $init
    python rotate_uv_gsmtl959.py $exp $init
    python create_netcdf_gsmtl959.py $exp $init
    rm np_fcst*_asia_${init}.nc
done