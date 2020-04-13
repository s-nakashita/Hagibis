#!/bin/bash
palette=blue,red,forestgreen,greenyellow,blue,red,forestgreen,greenyellow
cptfile=track_tigge.cpt
gmt makecpt -C$palette -T0/8/1 > ${cptfile}
#gmt makecpt -A10 -C$palette -T0/4/1 >> ${cptfile}
#for center in jma ncep ukmo; do
#    gmt makecpt -C${palette} -T${t0}/${t1}/${dt} >> ${cptfile}
#done

