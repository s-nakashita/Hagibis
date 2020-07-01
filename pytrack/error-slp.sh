#!/bin/sh
init=${init:-2019101212}

bstfile=bst_hagibis.txt

for orig in ecmf egrr kwbc rjtd;do
awk -f error-slp.awk ${bstfile} $orig/track${init}.txt > $orig/error-slp${init}.txt
done
