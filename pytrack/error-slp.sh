#!/bin/sh
init=${init:-2019101212}

bstfile=bst_hagibis.txt

for orig in ecmf egrr kwbc rjtd;do
awk -f error-slp.awk ${bstfile} $orig/track${init}_mod.txt > $orig/error-slp${init}_mod.txt
done
