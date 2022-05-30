#!/bin/sh
init=${init:-2019101312}

bstfile=bst_hagibis.txt

for orig in ecmf egrr kwbc rjtd;do
awk -f error.awk ${bstfile} $orig/track${init}_mod.txt > $orig/error${init}_mod.txt
done
