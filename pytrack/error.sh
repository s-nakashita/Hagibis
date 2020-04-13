#!/bin/sh
init=${init:-2019101312}

bstfile=bst_hagibis.txt

for orig in ecmf egrr kwbc rjtd;do
awk -f error.awk ${bstfile} $orig/track${init}.txt > $orig/error${init}.txt
done
