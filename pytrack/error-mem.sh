#!/bin/sh
init=${init:-2019100912}

bstfile=bst_hagibis.txt
orig=ukmo
case $orig in
    ecmwf ) M=50 ;;
    jma   ) M=26 ;;
    ncep  ) M=20 ;;
    ukmo  ) M=17 ;;
esac

for m in `seq $M`;do
if [ $m -lt 10 ]; then
    m=0$m
fi
awk -f error.awk ${bstfile} $orig/track${init}_${m}.txt > $orig/error${init}_${m}.txt
done
