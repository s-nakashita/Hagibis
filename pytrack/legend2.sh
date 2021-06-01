#/bin/sh
cptfile=${1}
awk '$1!~/[BFN]/{print $2}' ${cptfile} > cpt_tigge.txt
paste cpt_tigge.txt init_c.txt > cpt_init.txt
#echo 'N 2' > legend2.txt # legend column set to 2
echo 'S 0c c 0.1i 255/255/255/0 1p,0/0/0 0.5c best track' >> legend2.txt
#cat cpt_init.txt >> legend2.txt
awk '{print "S 0c c 0.1i "$1" 1p,"$1" 0.5c "$2}' cpt_init.txt >> legend2.txt
#rm -f cpt.txt cpt_init.txt
