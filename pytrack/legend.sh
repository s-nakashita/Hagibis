#/bin/sh
if [ $# -lt 1 ]; then
  echo "$0 cptfile initfile"
  exit
fi
cptfile=${1}
initfile=${2}
awk '$1!~/[BFN]/{print $2}' ${cptfile} > cpt.txt
paste cpt.txt ${initfile} > cpt_init.txt
echo 'S 0c c 0.1i 255/255/255/0 1p,0/0/0 0.5c best track' > legend.txt
awk '{print "S 0c c 0.1i "$1" 1p,"$1" 0.5c "$2}' cpt_init.txt >> legend.txt
###rm -f cpt.txt cpt_init.txt
