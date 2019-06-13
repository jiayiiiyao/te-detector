chr=$1
start=$2
res=L1-rmsk-$1-$2

grep $chr L1-altered-rmsk | grep -C 10 $start > $res

