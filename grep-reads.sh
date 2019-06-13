clusters=$1
reads=$2
selected=$3


awk '{ print $6 }' $clusters | uniq > cluster-reads
while read LINE;do
	grep -A 2 $LINE $reads >> $selected
done < cluster-reads