cluster=$1

shift || {
    cat <<EOF
Usage: $0 pairwise-alignments.maf sequences.fasta > out.fasta
Get the FASTA sequences whose names appear in bottom MAF lines.
EOF
    exit
}

sed 's/^>/> /' "$@" |

awk '
BEGIN {while (getline < "'"$cluster"'") a[$1] = 1}
/^>/ {i = $1 in a} i
' |

sed 's/^> />/'

