insertion=$1

shift ||{
    exit
}


sed 's/^>/> /' "$@" |

awk '
BEGIN {while (getline < "'"$insertion"'") a[$1] = 1}
/^>/ {i = $2 in a} i
' |

sed 's/^> />/'
