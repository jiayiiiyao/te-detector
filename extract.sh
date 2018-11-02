insertion=$1

shift ||{
    exit
}


sed 's/^>/> /' "$@" |

awk '
BEGIN {while (getline < "'"$insertion"'") a[$4] = 1}
/^>/ {i = $4 in a} i
' |

sed 's/^> />/'
