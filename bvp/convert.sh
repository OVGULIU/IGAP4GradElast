#! /bin/sh

# converts () -> []
# sh ./convert.sh ndofxnddim npar ne ncomp eval_residual.c

ndofxnddim=$1
npar=$2
ne=$3
ncomp=$4
fname=$5

i=0
while [ "$i" -lt "$ndofxnddim" ]
do
    sed -i $(echo 's/u('$i')/(*(u+'$i'))/g') $fname
    i=$(($i + 1))
done

l=0
while [ "$l" -lt "$npar" ]
do
    sed -i $(echo 's/par('$l')/(*(par+'$l'))/g') $fname
    l=$(($l + 1))
done

l=0
while [ "$l" -lt "$ne" ]
do
    sed -i $(echo 's/e('$l')/(*(e+'$l'))/g') $fname
    l=$(($l + 1))
done

l=0
while [ "$l" -lt "$ncomp" ]
do
    sed -i $(echo 's/comp('$l')/(*(comp+'$l'))/g') $fname
    l=$(($l + 1))
done

sed -i $(echo 's/Power/pow/g') $fname
sed -i $(echo 's/Sqrt/sqrt/g') $fname
sed -i $(echo 's/Log/log/g') $fname

exit 0

