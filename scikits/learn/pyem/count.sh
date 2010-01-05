#! /bin/sh
# Last Change: Mon Sep 11 08:00 PM 2006 J

n=0
np=0
nc=0

files=`find . -name '*.py'`

for i in $files; do
    tp=`wc "$i" | tr -s " " | cut -f 2 -d" "`
    let np="$np + $tp"
done

files=`find . -name '*.[ch]'`

for i in $files; do
    tp=`wc "$i" | tr -s " " | cut -f 2 -d" "`
    let nc="$nc + $tp"
done

echo $nc
echo $np
