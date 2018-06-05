#!/bin/bash
f=0
for i in `seq 1000 0`;
do
if [ ! -f ${i}.dat.tga ]; then
if [ $f == 0 ]; then
k=`expr $i + 1`
cp ${k}.dat.tga ${i}.dat.tga
echo "do"
f=1
fi
else
f=0
fi
done
