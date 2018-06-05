#!/bin/bash
rm run_sw_vmd.sh
for i in `seq 0 19`;
do
first=` expr $i \* 50 `
last=` expr $first + 49 `
echo "./run_vmd_batched.sh $first $last > log$i.out &" >> run_sw_vmd.sh
done
