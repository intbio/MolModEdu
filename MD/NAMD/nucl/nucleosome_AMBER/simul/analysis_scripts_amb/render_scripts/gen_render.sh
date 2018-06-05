#!/bin/bash
rm run_sw_tach.sh
cd dat;
f=`ls`
cd ..
for i in $f
do
echo "./tachyon -trans_vmd -shade_blinn -fullshade -aasamples 10 -rescale_lights 0.5 -add_skylight 0.8 dat/$i -format TARGA -o tga/$i.tga" >> run_sw_tach.sh
done
