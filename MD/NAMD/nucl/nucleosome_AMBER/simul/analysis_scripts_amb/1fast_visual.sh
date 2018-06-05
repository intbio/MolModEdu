#!/bin/bash
cd fast_render
rm -f ../../analysis_data/movie.mov
vmdpy -e movie.tcl
./make_movie.sh
rm -rf dat

cd ../overlayed_img
./make_overlayed_img.sh
