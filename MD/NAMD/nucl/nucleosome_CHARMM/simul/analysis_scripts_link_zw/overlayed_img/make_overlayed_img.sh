#!/bin/bash

vmdpy -e over_front.tcl
vmdpy -e over_back.tcl
vmdpy -e over_right.tcl
vmdpy -e over_left.tcl
vmdpy -e over_bottom.tcl

composite -blend 70 ../../analysis_data/front.tga ../../analysis_data/front_overlayed.tga ../../analysis_data/front_overlayed_1.png
composite -blend 70 ../../analysis_data/back.tga ../../analysis_data/back_overlayed.tga ../../analysis_data/back_overlayed_1.png
composite -blend 70 ../../analysis_data/left.tga ../../analysis_data/left_overlayed.tga ../../analysis_data/left_overlayed_1.png
composite -blend 70 ../../analysis_data/right.tga ../../analysis_data/right_overlayed.tga ../../analysis_data/right_overlayed_1.png
composite -blend 70 ../../analysis_data/bottom.tga ../../analysis_data/bottom_overlayed.tga ../../analysis_data/bottom_overlayed_1.png

rm -rf dat
