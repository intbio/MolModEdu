#!/bin/bash
ffmpeg -i dat/%d.dat.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  ../../analysis_data/movie.mov
ffmpeg -i dat/%d.dat.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  ../../analysis_data/movie.wmv

composite -blend 70 ../../analysis_data/initial.tga ../../analysis_data/overlayed.tga ../../analysis_data/overlayed_1.png
