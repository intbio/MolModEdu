#!/bin/bash


# We need here to characterize water in nucleosome
# we start with a volmap


vmdpy -e water_volmap.tcl # will calculate ion volmaps and make a movie
vmdtpy -e water_volmap_msms.tcl # will calculate ion volmaps and make a movie

#TODO: characterize residence times of water (probably this will be in 6contacts.sh)
# diffusion, water at interfaces, compare with crystallographic water



