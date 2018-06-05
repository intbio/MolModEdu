#!/bin/bash


# We need here to characterize ion distribution and behavior
# RDF from DNA, distance from surface
# volmap
# condensation events !!! - still need to do, cf. Papoyan in JACS
# Manning radius ? Papoyan said is around 1 nm for DNA(?)
# major vs minor groove - still need to do - compare with literature and only DNA sumulations
# Sodium is in minor groove - in accordance with previous studies.
#Let's start with rdf and also check for equlibration

vmdtpy -e ion_rdf.tcl # will calculate plot dna-ion rdfs
vmdtpy -e ion_atm.tcl # will calculate number of ions as a function of distance from nucleosome 
vmdtpy -e ion_volmap.tcl # will calculate ion volmaps and make a movie


#TODO: characterize condenstaion, probably by residence times of distance between atoms (see Papoyan work)
# TODO: trace profiles of ions along the DNA (?)



