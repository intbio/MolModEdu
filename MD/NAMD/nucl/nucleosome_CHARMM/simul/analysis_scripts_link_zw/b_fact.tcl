mol load psf ../analysis_data/only_nucl_init.psf
#mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol addfile ../analysis_data/md_nucl.dcd first 600 waitfor all

source libs/bfactor.tcl
set nframes [expr  [molinfo top get numframes] ]
#tweek to reduce load
#set nframes 100
#set sel [atomselect top "all and noh"]
bfactor "all and noh" 0 $nframes "../analysis_data/nucl_bfactors.pdb" 10
#the scaling is needed to get normal values, that can be written in pdb.
#you have to correct during visualization for that!!!
#hbonds -sel1 $chi -sel2 $chj -dist 3.0 -ang 30 -type unique -writefile yes -outfile ../analysis_data/dna_hbonds.dat

exit


