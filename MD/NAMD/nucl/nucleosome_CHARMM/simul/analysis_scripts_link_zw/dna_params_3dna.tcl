#we calculate DNA params
# you need to make sure that in first fram all pairs are identified.
# for 1kx5 we had to adjust parameter <hb_dist1>4.0</hb_dist1> to <hb_dist1>4.2</hb_dist1>
# in misc_3dna.par


exec mkdir -p ../analysis_data/dna_params

mol load psf ../analysis_data/only_nucl_init.psf

mol addfile ../analysis_data/only_nucl_init.pdb waitfor all
mol addfile ../analysis_data/md_nucl.dcd first 100 step 1 waitfor all

#

set dna [atomselect top "nucleic"]


set nframes [expr  [molinfo top get numframes] - 1 ]

$dna frame 0

$dna writepdb ../analysis_data/dna_params/cryst.pdb
exec ./run_x3dna.sh cryst cryst 2> ../analysis_data/dna_params/stderr.log


for { set i 1 } { $i<=$nframes } { incr i } {

$dna frame $i

$dna writepdb ../analysis_data/dna_params/step$i.pdb

exec ./run_x3dna.sh cryst step$i 2> ../analysis_data/dna_params/stderr.log



#This line is deprecated
#exec python2.7 myplot_dna_par.py ../analysis_data/dna_params/step_$i.dat 2> ../analysis_data/dna_params/stderr.log

}

exec /opt/local/bin/python2.7 dna_par_3dna_recombine.py

#Deprecated
# exec rm -f ../analysis_data/dna_param.mov
# exec ffmpeg -i ../analysis_data/dna_params/step_%d.png -s 1000x750 -q:v 0 -pix_fmt yuv420p  ../analysis_data/dna_param.mov




#exec rm -rf ../analysis_data/dna_params

exit


