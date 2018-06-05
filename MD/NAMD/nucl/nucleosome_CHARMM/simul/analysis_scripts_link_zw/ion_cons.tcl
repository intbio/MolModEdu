# we calculate ion volume map and visualize it
exec mkdir -p dat
package require colorscalebar
axes location off
source fast_render/input_param.tcl
source fast_render/add_text_layer.tcl



color scale method BWR
color scale midpoint 0.33

set minsc 0.0
set maxsc 0.05

mol load psf ../input/1kx5_ready.psf
mol addfile ../analysis_data/system_aligned.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl_solv.dcd first 500 waitfor all

# volmap density [atomselect top "resname TIP3 and name OH2"] -res 0.5 -minmax {{-52.0 -53.0 -35.0} {52.0 53.0 35.0}} -o ../analysis_data/water_volmap.dx -mol top -allframes -combine avg -checkpoint 10000
#volmap density [atomselect top "resname CLA"] -res 0.5 -minmax {{-80.0 -80.0 -50.0} {80.0 80.0 50.0}} -o ../analysis_data/ion_cla_volmap.dx -mol top -allframes -combine avg -checkpoint 10000

set scale 3.0
set nf [molinfo top get numframes]	 

display rendermode GLSL
mol delrep 0 top

#  H3 H3 Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHA CHE}
mol material AOShiny
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0
mol smoothrep top 0 5
mol scaleminmax top  0 $minsc $maxsc

# H4 H4 Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHB CHF}
mol material AOShiny
mol addrep top
mol selupdate 1 top 0
mol colupdate 1 top 0
mol smoothrep top 1 5
mol scaleminmax top  1 $minsc $maxsc

# H2a H2a Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHC CHG}
mol material AOShiny
mol addrep top
mol selupdate 2 top 0
mol colupdate 2 top 0
mol smoothrep top 2 5
mol scaleminmax top  2 $minsc $maxsc

# H2b H2b Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHD CHH}
mol material AOShiny
mol addrep top
mol selupdate 3 top 0
mol colupdate 3 top 0
mol smoothrep top 3 5
mol scaleminmax top  3 $minsc $maxsc

# Dna sugar-phosphate backbone representation
mol representation NewCartoon 1.090000 10.000000 2.070000 1
mol color ColorID 6
mol selection {nucleic and backbone}
mol material AOEdgy
mol addrep top
mol selupdate 4 top 0
mol colupdate 4 top 0
mol smoothrep top 4 5
mol scaleminmax top  4 $minsc $maxsc

# Dna nucleobases representation
# mol representation Licorice 0.300000 10.000000 10.000000
# mol color ColorID 17
# mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname GUA CYT}
# mol material AOEdgy
# mol addrep top
# mol selupdate 5 top 0
# mol colupdate 5 top 0
# mol smoothrep top 5 5
# mol scaleminmax top  5 $minsc $maxsc

# mol representation Licorice 0.300000 10.000000 10.000000
# mol color ColorID 25
# mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname THY ADE}
# mol material AOEdgy
# mol addrep top
# mol selupdate 6 top 0
# mol colupdate 6 top 0
# mol smoothrep top 6 5
# mol scaleminmax top 6 $minsc $maxsc




# sets a variable for knowing were real molecules ends if there are more than one
 set molnum [expr [molinfo num]]

# #--------set start txt position in px-----

set fc 0
for {set i 0} {$i <= 99} {incr i 1} {

set pos [expr $i./100]
mol modstyle 5 0 VolumeSlice $pos 0.000000 2.000000 2.000000
#rotate y by 1
display update
render snapshot dat/$fc.tga

incr fc
}

for {set i 99} {$i > 0} {incr i -1} {

set pos [expr $i./100]
mol modstyle 5 0 VolumeSlice $pos 0.000000 2.000000 2.000000
#rotate y by 1
display update
render snapshot dat/$fc.tga

incr fc
}

exec rm -f ../analysis_data/nucl_water_slice.mov
exec ffmpeg -i dat/%d.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  ../analysis_data/nucl_water_slice.mov


exec rm -rf dat




#exit


