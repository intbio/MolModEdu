# we calculate ion volume map and visualize it
exec mkdir -p dat
package require colorscalebar
axes location off
source fast_render/input_param.tcl
source fast_render/add_text_layer.tcl

display resize 2000 1500

color scale method BWR
color scale midpoint 0.5

set minsc 0.3
set maxsc 2.0

mol load psf ../input/1kx5_ready.psf
mol addfile ../analysis_data/system_aligned.pdb waitfor all
mol ssrecalc top
mol addfile ../analysis_data/md_nucl_solv.dcd first 100 waitfor all

volmap density [atomselect top "resname SOD"] -res 0.5 -minmax {{-80.0 -80.0 -50.0} {80.0 80.0 50.0}} -o ../analysis_data/ion_sod_volmap.dx -mol top -allframes -combine avg -checkpoint 10000
volmap density [atomselect top "resname CLA"] -res 0.5 -minmax {{-80.0 -80.0 -50.0} {80.0 80.0 50.0}} -o ../analysis_data/ion_cla_volmap.dx -mol top -allframes -combine avg -checkpoint 10000

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
# mol smoothrep top 0 5

# H4 H4 Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHB CHF}
mol material AOShiny
mol addrep top
mol selupdate 1 top 0
mol colupdate 1 top 0
# mol smoothrep top 1 5

# H2a H2a Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHC CHG}
mol material AOShiny
mol addrep top
mol selupdate 2 top 0
mol colupdate 2 top 0
# mol smoothrep top 2 5

# H2b H2b Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHD CHH}
mol material AOShiny
mol addrep top
mol selupdate 3 top 0
mol colupdate 3 top 0
# mol smoothrep top 3 5

# Dna sugar-phosphate backbone representation
mol representation NewCartoon 1.090000 10.000000 2.070000 1
mol color ColorID 6
mol selection {nucleic and backbone}
mol material AOEdgy
mol addrep top
mol selupdate 4 top 0
mol colupdate 4 top 0
# mol smoothrep top 4 5/

# Dna nucleobases representation
mol representation Licorice 0.300000 10.000000 10.000000
mol color ColorID 6
mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname GUA CYT}
mol material AOEdgy
mol addrep top
mol selupdate 5 top 0
mol colupdate 5 top 0
# mol smoothrep top 5 5

mol representation Licorice 0.300000 10.000000 10.000000
mol color ColorID 6
mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname THY ADE}
mol material AOEdgy
mol addrep top
mol selupdate 6 top 0
mol colupdate 6 top 0
# mol smoothrep top 6 5


#{{length 0.5} {width 0.05} {auto_scale 1} {fixed 1} {min 0} {max 100} {label_num 5} {text 16} {fp_format 0} {x_pos -1.0} {y_pos -1.0} {replacebar 1} {molid top} {repid 0} {showlegend 0} {legend "Color Scale:"}} {
# ::ColorScaleBar::color_scale_bar 1.5 0.1 1 1 0 100 5 16 0 -1.0 -1.0 1 top 0 1 "RMSF, A" 

# mol fix 0
# mol free 1

# translate by -0.6 0.5 0

# mol free 0
# mol fix 1

#color scheme  for histones and DNA
color Segname CHA 6
 #blue3
color Segname CHB 6 
#green
color Segname CHC 6 
#yellow2
color Segname CHD 6 
#red3
color Segname CHE 6 
#blue3
color Segname CHF 6 
#reen
color Segname CHG 6 
#ellow2
color Segname CHH 6 
#red3


color Highlight Nonback 6
color Highlight Nucback 2
color Display Background white
color scale method RWB
color change rgb  0 0.1 0.2 0.7 ;# blue
color change rgb  1 0.7 0.2 0.1 ;# red
color change rgb  3 0.7 0.4 0.0 ;# orange
color change rgb  4 0.8 0.7 0.1 ;# yellow
color change rgb  7 0.1 0.7 0.2 ;# green
color change rgb 10 0.1 0.7 0.8 ;# cyan
color change rgb 11 0.6 0.1 0.6 ;# purple
color change rgb 23 0.550000011920929 0.7099999785423279 0.9800000190734863



# ================= movie starts ==================== 

# sets a variable for knowing were real molecules ends if there are more than one
 set molnum [expr [molinfo num]]

# #--------set start txt position in px-----
 set txtx [expr -($dispw/2) + 20 ]
 set txty [expr ($disph/2) - 20 ]

 # set step of txt lines
 set txtstep  [expr $disph/30 ]

# # set counter for lines
 set txtlncount 0

# #prepare scene
 scale by $scale
 translate by $transx $transy $transz
 axes location off
 display update ui

# #add some text on new layer
 add_text_layer HEADNAME
 draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Ion distribution around nucleosome" size 1.5 thickness 6
 incr txtlncount 1
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Isosurfaces for average number density of ions" size 1.5 thickness 6
 incr txtlncount 2


 add_text_layer HISTONES

# #Turn off all
set num_of_reps [molinfo 0 get numreps]
for {set i 0} { $i < $num_of_reps } { incr i 1 } {
mol showrep 0 $i off }

#set framecounter
set fc 0

# #render the first frame
# #render Tachyon dat/$fc.dat
# incr fc 1

# # h3_1-2 = 1------------------------------------
# # we'll copy materials and should know wich is the last
# #so $matstartnum is a numbet of starting material name

set matnum 22

set text(0) "Histones H3"
set text(1) "Histones H4"
set text(2) "Histones H2A"
set text(3) "Histones H2B"
set color_text(0) 6 
#blue3
set color_text(1) 6 
#green
set color_text(2) 6 
#yellow2
set color_text(3) 6 
#red3

# #First four pair of histones
for {set r 0} {$r < 4} {incr r 1} {
	#text
	draw color $color_text($r)
	draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "$text($r)" size 1.5 thickness 6
	incr txtlncount 1
	#molecule
	mol showrep 0 $r on
	#copy a material for them
	material add copy AOShiny
	mol modmaterial $r 0 Material$matnum
	
	
# 	#make a smooth transparent appearance for 30 frames with rotation
	
		material change opacity Material$matnum 1
		#render Tachyon dat/$fc.dat
		
	

# 	

	incr matnum 1
}

# # Adding DNA
# #text
add_text_layer DNA
draw color 6
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "DNA" size 1.5 thickness 6
incr txtlncount 1

# draw color 3
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "AT pairs" size 1.5 thickness 6
# incr txtlncount 1

# draw color 6
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "GC pairs" size 1.5 thickness 6
# incr txtlncount 1

draw color 27
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Sodium ions" size 1.5 thickness 6
incr txtlncount 1

draw color 21
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Chloride ions" size 1.5 thickness 6
incr txtlncount 1


#molecule
mol showrep 0 4 on
mol showrep 0 5 on
mol showrep 0 6 on
mol showrep 0 7 on
mol showrep 0 8 on
mol modmaterial 4 0 AOEdgy
mol modmaterial 5 0 AOEdgy
mol modmaterial 6 0 AOEdgy
mol modmaterial 7 0 AOEdgy
mol modmaterial 8 0 AOEdgy

	display update
	#material change opacity AOEdgy 0.8
	
mol top 0


mol fix 1

display update


mol material AOEdgy
mol addrep 0
mol modcolor 7 0 ColorID 27
mol modstyle 7 0 Isosurface 0.001120 0 0 0 1 1

mol addrep 0
mol modcolor 8 0 ColorID 21
mol modstyle 8 0 Isosurface 0.000800 1 0 0 1 1

# mol color ColorID 1
# mol representation Isosurface 0.001432 0 2 0 1 1
# mol selection all
# mol material Opaque
# mol addrep 0

# mol color ColorID 7
# mol representation Isosurface 0.000651 1 2 0 1 1
# mol selection all
# mol material Opaque
# mol addrep 0

render snapshot ../analysis_data/nucl_iso_ions.tga

for {set i 0} {$i <= 360} {incr i 1} {

rotate y by 1
display update
render snapshot dat/$i.tga

}

exec rm -f ../analysis_data/nucl_iso_ions.mov
exec ffmpeg -i dat/%d.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  ../analysis_data/nucl_iso_ions.mov


# exec rm -rf dat




#exit


