

source fast_render/input_param.tcl
source fast_render/add_text_layer.tcl


####
#Display the reference crystal structure
# mol load psf ../analysis_data/only_nucl_init.psf
mol load pdb ../analysis_data/only_nucl_average.pdb

mol modstyle 0 0 NewCartoon 0.840000 20.000000 2.630000 0

###

 # mol load psf ../analysis_data/only_nucl_init.psf
 mol load pdb ../analysis_data/only_nucl_init.pdb
 mol ssrecalc top
 # mol addfile ../analysis_data/md_nucl.dcd step 10 waitfor all 

set molid [ molinfo top ]
 
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

# H4 H4 Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHB CHF}
mol material AOShiny
mol addrep top
mol selupdate 1 top 0
mol colupdate 1 top 0

# H2a H2a Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHC CHG}
mol material AOShiny
mol addrep top
mol selupdate 2 top 0
mol colupdate 2 top 0

# H2b H2b Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Segname
mol selection {segname CHD CHH}
mol material AOShiny
mol addrep top
mol selupdate 3 top 0
mol colupdate 3 top 0

# Dna sugar-phosphate backbone representation
mol representation NewCartoon 1.090000 10.000000 2.070000 1
mol color ColorID 6
mol selection {nucleic and backbone}
mol material AOEdgy
mol addrep top
mol selupdate 4 top 0
mol colupdate 4 top 0

# Dna nucleobases representation
mol representation Licorice 0.300000 10.000000 10.000000
mol color ColorID 6
mol selection { nucleic and noh and not name P O1P O2P O3' O5' C5' and resname GUA CYT}
mol material AOEdgy
mol addrep top
mol selupdate 5 top 0
mol colupdate 5 top 0

mol representation Licorice 0.300000 10.000000 10.000000
mol color ColorID 3
mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname THY ADE}
mol material AOEdgy
mol addrep top
mol selupdate 6 top 0
mol colupdate 6 top 0


#color scheme  for histones and DNA
#old scheme
# color Segname CHA orange2
# color Segname CHB yellow3
# color Segname CHC mauve
# color Segname CHD cyan3
# color Segname CHE orange2
# color Segname CHF yellow3
# color Segname CHG mauve
# color Segname CHH cyan3

color Segname CHA blue3
color Segname CHB green
color Segname CHC yellow2
color Segname CHD red3
color Segname CHE blue3
color Segname CHF green
color Segname CHG yellow2
color Segname CHH red3

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
 draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "MD simulations of NCP without tails (1000 ns)" size 1.5 thickness 3
 incr txtlncount 1
 draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Crystal structure vs average simulation structure (in cyan)" size 1.5 thickness 3
 incr txtlncount 2

 add_text_layer HISTONES

# #Turn off all
set num_of_reps [molinfo $molid get numreps]
#for {set i 0} { $i < $num_of_reps } { incr i 1 } {
#mol showrep $molid $i off }

set matnum 22

set text(0) "Histones H3"
set text(1) "Histones H4"
set text(2) "Histones H2A"
set text(3) "Histones H2B"
set color_text(0) blue3
set color_text(1) green
set color_text(2) yellow2
set color_text(3) red3

# #First four pair of histones
for {set r 0} {$r < 4} {incr r 1} {
	#text
	draw color $color_text($r)
	draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "$text($r)" size 1.5 thickness 3
	incr txtlncount 1
	#molecule
	mol showrep $molid $r on
	#copy a material for them
	material add copy AOShiny
	# mol modmaterial $r top Material$matnum
	mol modmaterial $r top AOShiny
	
	
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
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "DNA" size 1.5 thickness 3
incr txtlncount 1

draw color 3
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "AT pairs" size 1.5 thickness 3
incr txtlncount 1

draw color 6
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "GC pairs" size 1.5 thickness 3
incr txtlncount 1

#ion lables
# draw color 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Sodium ions" size 1.5 thickness 3
# incr txtlncount 1

# draw color 7
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Chloride ions" size 1.5 thickness 3
# incr txtlncount 1

#molecule
mol showrep $molid 4 on
mol showrep $molid 5 on
mol showrep $molid 6 on
mol showrep $molid 7 on
mol showrep $molid 8 on
mol modmaterial 4 $molid AOEdgy
mol modmaterial 5 $molid AOEdgy
mol modmaterial 6 $molid AOEdgy
mol modmaterial 7 $molid AOEdgy
mol modmaterial 8 $molid AOEdgy

	display update
	#material change opacity AOEdgy 0.8
	
display update

add_text_layer TIME
draw color 0
set time [format "Front"]
draw text " $txtx [expr $txty-(27*$txtstep)] 0 " $time size 1.5 thickness 3

render snapshot ../analysis_data/cryst_vs_aver_front.tga
rotate y by 90

mol delete top
add_text_layer TIME
draw color 0
set time [format "Right"]
draw text " $txtx [expr $txty-(27*$txtstep)] 0 " $time size 1.5 thickness 3

render snapshot ../analysis_data/cryst_vs_aver_right.tga
rotate y by 90

mol delete top
add_text_layer TIME
draw color 0
set time [format "Back"]
draw text " $txtx [expr $txty-(27*$txtstep)] 0 " $time size 1.5 thickness 3

render snapshot ../analysis_data/cryst_vs_aver_back.tga
rotate y by 90

mol delete top
add_text_layer TIME
draw color 0
set time [format "Left"]
draw text " $txtx [expr $txty-(27*$txtstep)] 0 " $time size 1.5 thickness 3

render snapshot ../analysis_data/cryst_vs_aver_left.tga

mol top 1
exit

