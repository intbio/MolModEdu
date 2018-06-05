

source fast_render/input_param.tcl
source fast_render/add_text_layer.tcl
set scale 1.5

 mol load psf ../analysis_data/only_nucl_init.psf
 mol addfile ../analysis_data/only_nucl_init.pdb
 mol ssrecalc top
 mol addfile ../analysis_data/md_nucl.dcd waitfor all



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
mol color ColorID 17
mol selection { nucleic and noh and not name P O1P O2P O3' O5' C5' and resname GUA CYT}
mol material AOEdgy
mol addrep top
mol selupdate 5 top 0
mol colupdate 5 top 0

mol representation Licorice 0.300000 10.000000 10.000000
mol color ColorID 25
mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname THY ADE}
mol material AOEdgy
mol addrep top
mol selupdate 6 top 0
mol colupdate 6 top 0


#color scheme  for histones and DNA
color Segname CHA orange2
color Segname CHB yellow3
color Segname CHC mauve
color Segname CHD cyan3
color Segname CHE orange2
color Segname CHF yellow3
color Segname CHG mauve
color Segname CHH cyan3
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
 draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "MD simulations of NCP with truncated tails (100 ns)" size 1.5 thickness 3
 incr txtlncount 2

 add_text_layer HISTONES

# #Turn off all
set num_of_reps [molinfo 0 get numreps]
for {set i 0} { $i < $num_of_reps } { incr i 1 } {
mol showrep 0 $i off }

set matnum 22

set text(0) "Histones H3"
set text(1) "Histones H4"
set text(2) "Histones H2A"
set text(3) "Histones H2B"
set color_text(0) 31
set color_text(1) 18
set color_text(2) 13
set color_text(3) 22

# #First four pair of histones
for {set r 0} {$r < 4} {incr r 1} {
	#text
	draw color $color_text($r)
	draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "$text($r)" size 1.5 thickness 3
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
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "DNA" size 1.5 thickness 3
incr txtlncount 1

draw color 25
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "AT pairs" size 1.5 thickness 3
incr txtlncount 1

draw color 17
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "GC pairs" size 1.5 thickness 3
incr txtlncount 1

draw color 1
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Sodium ions" size 1.5 thickness 3
incr txtlncount 1

draw color 7
draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Chloride ions" size 1.5 thickness 3
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
	
display update

