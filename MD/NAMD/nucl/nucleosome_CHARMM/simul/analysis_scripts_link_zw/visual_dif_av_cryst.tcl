exec mkdir -p dat

package require colorscalebar

axes location off
source fast_render/input_param.tcl
#source aa_hydro_scale.tcl
#source colorbar.tcl
source fast_render/add_text_layer.tcl



set scale 1.7

color scale method BWR
color scale midpoint 0.5
color scale offset -0.05

set minsc 0.0
set maxsc 5.0



#loading molecule
mol load psf ../analysis_data/only_nucl_init.psf
mol addfile ../analysis_data/dif_av_cryst.pdb 
display rendermode GLSL
# mol load psf 1kx5nt_ready.psf
# mol addfile md_nucl_solv.dcd step 1 first 0 last 1000  waitfor all



#mol new 1KX5.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top

#  H3 H3 Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Beta
mol selection {segname CHA CHE}
mol material AOShiny
mol addrep top
mol selupdate 0 top 0
mol colupdate 0 top 0
mol smoothrep top 0 5
mol scaleminmax top  0 $minsc $maxsc

# H4 H4 Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Beta
mol selection {segname CHB CHF}
mol material AOShiny
mol addrep top
mol selupdate 1 top 0
mol colupdate 1 top 0
mol smoothrep top 1 5
mol scaleminmax top  1 $minsc $maxsc

# H2a H2a Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Beta
mol selection {segname CHC CHG}
mol material AOShiny
mol addrep top
mol selupdate 2 top 0
mol colupdate 2 top 0
mol smoothrep top 2 5
mol scaleminmax top  2 $minsc $maxsc

# H2b H2b Dimer representation
mol representation NewCartoon 0.840000 20.000000 2.630000 0
mol color Beta
mol selection {segname CHD CHH}
mol material AOShiny
mol addrep top
mol selupdate 3 top 0
mol colupdate 3 top 0
mol smoothrep top 3 5
mol scaleminmax top  3 $minsc $maxsc

# Dna sugar-phosphate backbone representation
mol representation NewCartoon 1.090000 10.000000 2.070000 1
mol color Beta
mol selection {nucleic and backbone}
mol material AOEdgy
mol addrep top
mol selupdate 4 top 0
mol colupdate 4 top 0
mol smoothrep top 4 5
mol scaleminmax top  4 $minsc $maxsc

# Dna nucleobases representation
mol representation Licorice 0.300000 10.000000 10.000000
mol color Beta
mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname GUA CYT}
mol material AOEdgy
mol addrep top
mol selupdate 5 top 0
mol colupdate 5 top 0
mol smoothrep top 5 5
mol scaleminmax top  5 $minsc $maxsc

mol representation Licorice 0.300000 10.000000 10.000000
mol color Beta
mol selection {nucleic and noh and not name P O1P O2P O3' O5' C5' and resname THY ADE}
mol material AOEdgy
mol addrep top
mol selupdate 6 top 0
mol colupdate 6 top 0
mol smoothrep top 6 5
mol scaleminmax top 6 $minsc $maxsc


#{{length 0.5} {width 0.05} {auto_scale 1} {fixed 1} {min 0} {max 100} {label_num 5} {text 16} {fp_format 0} {x_pos -1.0} {y_pos -1.0} {replacebar 1} {molid top} {repid 0} {showlegend 0} {legend "Color Scale:"}} {
::ColorScaleBar::color_scale_bar 1.5 0.1 1 1 0 100 5 16 0 -1.0 -1.0 1 top 0 1 "Deviation, A" 

mol fix 0
mol free 1

translate by -0.6 0.5 0

mol free 0
mol fix 1

# mol color ColorID 1
# mol selection "resname SOD and within 10 of (protein or nucleic)"
# mol representation VDW
# mol addrep top
# mol selupdate 7 top on
# mol colupdate 7 top 0
# mol smoothrep top 7 5

# #Put the "charmm r_min/2" / 2^(1/6)
# set atsod [atomselect top "resname SOD"]
# $atsod set radius 1.2568

# mol color ColorID 7
# mol selection "resname CLA and within 10 of (protein or nucleic)"
# mol representation VDW
# mol addrep top
# mol selupdate 8 top on
# mol colupdate 8 top 0
# mol smoothrep top 8 5
# set atcla [atomselect top "resname CLA"]
# $atcla set radius 2.0223

# histones surface Resname - modificated for hydrophobycity with aa_hydro_scale.tcl
# mol representation MSMS 1.400000
# mol color ResName
# mol selection {protein}
# mol material AOEdgy
# mol addrep top
# mol selupdate 6 top 0
# mol colupdate 6 top 0
# mol scaleminmax top 0 0.000000 0.000000

# # histones surface colored as histones
# mol representation MSMS 1.400000 0.000000
# mol color Segname
# mol selection {protein}
# mol material AOShiny
# mol addrep top
# mol selupdate 7 top 0
# mol colupdate 7 top 0
# mol scaleminmax top 0 0.000000 0.000000

# this is a representation for white glow outside of surface
 
# mol representation VDW 1.300000 20.000000
# mol color ColorID 8
# mol selection {protein }
# mol material Transparent
# mol addrep top
# mol selupdate 8 top 0
# mol colupdate 8 top 0
# mol scaleminmax top 1 0.000000 0.000000
# material change ambient Transparent 1.000000
# material change diffuse Transparent 0.000000
# material change specular Transparent 0.000000
# material change shininess Transparent 0.000000
# material change opacity Transparent 0.150000
# material change outline Transparent 0.000000
# material change outlinewidth Transparent 0.000000
# material change transmode Transparent 0.000000

# # positive charges
# mol representation VDW 1.100000 20.000000
# mol color ColorID 0
# mol selection {resname LYS  ARG and name  NH2 NZ}
# mol material AOChalky
# mol addrep top
# mol selupdate 9 top 0
# mol colupdate 9 top 0
# mol scaleminmax top 2 0.000000 0.000000

# # negative charges
# mol representation VDW 1.100000 20.000000
# mol color ColorID 30
# mol selection {resname GLU ASP and name  OD1  OE1}
# mol material AOChalky
# mol addrep top
# mol selupdate 10 top 0
# mol colupdate 10 top 0
# mol scaleminmax top 3 0.000000 0.000000

# #DNA surface
# mol representation MSMS 1.400000 0.000000
# mol color Backbone
# mol selection {nucleic }
# mol material AOChalky
# mol addrep top
# mol selupdate 11 top 0
# mol colupdate 11 top 0
# mol scaleminmax top 4 0.000000 0.000000

# #DNA charges
# mol representation VDW 1.100000 12.000000
# mol color Name
# mol selection {name OP1 }
# mol material AOChalky
# mol addrep top
# mol selupdate 12 top 0
# mol colupdate 12 top 0
# mol scaleminmax top 5 0.000000 0.000000

#mol rename top 1KX5.pdb

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
 draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Difference between crystal and average simulations structures" size 1.5 thickness 3
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

# draw color 25
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "AT pairs" size 1.5 thickness 3
# incr txtlncount 1

# draw color 17
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "GC pairs" size 1.5 thickness 3
# incr txtlncount 1

# draw color 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Sodium ions" size 1.5 thickness 3
# incr txtlncount 1

# draw color 7
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Chloride ions" size 1.5 thickness 3
# incr txtlncount 1

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
	material change opacity AOEdgy 0.8
	
mol top 0
#set nframes [expr  [molinfo top get numframes] - 1 ]
#add_text_layer TIME
# for {set i $first} {$i <= $last} {incr i 1} {
#animate goto $i

#mol delete top
#add_text_layer TIME
#draw color 0
#set time [format "Time: %5.1f ns" [expr $i * 1.0]]
#draw text " $txtx [expr $txty-(27*$txtstep)] 0 " $time size 1.5 thickness 3



mol fix 1

display update
#if { $i == 500 } {
#mol showrep 0 7 off
#mol showrep 0 8 off
#}


render snapshot ../analysis_data/dif_av_cryst.tga

for {set i 0} {$i <= 360} {incr i 1} {

rotate y by 1
display update
render snapshot dat/$i.tga

}

exec rm -f ../analysis_data/dif_av_cryst.mov
exec ffmpeg -i dat/%d.tga -s 1000x750 -q:v 0 -pix_fmt yuv420p  ../analysis_data/dif_av_cryst.mov


exec rm -rf dat

#render Tachyon dat/$fc.dat
#draws a colorbar with some text for hydrophob.surface 
# add_text_layer COLORBAR

#adds a colorbar for hydrophobycity
# colorbar 200 20 $txtx -120 $r0 $r1 $g0 $g1 $b0 $b1 11 -2.5 3.0

# incr txtlncount 2
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Amino acids" size 1.0 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Hydrophobicity scale*" size 1.0 thickness 1
# draw text " $txtx [expr - ($disph / 2) + 20 ] 0 " "*) Shaytan et al. DOI: 10.1021/bm8015169  Table 2d" size 1.0 thickness 1

# turn on a hydrophob. surface
# mol showrep 0 6 on
# material add copy AOShiny
# mol modmaterial 6 0 Material$matnum
# for {set r 1} {$r < 31} {incr r 1} {
# 	display update
# 	rotate y by -1
# 	material change opacity Material$matnum [expr {double($r)/30}]
# 	#make cartoon histones thiner
# 	mol modstyle 0 0 NewCartoon [expr {0.84 - double($r)/30*0.74}] 20.000000 2.630000 0
# 	mol modstyle 1 0 NewCartoon [expr {0.84 - double($r)/30*0.74}] 20.000000 2.630000 0
# 	mol modstyle 2 0 NewCartoon [expr {0.84 - double($r)/30*0.74}] 20.000000 2.630000 0
# 	mol modstyle 3 0 NewCartoon [expr {0.84 - double($r)/30*0.74}] 20.000000 2.630000 0
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }
# incr matnum 1

#turn off cartoon histones and histones  text
# for {set i 0} { $i < 4} { incr i 1} {
# 	mol showrep 0 $i off
# }
# mol off 2

# a little bit big construction just for moving 3 letters =)
# mol top 3
# for {set i 0} {$i < 40} {incr i 1} {
# 	display update
# 	draw delete all
# 	draw color 6
# 	draw text " $txtx [expr $txty-(5*$txtstep)+($i*$disph/300)] 0 " "DNA" size 1.5 thickness 3
# 	rotate y by -1
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }
#add additional text
# mol top 4
# draw color 16
# set txtlncount 3
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Protein representation:" size 1.2 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "surface, probe radius 1.4 A" size 1.2 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "color - hydrophobicity scale" size 1.2 thickness 1

# for {set i 0} {$i < 290} {incr i 1} {
# 	display update
# 	rotate y by -1
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }

# fade off the DNA
# for {set r 1} {$r < 31} {incr r 1} {
# 	display update
# 	material change opacity AOEdgy [expr {0.8-double($r)/30*0.8}]
# 	rotate y by -1
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }
# turn off the DNA
# mol showrep 0 5 off
# mol showrep 0 4 off
# #turn off "DNA" text
# mol off 3

# for {set i 0} {$i < 330 } {incr i 1} {
# 	display update
# 	rotate y by -1
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }

# turn on a new surface
# material add copy AOShiny
# mol modmaterial 7 0 Material$matnum
# material change opacity Material$matnum 0.0000
# for {set r 1} {$r < 16} {incr r 1} {
# 	display update
# 	rotate y by -1
# 	material change opacity Material[expr $matnum - 1] [expr {1 - double($r)/15}]
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }
# mol showrep 0 6 off
# mol showrep 0 7 on
# mol showrep 0 8 on

# for {set r 1} {$r < 16} {incr r 1} {
# 	display update
# 	rotate y by -1
# 	material change opacity Material$matnum [expr {double($r)/15}]
# 	material change opacity Transparent [expr {double($r)/15*0.15}]
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }
# incr matnum 1

#change text
# mol off 4 
# mol on 2

# add_text_layer SURFACE1
# set txtlncount 5
# draw color 16
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "Protein representation:" size 1.2 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 0 " "surface, probe radius 1.4A" size 1.2 thickness 1

# for {set i 0} {$i < 330 } {incr i 1} {
# 	display update
# 	rotate y by -1
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }

# material add copy AOChalky
# mol modselect 8 0 {protein and not (resname LYS  ARG and name  NH2 NZ) and not (resname GLU ASP and name  OD1  OE1)}
# mol showrep 0 9 on
# mol showrep 0 10 on
# mol modmaterial 9 0 Material$matnum
# mol modmaterial 10 0 Material$matnum
# for {set r 1} {$r < 31} {incr r 1} {
# 	display update
# 	rotate y by -1
# 	material change opacity Material$matnum [expr {double($r)/30}]
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }
# incr matnum 1

#charges text
# add_text_layer CHARGES
# set txtlncount 8
# draw color 16
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 100 " "Charged atoms" size 1.2 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 100 " "red - negative" size 1.2 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 100 " "(O in ASP, GLU)" size 1.2 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 100 " "(O in phosphate)" size 1.2 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 100 " "blue - positive" size 1.2 thickness 1
# incr txtlncount 1
# draw text " $txtx [expr $txty-($txtlncount*$txtstep)] 100 " "(N in LYS, ARG)" size 1.2 thickness 1
# incr txtlncount 1

# for {set i 0} {$i < 330 } {incr i 1} {
# 	display update
# 	rotate y by -1
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }

#enables DNA surface with charges
# material add copy AOShiny
# material add copy AOChalky
# mol modmaterial 11 0 Material$matnum
# mol modmaterial 12 0 Material[expr $matnum + 1]
# mol showrep 0 11 on
# mol showrep 0 12 on

# for {set r 1} {$r < 31} {incr r 1} {
# 	display update
# 	rotate y by -1
# 	material change opacity Material$matnum [expr {double($r)/30}]
# 	material change opacity Material[expr $matnum + 1] [expr {double($r)/30}]
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }

# for {set i 0} {$i < 330 } {incr i 1} {
# 	display update
# 	rotate y by -1
# 	render Tachyon dat/$fc.dat
# 	incr fc 1
# }
# incr matnum 2
exit
