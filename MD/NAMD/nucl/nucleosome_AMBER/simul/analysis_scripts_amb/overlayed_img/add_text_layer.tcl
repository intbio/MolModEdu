#############################################################################################
#File version 1.0 07.03.12
#
#This file is a part of 
# Molecular Modeling Visualization Scripts Library
#
#Developed by Grigoriy Armeev and Alexey Shaytan
#Copyright (c) 2012 Molecular Simulations Group, Lomonosov Moscow State University
#
#for more information visit http://www.molsim.org
#mail at info@molsim.org
#############################################################################################
#
#
#
#Needs vars from input_param.tcl, needs vars molnum dispw transx transy transz
#
#---------------proc prepares a new molecule layer for text and graph----------------
#
#
#
proc add_text_layer {name} {
global molnum
global dispw
global transx
global transy
global transz
mol fix all
mol new
mol rename top $name
set top [molinfo top]
set sc [expr double(4)/$dispw]
set mat1 "$sc 0 0 0"
set mat2 "0 $sc 0 0"
set mat3 "0 0 $sc 0"
set mat4 "0 0 0 1"
lappend matrix $mat1 $mat2 $mat3 $mat4
lappend matrix1 $matrix
molinfo $top set rotate_matrix {{{1 0 0 0}{0 1 0 0}{0 0 1 0}{0 0 0 1}}}
molinfo $top set scale_matrix "$matrix1"
translate to 0 0 0
mol fix all
#translate by $transx $transy $transz
mol free all
for {set i [expr $molnum] } {$i < [molinfo num] } {incr i 1} {
mol fix $i
draw materials on
}
}
