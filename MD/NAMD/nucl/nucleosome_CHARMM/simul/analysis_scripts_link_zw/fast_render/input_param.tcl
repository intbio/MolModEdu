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
# Initiate vars: $dispw $disph
#
#
#Input parameters for movie.tcl
# Here we put some initial variables such as display width&height in pixels scaling for molecule and parts of translation vector.

display rendermode GLSL


set scale 1.75
set dispw 1000
set disph 750
set transx 0.3
set transy -0.40
set transz 0
# Display settings

display antialias on


display resize $dispw $disph

display eyesep       0.960000
display focallength  6.000000
display height       6.000000
display distance     -2.000000
display projection   Orthographic
display nearclip set 0.500000
display farclip  set 10.000000
display depthcue off

# I wish this lines will make some difference for tachyon
display ambientoclussion on
display shadows on
display antialias on
display aoambient 0.75
display aodirect 0.3

color change rgb 24 0.15 0.25 0.93
color change rgb 30 0.81 0.18 0.18
