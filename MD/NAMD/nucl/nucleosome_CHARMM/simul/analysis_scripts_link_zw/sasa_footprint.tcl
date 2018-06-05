set radius 1.4
set nuclchain1 "segname CHI "
set nuclchain2 "segname CHJ "

 mol load psf ../analysis_data/only_nucl_init.psf
 mol addfile ../analysis_data/only_nucl_init.pdb
 mol ssrecalc top
 mol addfile ../analysis_data/md_nucl.dcd step 100 waitfor all 


# mol new ../input/1kx5nt_ready.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
# mol addfile md.dcd type dcd first 0 last -1 step 100 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation Lines

set filename "../analysis_data/footprint.txt"

set fileId [open $filename "w"]

set sel1 [atomselect top "$nuclchain1" frame 0]
set sel2 [atomselect top "$nuclchain2" frame 0]
set ids1 [[atomselect top "$nuclchain1 and name P"] get resid]
set ids2 [[atomselect top "$nuclchain2 and name P"] get resid]
set result1 [dict create]
set result2 [dict create]
foreach a $ids1 {
dict set result1 $a 0
}
foreach a $ids2 {
dict set result2 $a 0
}

for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
    puts "Processing $i frame"
    foreach a $ids1 {
        set tempsel [atomselect 0 "$nuclchain1 and resid '$a' and name O3' P" frame $i]
        set sasa [measure sasa $radius $sel1 -restrict $tempsel]
        dict set result1 $a [expr [dict get $result1 $a] + $sasa]
        $tempsel delete
    }
    puts "$nuclchain1 done"
    

    foreach a $ids2 {
        set tempsel [atomselect 0 "$nuclchain2 and resid '$a' and name O3' P" frame $i]
        set sasa [measure sasa $radius $sel2 -restrict $tempsel]
        dict set result2 $a [expr [dict get $result2 $a] + $sasa]
        $tempsel delete
    }
    puts "$nuclchain2 done"
}

foreach a $ids1 {
dict set result1 $a [expr [dict get $result1 $a] / double([molinfo top get numframes])]
}
foreach a $ids2 {
dict set result2 $a [expr [dict get $result2 $a] / double([molinfo top get numframes])]
}

puts $fileId "# resid $nuclchain1 $nuclchain2"

foreach a $ids1 {
puts $fileId "$a [dict get $result1 $a] [dict get $result2 $a]"
}

close $fileId

