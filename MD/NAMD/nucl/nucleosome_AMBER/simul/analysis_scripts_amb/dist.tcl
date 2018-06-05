proc mdsm { selstring1 selstring2 } { 

set sel1 [atomselect top $selstring1] 
set sel2 [atomselect top $selstring2] 
 
  
set cent_sel1 [measure center $sel1] 
set cent_sel2 [measure center $sel2] 
  
 
set dist [veclength [vecsub $cent_sel1 $cent_sel2]] 
 
$sel1 delete 
$sel2 delete 
return $dist 
}
