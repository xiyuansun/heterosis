#!/bin/awk -f
BEGIN {}
{
   a[FNR] = (a[FNR] ? a[FNR] FS : "") $i } END { for(i=1;i<=FNR;i++) print a[i] 
}
END { }

#BEGIN{}

#{out=""; 
#  for(i=1;i<=NF;i++){
#    out=$out" "$i}; 
#    print $out}
#    
#END{}