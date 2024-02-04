#!/bin/bash

builder builder.txt
diff $1.geo $1.geo.golden > junk
diff $1_modes.txt $1_modes.txt.golden >> junk
diff $1.proj $1.proj.golden >> junk

size="$(wc -c < junk)"

if [[ $size = 0 ]]; then
   echo "$1: pass"
else
   echo "$1: FAIL"
fi 

rm junk

