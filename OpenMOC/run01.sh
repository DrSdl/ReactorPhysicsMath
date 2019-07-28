#!/bin/sh

max=100

for((i=1; i<=$max; ++i ))
do
    echo "$i"
    cp uox-assembly-ensemble-calculation_02.py template.py
    sed -i 's/drsdlID/'"$i"'/g' template.py
    python template.py
done
