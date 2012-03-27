#!/usr/bin/env bash

mkdir -p mesh
MESH=mesh/world2d
QUALITY=34

polyCreateWorld -d2 -m1 $MESH
#polyTranslate -z -0.25 $MESH
polyScale -x 30 -y 10 $MESH

nSegment=24
rm -r circ.txt
for i in `seq 1 $nSegment`; do 
    echo $nSegment $i |awk '{ print (sin($2*(2.*3.141592)/$1) " " cos($2*(2.*3.141592)/$1)) }' >> circ.txt
done

polyCreateWorld -d2 -C -t circ.txt -s 1 -m2 circ

# Radius 2m
polyScale -x 2 -y 2 circ
# Mittelpunktstiefe 5m
polyTranslate -y -5 circ

polyMerge $MESH circ $MESH

polyTranslate -x 5 circ

polyMerge $MESH circ $MESH

polyConvert -V -o $MESH-poly.vtk $MESH.poly

dctriangle -v -'q'$QUALITY -S $MESH

