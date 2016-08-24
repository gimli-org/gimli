#!/usr/bin/env bash

mkdir -p mesh

MESH=mesh/test2d

polyCreateWorld -d2 -x200 -y 100 -m1 $MESH

echo "-100 -20" > line; echo "100 -20" >> line
polyAddProfile -i line $MESH

polyAddVIP -R -z -25 -m2 $MESH

dctriangle -S -a 10 -v $MESH

