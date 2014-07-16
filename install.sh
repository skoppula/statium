#!/bin/bash -ex

#Unzipping library's interacting pair data
tar -zxvf data/ip_90_wGLY.tar.gz

#Unzipping library's PDB file data into culled_90 folder
mkdir culled_90
i=0; for i in {0..9}; do tar -zxf culled_90_$i.tar.gz; mv culled_90_$i/* culled_90/; done
