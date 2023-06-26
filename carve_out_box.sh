#!/bin/bash

echo "syntax: ./carve_out_box.sh filename.map x y z"


##A=(26 50 17)		#Subunit A: (26 50 17) y,x,z, len=(6 10 7) 
##B=(15 52 58)		#Subunit B: (15 52 58) y,x,z, len=(8 10 6)
##C=(-94 51 73)		#Subunit C: (-94 51 73) y,x,z, len=(7 9 6)
##D=(116 49 31)		#Subunit D: (-82 49 33) y,x,z, len=(7 9 7)

##var=$2[@]

##y=$(echo ${!var} | awk '{ print $1 }')
##x=$(echo ${!var} | awk '{ print $2 }')
##z=$(echo ${!var} | awk '{ print $3 }')


#coordinates of the central atom
#x=$(grep $3 $2 | tail -1 | awk '{ print $3 }')
#y=$(grep $3 $2 | tail -1 | awk '{ print $4 }')
#z=$(grep $3 $2 | tail -1 | awk '{ print $5 }')

###Fractional coordinate of a atom with respect to which the box will be carved out####

x=$2
y=$3
z=$4

#box half length or radius
ry=12
rx=12
rz=12


sftools <<eof
MAPIN $1
MAPLIMIT GRID  $(($x-$rx)) $(($x+$rx)) $(($y-$ry)) $(($y+$ry)) $(($z-$rz)) $(($z+$rz))
MAPOUT carved_$1
END
eof
