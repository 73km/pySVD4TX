#!/bin/bash

#"usage: ./find_coordinates pdbfile resname subunit resid corresponding_ccp4_map"

#outfile=$(basename -s .pdb $1)_$3_$4.coordinate
#echo "$2 $3 $4" > $outfile
#echo "Orthogonal coordinates" >> $outfile

if [ $4 -gt 99 ]
then 	
	if grep "HETATM " $1 | grep -q "$2 $3 $4"; then grep "$2 $3 $4" $1 | grep "HETATM " | awk '{print $4, $3, $7, $8, $9}' >> dummy.ortho_coordinate
	elif grep "HETATM1" $1 | grep -q "$2 $3 $4"; then grep "$2 $3 $4" $1 | grep "HETATM1" | awk '{print $3, $2, $6, $7, $8}' >> dummy.ortho_coordinate
	else grep "ATOM"  $1 | grep "$2 $3 $4"| awk '{print $4, $3, $7, $8, $9}' >> dummy.ortho_coordinate
	fi
elif [ $4 -gt 9 ]
then
	if grep "HETATM " $1 | grep -q "$2 $3  $4"; then grep "$2 $3  $4" $1 | grep "HETATM " | awk '{print $4, $3, $7, $8, $9}' >> dummy.ortho_coordinate
	elif grep "HETATM1" $1 | grep -q "$2 $3  $4"; then grep "$2 $3  $4" $1 | grep "HETATM1" | awk '{print $3, $2, $6, $7, $8}' >> dummy.ortho_coordinate
	else grep "ATOM" $1 | grep "$2 $3  $4" | awk '{print $4, $3, $7, $8, $9}' >> dummy.ortho_coordinate
	fi
	#grep "$2 $3  $4" $1 | awk '{print $4, $3, $7, $8, $9}' >> dummy.ortho_coordinate 	
else
	if grep "HETATM " $1 | grep -q "$2 $3   $4"; then grep "$2 $3   $4" $1 | grep "HETATM " | awk '{print $4, $3, $7, $8, $9}' >> dummy.ortho_coordinate
	elif grep "HETATM1" $1 | grep -q "$2 $3   $4"; then grep "$2 $3   $4" $1 | grep "HETATM1" | awk '{print $3, $2, $6, $7, $8}' >> dummy.ortho_coordinate
	else grep "ATOM"  $1 | grep "$2 $3   $4"| awk '{print $4, $3, $7, $8, $9}' >> dummy.ortho_coordinate
	fi
	#grep "$2 $3   $4" $1 | awk '{print $4, $3, $7, $8, $9}' >> dummy.ortho_coordinate
fi


awk '{print $3, $4, $5}' dummy.ortho_coordinate > $(basename -s .pdb $1)_$2_$3_$4.ocoordinates
rm dummy.ortho_coordinate
