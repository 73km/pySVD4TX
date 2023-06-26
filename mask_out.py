#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 16:37:36 2021
@author: tmalla
"""
from shutil import copy
import mrcfile
import numpy as np
import math as m
import os
import sys

###############################input parameters here#################

input_map_file=sys.argv[1]       #original ccp4 map file for the required timepoint
input_coordinate_file=sys.argv[2] #orthogonal coordinate file for residue of interest
masked_out_map=input("Give name for output map:\n")        

r = 3  # mask radius around atom

#####################################################################
input_map=mrcfile.open(input_map_file, mode='r')
coordinates=np.loadtxt(input_coordinate_file,dtype=np.float32)

if not os.path.exists(masked_out_map):
    copy(input_map_file,masked_out_map)

(a,b,c)=np.array([input_map.header.cella.x,input_map.header.cella.y,input_map.header.cella.z])
(alpha,beta,gamma)=np.array([input_map.header.cellb.alpha,input_map.header.cellb.beta,input_map.header.cellb.gamma])
mx = input_map.header.mx
my = input_map.header.my
mz = input_map.header.mz
#print(mx,my,mz)
input_mapdata=input_map.data


##############################################################3
def Orthogonalization_Matrix(a,b,c,alpha,beta,gamma):
    #x = M*X; x-fractional coordinate system x= X/a, X is Cartesian coordinate vectors
    #M - Deorthogonalization matrix; 
    #[M^-1 - Orthogonalization Matrix]; X = [M^-1]x
   
    alpha = alpha*(m.pi/180); beta = beta*(m.pi/180); gamma = gamma*(m.pi/180)
    V = a*b*c*m.sqrt(1-(m.cos(alpha))**2-(m.cos(beta))**2-(m.cos(gamma))**2 
                     + 2*m.cos(alpha)*m.cos(beta)*m.cos(gamma))
    M11 = a
    M12 = b*m.cos(gamma)
    M13 = c*m.cos(beta)
    M21 = 0.0 
    M22 = b*m.sin(gamma)
    M23 = c*(m.cos(alpha)-m.cos(beta)*m.cos(gamma))/(m.sin(gamma))
    M31 = 0.0
    M32 = 0.0
    M33 = V/(a*b*m.sin(gamma))
    M = np.array(([M11,M12,M13],[M21,M22,M23],[M31,M32,M33]))
    return np.around(M,decimals=5)

def O2F(X,Y,Z): #Converts orthogonal coordinates to fractional coordinates
    XF,YF,ZF=np.dot(M_inverse,np.array([X,Y,Z])) 
    return XF,YF,ZF

M=Orthogonalization_Matrix(a,b,c,alpha,beta,gamma)
M_inverse=np.linalg.inv(M)

n=(len(coordinates[:,0]))  #This is number of atoms in coordinates file.
#zero matrix that is same shape as coordinate and where new fractional coordinates will be stored
Frac_coordinates=np.zeros(np.shape(coordinates)) 

for i in range(n):
    Frac_coordinates[i,0],Frac_coordinates[i,1],Frac_coordinates[i,2]=O2F(coordinates[i,0],coordinates[i,1],coordinates[i,2])
        

#Multiply by mx, my and mz to obtain the exact coordinate inside the unitcell box. 
#Since these are positions, I convert them into integers. And these values could 
#be negative which are dealt later.     
mxes=np.array([mx,my,mz])
map_coordinates=np.round(Frac_coordinates * mxes).astype(int)

#From given input r, corresponding length in volxels is calculated.
cons = np.round(my/a) * r
boxlen = (2 * cons + 1).astype(int)


#Iterating through each atom, coordinates of cube with dimension boxlen is calculated
#with the atom coordinate at the center of cube.
def func_box_coordinates(points):
    point = points.tolist()
    x = np.arange(point[0]-cons,point[0]+cons+1)
    y = np.arange(point[1]-cons,point[1]+cons+1)
    z = np.arange(point[2]-cons,point[2]+cons+1)
    ind = []
    for i in x:
        for j in y:
            for k in z:
                xyz=np.array((i,j,k))
                if (m.dist(point,xyz)) <= cons: #this line carves out a sphere around atoms, else a box
                    ind.append(xyz)
    return(np.asarray(ind))


box_coordinate_list = []

for i in range(n):
    point=np.array(map_coordinates[i,:]) 
    box_coordinate = func_box_coordinates(point)
    box_coordinate_list.append(box_coordinate)

total_grid_points=np.prod(np.asarray(box_coordinate_list).shape)/3
#The box coordinates come out as column array. All of them are combined and the
#duplicates are removed. Since the atoms are so close, of course the voxels will repeat.
#bbb array is a tall array with each row being a coordinate. ccc array removes duplicates.    
bbb=np.asarray(box_coordinate_list).reshape(total_grid_points.astype(int),3)
ccc=np.unique(bbb.astype(int),axis=0)
ddd=np.copy(ccc)

#Negative values mean the coordinates are outside of the box. But because of crystallographic symmetry, 
#adding mx,my or mz to the coorresponding coordinate value shifts the position inside the box. Because the 
#map is everything that's inside the unit cell, this is important.
for i in range(len(ccc[:,0])):
    if ddd[i,0] < 0:
        ddd[i,0] += mx
    if ddd[i,1] < 0:
        ddd[i,1] += my
    if ddd[i,2] < 0:
        ddd[i,2] += mz
  

#A is an array whose shape is same as the original map data array with all elements zero
A = np.zeros(input_mapdata.shape,dtype=np.float32)

#Change the order of (x,y,z) here to get the correct axes. Crystallographic ccp4 maps are 
#notorious and they change all the time.
y=np.array(ddd[:,0],dtype=np.intp)
x=np.array(ddd[:,1],dtype=np.intp)
z=np.array(ddd[:,2],dtype=np.intp)

#Now only the coordinates falling inside our mask box are updated with original values from 
#corresponding position in original map file. The result is a map with values zero everywhere
#except inside our desired box. 
A[x,y,z] = input_mapdata[x,y,z]
#A[x,y,z] = 1   #use this value instead to generate the mask for display.

input_map.close

#Writing a new map file with these values. 
masked_map = mrcfile.open(masked_out_map, mode='r+')
masked_map.set_data(A)
masked_map.header.mx = mx
masked_map.header.my = my
masked_map.header.mz = mz
masked_map.close

print('done')

