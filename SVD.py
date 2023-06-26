#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 12:14:46 2021

@author: tmalla
"""
import mrcfile
import numpy as np
import os
from scipy.linalg import svd
from scipy.optimize import curve_fit
from shutil import copy2
import time
import matplotlib.pyplot as plt 
import warnings
import csv

warnings.filterwarnings('ignore')

  

##################################################################

def SVD_TEK(svdA):
    # I use QR decomposition method. There are other methods available.

    Q, R = np.linalg.qr(svdA)   
    svdu , svds, svdvt = svd(R,lapack_driver='gesvd')
    svdU=Q@svdu
    svdV=np.transpose(svdvt)
    return svdU,svds,svdV

#####################################################################
def write_map(refmap,output_map,data):
    copy2(refmap,output_map)
    outmap = mrcfile.open(output_map,mode='r+')
    inmap = mrcfile.open(refmap)
    outmap.set_data(np.array(data).reshape(inmap.data.shape))
    # mmean = mrcfile.open('input_files/'+str(map_files[i])+maptype+mapformat).header.dmean
    # rrms = mrcfile.open('input_files/'+str(map_files[i])+maptype+mapformat).header.rms
    # mx = inmap.header.mx
    # my = inmap.header.my
    # mz = inmap.header.mz
    outmap.header.mx = inmap.header.mx
    outmap.header.my = inmap.header.my
    outmap.header.mz = inmap.header.mz
    # modifiedfile.header.dmean = mmean
    # modifiedfile.header.rms = rrms
    outmap.close
    inmap.close
    
#######################################################################

    
t=time.time()

carved='carved_'
masked='masked_'
identifier='D_'
mapformat='.map'

map_files = [
    'water',
    '3ms',
    '6ms',
    '15ms',
    '30ms',
    '66ms',
    '240ms',
    '700ms'
    ]

print('\nInput files:\n' + str([carved + masked + identifier + sub + mapformat for sub in map_files]))

n=len(map_files)


for i in range(n):
    vars()["map_"+str(i)] = mrcfile.open('input_files/'+carved+masked+identifier+str(map_files[i])+mapformat,mode='r')
    #shapes=mrcfile.open('input_files/'+carved+masked+identifier+str(map_files[i])+mapformat,mode='r').data.shape
    shapes=(locals()['map_'+str(0)]).data.shape
    if (locals()['map_'+str(0)]).data.shape == (locals()['map_'+str(i)]).data.shape :
        m = np.prod(shapes)
        A = np.zeros((m,n),dtype=np.float32)
    else :
        raise(Exception("\nDimension of map #:"+str(i+1), (locals()['map_'+str(i)]).data.shape,", does not match with that of previous ones:", (locals()['map_'+str(i-1)]).data.shape))
    
if not os.path.exists('output_files'):
    os.makedirs('output_files')

ind = []
C = np.copy(A)
print("\nChecking ±3 rmsd condition")
for i in range(n):
    wholemap = mrcfile.open('input_files/'+str(map_files[i])+mapformat, mode='r')
    mean_wholeA = wholemap.header.dmean
    rms_wholeA = wholemap.header.rms
    A[:,i] = (locals()['map_'+str(i)]).data.reshape(-1)
    dummyA = np.copy(A[:,i])
    vars()["map_"+str(i)].close()
    cond=(np.where((dummyA <= (mean_wholeA - 3 * rms_wholeA)) | (dummyA >= (mean_wholeA + 3 * rms_wholeA))))
    ind = np.append(ind,np.asarray(cond))
    wholemap.close

iind = np.unique(np.asarray(ind).astype(int))

C[iind,:] = A[iind,:]
     
#Uncomment the following lines to save these modified maps.

# print("Saving evolved maps as Modified* inside output_files.")
# for i in range(n):
#     modified_inmap = 'input_files/'+carved+masked+identifier+str(map_files[i])+mapformat
#     modified_outmap = 'output_files/Modified-'+carved+masked+identifier+str(map_files[i])+mapformat
#     modified_data = np.array(C[:,i])
#     write_map(modAified_inmap,modified_outmap,modified_data)

    
######################################Computing SVD###########################################################

print("\nNow computing SVD:")

(U1,S1,V1) = SVD_TEK(C)
print(S1)

#Uncomment the following lines to print the left singular vector maps

# print("Writing LSVs into maps")
# for i in range(n):
#     svd_inmap = 'input_files/'+carved+masked+identifier+str(map_files[0])+mapformat
#     svd_outmap = 'output_files/LSV'+str(i+1)+'_'+identifier+mapformat
#     svd_data = np.array(U1[:,i])
#     write_map(svd_inmap,svd_outmap,svd_data)


##############Recreating original maps from first 3 significant values#########################################

print('\nRecreating original map from 3 significant vectors')
uu=U1[:,[0,1,2]]
ss=np.diag(S1[[0,1,2]])
vv=np.transpose(V1[:,[0,1,2]])
AA=uu@ss@vv

#Uncomment the following lines to save the recreated maps from significant values

# print("Writing recreated maps from 3 SVs")
# for i in range(n):
#     recreated_inmap = 'input_files/'+carved+masked+identifier+str(map_files[i])+mapformat
#     recreated_outmap = 'output_files/Recreated-'+carved+masked+identifier+str(map_files[i])+mapformat
#     recreated_data = np.array(AA[:,i])
#     write_map(recreated_inmap,recreated_outmap,recreated_data)

##################################################################

(U,S,V) = SVD_TEK(AA)

with open(identifier+'.txt','w') as rsvfile:
    csv.writer(rsvfile,delimiter=' ').writerows(V)
rsvfile.close
print(V)


print("\nProgram completed successfully in %.2f seconds" %(time.time()-t))     
