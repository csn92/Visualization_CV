from greedy_initialisation1 import greedyROI
from spatial_components import update_spatial_components
from temporal_components import update_temporal_components
from sima import Sequence
# import spams
import numpy as np
import matplotlib.pyplot as plt
import nimfa
import scipy
import scipy.sparse
from nitime import utils 
from nitime import algorithms as alg
from nitime.timeseries import TimeSeries
#
from PIL import Image

# path ='/Users/Chithra/Documents/Columbia/Semester2/Test/z-series-004/z-series-004_Cycle00001_Element00001.h5'
# path ='/Users/Chithra/Documents/Columbia/Semester2/Test/TSeries-08242014-Day1-SessionA-000/TSeries-08242014-Day1-SessionA-000_Cycle00001_Element00001.h5'
path = '/Users/Chithra/Documents/Columbia/Semester2/NMF/ctxA-002/ctxA-002_Cycle00001_Element00001.h5'
# NORMING_VAL = 2194
options={}
sequence = Sequence.create('HDF5',path,'tzxyc')
print sequence.shape
(T, P, R, C, channel) = sequence.shape
# For testing, only use the first 3 
P=1
# R=150
# C=150
T=500
d = R*C #Number of pixels
Y = np.zeros((R,C,T),dtype='float64') #Number of observations
for i in range(0,T):
    vol = sequence._get_frame(i).astype('float32')
    vol = np.nanmean(vol,axis=3)[0,:,:]
    print i
    Y[:,:,i] = vol

# nr = raw_input("Enter the number of neurons to be found?")
nr=1000
nr_size = 15
nr_std = 5

basis, Cin, center, data = greedyROI(Y[:,:,0:T],nr,nr_size,nr_std,C,R,T,P)
Ain = np.reshape(basis,(d,nr))
Cin = Cin.transpose()


# Correlation image
Cn = np.mean(Y, axis=2)
img = Image.fromarray(Cn)
img.show()

# Center of ROIs found from Initialization algorithm
implot = plt.imshow(Cn)
plt.scatter(x=center[:,0],y=center[:,1])
plt.show()

# Compute estimates of noise for every pixel and global time constant

# ff = np.nonzero(np.sum(Ain, axis=1))
# # Order of AR system
p=1
# options.update({'pixels':ff})
# # Y = Image.open('ctxA-002/')

# # img = Image.fromarray(Yr[:,:,20])
# # img.show()

Yr= np.reshape(Y,(d,T))
coefs_est = np.zeros((d,T))
for i in xrange(d):
    coefs_est[i], sigma_est = alg.AR_est_YW(Yr[i,:], p)
print coefs_est

print np.shape(np.dot(Ain,Cin))
nmf = nimfa.Nmf(np.maximum(Yr-np.dot(Ain,Cin),0), rank=1)
nmf_fit = nmf()
bin = nmf_fit.basis()
fin = nmf_fit.coef()
print bin,fin

print "Finished doing NMF",

# A, b = spams.lasso(Yr-np.dot(bin,fin),pos=True,D=Cin,lambda2=np.multily(coefs_est,np.power(T,2)))
# print A,b
# Updating spatial components
# A,b = update_spatial_components(Yr,Cin,fin,Ain,P,C,R,coefs_est)

# # Updating temporal components
# C,f,Y_res = update_temporal_components(Yr,A,b,Cin,fin,P)


# # Merge found components
# A_in = A
# C_in = C
# repeat = True
# A_ = A
# C_ = C
# P['method'] = 'project'
# P['merge_thr'] = 0.8
# while repeat:
# 	A,C,nr,merged_ROIS = merge_ROIs(Y_res, A_, b, C_,f,P)
# 	if not merged_ROIS():
# 		repeat = False
# 	print nr
# 	A_ = A
# 	C_ = C

# P['method'] = 'constrained_foopsi'
# C2, f2, Y_res = update_temporal_components(Yr,A_or,b,C_or,f,P)


# print A,C2


