import sima
import numpy as np
from sima import Sequence
import cv2
import math
import scipy.optimize.leastsq


NORMING_VAL = 2194
path ='/Users/Chithra/Documents/Columbia/Semester2/Test/z-series-004/z-series-004_Cycle00001_Element00001.h5'
# path ='/Users/Chithra/Documents/Columbia/Semester2/Test/TSeries-08242014-Day1-SessionA-000/TSeries-08242014-Day1-SessionA-000_Cycle00001_Element00001.h5'
sequence = Sequence.create('HDF5',path,'tzxyc')
print sequence.shape
(T, P, R, C, channel) = sequence.shape
d = P*R*C #Number of pixels
Y = np.zeros(d,T) #Number of observations
print T,d

#Initialisation of values
K = 10 #Number of neurons
p=1 #Order
c = np.zeros(K,T) #Concentration of neuron k at time t
gamma = 1 #Gamma values of neuron k of order p
s = np.zeros(K, T) #Spike Signal of kth neuron at time t

F = np.zeros(d,T) #Observed Fluorescence signal
a = np.zeros(K) #Spatial footprint of neuron k
B = np.zeros(T) #Baseline time-varying vector
noise = np.random.normal(0,1,T)

#Model calcium concentration
for i in range(1,K+1):
	for t in range(1, T+1):
		c[i,t] +=gamma*c[i,t-1] +s(i,t)

#Spatial calcium concentration
for t in range(1,T+1):
	for i in range(1,K+1):
		F[:,t] +=a[i]*c[i,t]
	F[:,t] += B[tau]
	
for t in range(1,T+1):
	Y[:,t] = F[:,t]+noise


def gaussian_kernel(a, tau, window,i):
	(taux,tauy,tayz) = tau.shape
	(wx,wy,wz) = window.shape
	for d in range(i-wx/2,i+wx/2):
		if d>0:
			Kx[d] = 1/a*exp(-pow(d,2)/(2*pow(taux,2)))
	for d in range(i-wy/2,i+wy/2):
		if d>0:
			Ky[d] = 1/a*exp(-pow(d,2)/(2*pow(tauy,2)))
	for d in range(i-wz/2,i+wz/2):
		if d>0:
			Kz[d] = 1/a*exp(-pow(d,2)/(2*pow(tauz,2)))

	Ky = np.concatenate((Kx,Ky.T),axis=1)
	return np.concatenate((Ky,Kz),axis=1)
	
def greedy_neuron_initialisation(Y,K, tau, window):
	(d,T) = size(Y)
	R = Y
	D = np.zeros((d,d))
	median = np.zeros(d)
	variance = np.zeros(d)
	center = np.zeros(K)

	for i in range(0,d):
		D[i] = gaussian_kernel(a, tau,window,i)

	for i in range(0,d):
		m = np.median(Y([i,:])
		for j in range(0,T):
			Y[i,j] = Y[i,j] - m

	

	for k in range(0,K):
		rho = D.T*R
		for i in range(0,d):
			variance[i] = sum(rho[i,:])
		R = reshape(variance,(,x,y,z))
		max1 = max(R)
		i=0
		for i in range(0,x):
			for j in range(0,y):
				for z in range(0,z):
					if R[i,j,z]==max1:
						arg = (i,j,z)

		center[k] = arg
		(x,y,z)  = arg
		set1 = set((0,0,0))
		for i in range(i-x/2,i+x/2):
			for j in range(j-y/2,j+y/2):
				for z in range(z-y/2,z+y/2):
					set1.update((i,j,k))
		set1.remove((0,0,0))
		center[k] = set1

		#Minimize value function

		return A,C,b,f

	



		



	