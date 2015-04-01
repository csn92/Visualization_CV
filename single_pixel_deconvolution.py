import sima
import numpy as np
from sima import Sequence
import cv2
import math
from sklearn.linear_model import LassoCV, LassoLarsCV, LassoLarsIC


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
gamma = np.zeros(p) #Gamma values of neuron k of order p
s = np.zeros(K, T) #Spike Signal of kth neuron at time t

F = np.zeros(d,T) #Observed Fluorescence signal
alpha = np.zeros(K) #Spatial footprint of neuron k
base = np.zeros(T) #Baseline time-varying vector
noise = np.random.normal(0,1,T)

for k in range(1,K+1):
	#Calcium concentration
	for t in range(1,T+1):
		for p in range(1, p+1):
			c[k][t] += gamma[p]*c[k][t-i]
		c[k][t] += s[k][t]

	#Observed Fluorescence
	for t in range(1,T+1):
		Y[k][t] = alpha*c[k][t] + base + noise

	for i in range(1,d):
		find_parameter(c[i][:],  )

	#Calculate Auto Covariance
	for i  in range(1,d):
		N = size(c[i][:])
		Xs = np.average(c[i][:])
		C[i] = autocovariance(c[i][:], N, p, Xs)

	#Autocovariance
	for y in range(1,d):
		for tau in range(1,T):
			for j in range(1,p):
				C[y][tau] = C[y][tau]+gamma[j]*auto[y][tau-j]
			if tau<=p and tau>=1:
				C[y][tau] -= variance*gamma[tau]

def find_parameter(X,Y):
	model_bic = LassoLarsIC(criterion='bic')
	t1 = time.time()
	model_bic.fit(X,Y)
	t_bic = time.time() - t1
	alpha_bic = model_bic.alpha_

	model_aic = LassoLarsIC(criterion='aic')
	model_aic.fit(X,y)
	alpha_aic_ = model_aic.alpha_

def autocovariance(Xi, N, k, Xs):
    autoCov = 0
    for i in np.arange(0, N-k):
        autoCov += ((Xi[i+k])-Xs)*(Xi[i]-Xs)
    return (1/(N-1))*autoCov

def autocorrelation():
    return autocovariance(Xi, N, k, Xs) / autocovariance(Xi, N, 0, Xs)			

