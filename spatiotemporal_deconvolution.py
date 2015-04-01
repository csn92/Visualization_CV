import numpy as np 

# Variables
	# c: Concentration
	# s: spike
	# F: Fluorescence
	# a: alpha
	# B: baseline
	# eta: Noise	
	# K: Number of neurons
	# T: Time Period
	# Gamma: decay constant

C = np.array((K,T+1)) 
S = np.array(K,T+1))
F = np.array((T+1,1))
a = np.array((K,1))
B = np.array((T,1))
Y = np.array((T,1))
noise = #White noise

# Autoregressive dynamics
for i in range(0,K):
	for t in range(1,T+1):
		c[i,t] = gamma*C[i,t-1] + s[i,t]

# Spatial Calcium Concentration
for t in range(1,T+1):
	F[t,1] = (a[:,1]*c[:,t]).sum() + B[t,1]

# corruption by additive Gaussian noise
for t in tange(1,T+1):
	Y[t,1] = F[t,1] + noise
	
# Estimating A, b





