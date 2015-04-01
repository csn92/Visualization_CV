import numpy as np

T=1
c_in = np.zeros((1,T))
c = np.zeros((1,T))
s = np.zeros((1,T))
gamma = np.zeros((1,p))

G = np.zeros((T,T))

#G for spike
for i in range(0,T):
	for j in range(i,0,-1):
		if i==j:
			G[i,j] = 1
		else:
			G[i,j] = gamma[1,j]

G*(c - c_in) = s
G = np.