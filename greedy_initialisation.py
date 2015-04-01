import scipy.ndimage as ndimage
import numpy as np
import scipy.optimize.minimize
import numpy.linalg.norm

def gaussian_kernel(X,tau,window,Y):
	ndimage.filters.gaussian_filter(X,tau,output=Y,truncate=window)

def get_frames(Y, K,tau,window):
	R = Y
	(d,T) = Y.size

	#Define Gaussian Blur Matrix D
	D = np.zeros((d,d))
	gaussian_kernel(Y(:,))
	m = np.zeros(d)
	v = np.zeros(d)
	centers = np.zeros(K)
	ak = np.zeros(d)
	ck = np.zeros(T)
	constraints={}
	for i in range(0,d):
		m[i]=np.median(Y[i,:])
	
	#Subtract median
	for t in range(0,T):
		for i in range(0,d):
			Y[i,t] = Y[i,t] - m[i]

	for k in range(0,K):
		rho = D.transpose()*R
		for i in range(0,d):
			v[i] = rho[i,:]

		#Finding centers
		max1 = max(v)
		for i in range(0,x):
			for j in range(0,y):
				if R[i,j,z]==max1:
					arg = (i,j,z)

		centers[k] = arg
		(x,y)  = arg
		set1 = set((0,0))
		for i1 in range(x-wx/2,x+wx/2):
			for j1 in range(y-wy/2,y+wy/2):
				if i1>0 and j1>0 :
					set1.update((i1,j1))
		set1.remove((0,0))
		

		func = R - ak*ck.transpose()
		func = norm*(func)
		constraints.update({ak:})

		minimize(func,)


		
































path ='/Users/Chithra/Documents/Columbia/Semester2/Test/z-series-004/z-series-004_Cycle00001_Element00001.h5'
#path ='/Users/Chithra/Documents/Columbia/Semester2/Test/TSeries-08242014-Day1-SessionA-000/TSeries-08242014-Day1-SessionA-000_Cycle00001_Element00001.h5'
NORMING_VAL = 2194
sequence = Sequence.create('HDF5',path,'tzxyc')]
(T, P, R, C, channel) = sequence.shape
d = R*C #Number of pixels
Y = np.zeros(d,T) #Number of observations