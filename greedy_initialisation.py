import scipy.ndimage as ndimage
import numpy as np
import numpy.linalg
from sima import Sequence
import scipy
import math
import gc
import scipy.optimize as optimize
import Utilities

R=0
C=0
P=0
T=0
# path ='/Users/Chithra/Documents/Columbia/Semester2/Test/z-series-004/z-series-004_Cycle00001_Element00001.h5'
path ='/Users/Chithra/Documents/Columbia/Semester2/Test/TSeries-08242014-Day1-SessionA-000/TSeries-08242014-Day1-SessionA-000_Cycle00001_Element00001.h5'
NORMING_VAL = 2194
sequence = Sequence.create('HDF5',path,'tzxyc')
print sequence.shape
(T, P, R, C, channel) = sequence.shape
# For testing, only use the first 3 
T=5
P=1
R=150
C=150
d = R*C #Number of pixels
Y = np.zeros((d,T)) #Number of observations

# initialising with values of signals
for i in range(0,T):
	vol = sequence._get_frame(i).astype('float32')

	# vol /= NORMING_VAL
	# vol = np.clip(vol,0,1)
	vol = np.nanmean(vol,axis=3)[0,0:150,0:150]
	# print "Shape of vol ", np.shape(vol),d
	tau_x = np.var(vol,axis=0)
	tau_y = np.var(vol,axis=1)
	# tau_z = np.var(vol,axis=2)
	vol = vol.flatten()
	# print "Shape of vol", np.shape(vol), d

	Y[:,i] = vol

gaussian_filter= np.zeros((d,d))
cells = np.zeros((d,d))
# window = (8,8,4)



# import np.meshgrid, np.linspace
# Estimating the variance of noise
def getSn(Y,range_ff,method):

	# Estimate noise level with power spectral density method
	L = Y.flatten()

	psd_Y,ff = scipy.signal.welch(L,window=np.round(L/8),fs=1000)
	ind=[0]*len(ff)
	for i, item in enumerate(ff):
		if item>range_ff[0] and item<range_ff[1]:
			ind[i]=1
		else:
			ind[i]=0
	
	if method=="mean":
		sn = np.sqrt(np.mean(psd_Y(ind)/2))
	elif method=="median":
		sn = np.sqrt(np.median(psd_Y(ind)/2))
	elif method=="logmexp":
		sn = np.sqrt(np.exp(np.mean(np.log([psd_Y[i]/2 for x in ind]))))

	print "Sigma is " + str(sn)
	return sn
# get_frames(Y, K,tau,window)

def gaussian_kernel(X,tau,window,Y):
	ndimage.filters.gaussian_filter(X,tau,output=Y,truncate=window)
# Define normal distribution kernel
def gaussian_blur(tau,window):
	(wx,wy,wz) = window
	(taux,tauy,tauz) = tau
	
	gau1 = np.random.normal(0.0,np.power(taux,2),wx)
	gau2 = np.random.normal(0.0,np.power(tauy,2),wy)
	gau3 = np.random.normal(0.0,np.power(tauz,2),wz)

	return [gau1,gau2,gau3]


def define_gaussian_blur(tau,window,d,sigma,i):
	# (wx,wy,wz) = window
	(wx,wy) = window
	# x,y,z = get_coordinates(i)
	x,y = get_coordinates(i)
	# print i, x,y

	for x_size in range(-wx/2,wx/2+1):
		for y_size in range(-wy/2,wy/2+1):
			# for z_size in range(-wz/2,wz/2):
				x_tot = x+x_size
				y_tot = y+y_size
				# z_tot = z+z_size
				# if i<5:
				# 	print "Window sizes ", x_size, y_size, x_tot,y_tot, x,y
				if x_tot>=0 and (x_tot<R and (y_tot<C and y_tot>=0)):
				# if x_tot>0 and (x_tot<R and (y_tot<C and (y_tot>0 and (z_tot>0 and z_tot<P)))):
					# dist = math.pow(x_size,2)+math.pow(y_size,2)+math.pow(z_size,2)
					
					dist = math.pow(x_size,2)+math.pow(y_size,2)
					gaussVal = math.exp(-dist/(2*sigma*sigma))
					# i = find_pixel(x+x_size,y+y_size,z+z_size)
					j = find_pixel(x_tot,y_tot)
					cells[i,j] = 1
					gaussian_filter[i,j] = gaussVal


def get_frames(Y, K,tau,window,sigma,T):
	R = Y

	(d,T) = np.shape(Y)
	global gaussian_filter
	#Define Gaussian Blur Matrix D
	# tau = (1,1,1)
	for i in range(0,d):
		define_gaussian_blur(tau,window,d,sigma,i)
		if i%1000==0:
			gc.collect()
		
	
	#Subtract and store median value for each pixel
	# m = np.zeros(d)
	# for i in range(0,d):
	m=np.median(Y,axis=1)
	# m1 = np.repeat(m[0,i],m[i],axis=1)
	print "finished finding the median value ", m[d-1], m[1],m[3]
	for i in range(0,d):
		R[i,:] = R[i,:] - m[i]
	print "Finished subtracting the median value"

	A = np.zeros(K)
	C = np.zeros(T)
	v = np.zeros(d)
	rho = np.zeros((d,T))
	# Calculate variance explained by each kernel
	for k in range(0,K):
		print "Trying to find rho"
		i=0
		
		# rho = np.dot(gaussian_filter,R)
		for i in range(0,d):
			
			rho[i,:] = np.dot(gaussian_filter[i,:],R)
			# print i,gaussian_filter.shape, gaussian_filter[0,:]
			# # gaussian_filter = np.reshape(gaussian_filter[1:d]).reshape
			# gaussian_filter = np.delete(gaussian_filter,0,axis=0)
			# print gaussian_filter.shape
			# if gaussian_filter.size==0:
			# 	flag=True
			# print i
			if not (i%1000):
				gc.collect()
			i+=1
		# for i in range(0,d):
		# 	print i
		# 	for j in range(0,T):
		# 		# print i,j
		# 		rho[i,j] = np.dot(gaussian_filter[i,:],R[:,j])
		

		print "Fininshed finding the rho(density) "
		for i in range(0,d):
			v[i] = sum(rho[i,:])

		#Finding centers
		max1 = max(v)
		# Get center pixel
		bp = np.where(v == max1)
		bp = int(bp[0])
		print bp
		# Get center coordinates
		cx,cy= get_coordinates(bp)

		ak = np.zeros(d)
		ck = np.zeros(T)

		# For the bounds
		bnds=[]
		for i in cells[bp,:]:
			if cells[bp,i]==0:
				bnds.append((0,0))
			else:
				bnds.append((0,None))
		# for i in range(0,T):
		# 	bnds.append((None,None))
		# for i in range(0,d):
				# cons=({'type' : 'eq',
				#  'fun' : lambda x: np.array([x[0]])})
		# cons=({'type' : 'ineq',
		#  	   'fun' : lambda x: np.array(x),
		#  	    })
		# cons=(lambda x: x[0])
		x0 = ak.ravel()
		# x0 = np.append(x0,ck.ravel(),axis=0)
		x1 = ck.ravel()
		result = optimize.fmin_slsqp(func, (x0) ,args=(R,x1),f_ieqcons=cons,bounds=bnds)

		R -= np.dot(ak,ck)
	A.append(ak)
	C.append(ck)

	R +=m
	cons1 = (
			{'type' : 'ineq',
				 'fun' : lambda x: x[0], 'args':(R)},

			{'type' : 'ineq',
			'fun' : lambda x : np.array(x[1])}
			)
	
	result = optimize.fmin_slsqp(func2, args=[[b,f],R], method='SLSQP', constraints=cons1, bounds=bnds)
	
	return A,C,b,f

def cons(x,*args):
	return x[0]

def func(x,*args):
	ak = x
	R = args[0]
	ck = args[1]
	result = np.zeros((d,T))
	for i in ck:
		result[i,:] = ck*i
	# print R-result

	return R-result

def func2(x,R):
	return R-(np.dot(x[0],x[1]))



# print window, sigma
window=(3,3)
sigma = getSn(Y,[0.25,0.5],'logmexp')
tau = (tau_x,tau_y)
K=10
A,C,b,f = get_frames(Y,K,tau,window,sigma,T)






























