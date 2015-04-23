
def greedyneuron1D(Y,K,tau,window):
	R=Y

	

K = raw_input("Enter the number of neurons")
print "The number of neurons entered is %s." %K

path ='/Users/Chithra/Documents/Columbia/Semester2/Test/z-series-004/z-series-004_Cycle00001_Element00001.h5'
#path ='/Users/Chithra/Documents/Columbia/Semester2/Test/TSeries-08242014-Day1-SessionA-000/TSeries-08242014-Day1-SessionA-000_Cycle00001_Element00001.h5'
NORMING_VAL = 2194
sequence = Sequence.create('HDF5',path,'tzxyc')
print sequence.shape
(T, P, R, C, channel) = sequence.shape

# For testing, only use the first 100 and 12 planes
T=100
P=12
d = P*R*C #Number of pixels
Y = np.zeros((d,T)) #Number of observations

# initialising with values of signals
for i in range(0,T):
	vol = sequence._get_frame(i).astype('float32')
	# vol /= NORMING_VAL
	# vol = np.clip(vol,0,1)
	vol = np.nanmean(vol,axis=3)[0:12,:,:]
	tau_x = np.var(vol,axis=0)
	tau_y = np.var(vol,axis=1)
	tau_z = np.var(vol,axis=2)
	vol = vol.flatten()
	Y[:,i] = vol

window = (25,25,P)
# sigma = getSn(Y,[0.25,0.5],'logmexp')

# Estmate the paramteres of variance
tau = (tau_x,tau_y,tau_z)
greedyneuron1D(Y,K,tau,window)
