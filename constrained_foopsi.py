#Spike inference using a constrained foopsi approach
def constrained_foopsi(y,b,c1,g,sn,options):

# Inputs:

# 1. A,C,b,f using greedy initialization approach

# One variable for the coordinates.
# Another variable for the intensity values
def compute_coms(A,window,tau):
	(wx, wy) = window
	for a in A:
		for x_size in range(-wx/2,wx/2+1):
			for y_size in range()
			define=1


def median_filtering(A):
	A = scipy.ndimage.filters


def lars(X,y,sigma):

	n = len(y)
	sigma = sigma
	mu = np.zeros_like(y)
	beta = np.zeros((len(X.T), ), float)

# Compute center of mass for each spatial component and define search region for each pixel

	# min sum(sp)
	#subject to: sp>=0, b>=0, G*c = sp, c1>=0

	#||y-b-b-c_in|| <= sn*sqrt(T)

# Variables:
	# y: raw fluorescence data (vector of length(T))
	# c: denoised calcium concentration (T*1 vector)
	# b: baseline concentration (scalar)
	# c1: initial concentration (scalar)
	# g: discrete time constant
	# sn: noise standard deviation (scalar)
	# sp: spike vector (T*1 vector)
	options.noise_range=[0.25,0.5]
	options.noise_method='logmexp'
	options.lags=5
	options.resparse=0
	options.method='lars'
	options.p=2
	options.bas_noneng=1

	Ginv = [np.full(np.divide(G, scipy.sparse.speye(T)), np.ones(T, bas_est), gd_vec*ones(1,c1_est))]

	if bas_est:
		b=0
	if c1_est:
		c1=0
	if not sn:
		sn = GetSn(y,options.noise_range, options.noise_method)
	(A,B, spikes, C, D) = larse_regression_noise(y- b_lb *bas_est - b -c1*gd_vec, Ginv, 1, sh^2*T)
	sp = spikes[1:T]
	b = (spikes[T+bas_est] + b_lb)* bas_est + b*(1-bas_est)
	c1 = spikes*c1_est + c1*(1-c1_est)
	c = np.divide(G,sp)

	if not g:
		g = estimate_time_constants(y,p,sn,lags)
		while max(abs(np.roots(1,-g[:].tranpose()))>1):
			print "no stable AR"+str(p)+"model found. Checking for "+str(p+1)+"model"
			p = p+1
			g = estimate_time_constants(y,p,sn,lags)

		print "Stable AR model found"+str(p)

	if not bas_nonneg:
		b_lb=0
	else:
		b_lb = min(y)

	return c,b,c1,g,sn,sp

	
# Estimating time constant
def estimate_time_constants(y,p,sn,lags):
	# Estimate time constants from function

	lags = lags+p
	xc = np.cov(y,lags,'biased')
	xc = xc[:]
	A = scipy.linalg.toeplitz(xc(lags+(1:lags)), xc(lags+(1:p))) - sn^2*np.eye(lags,p)
	g = numpy.linalg.pinv(A)*xc(lags+2:end)

	return g

 
