import scipy.ndimage

def spationtemporal_spike_inference(A,C,b,f):

	while convergence:

		# Calculate center of mass of spatial components
		measurements.center_of_mass(A)
		dual(X,Y,sn,b,f)
		larse_regression_noise()


def dual(A,y,sn,b,f,bas_est,c1_est):
	
	thr = sn*math.sqrt(T)
	lagrangian_temporal_gradient(Ald, math.pow(thr,2), y-A*c -b*f,bas_est,c1_est)

	sklearn.linear_model.Lasso(alpha=1.0, fit_intercept=True, normalize=False, precompute=False, copy_X=True, max_iter=1000, tol=0.0001, warm_start=False, positive=False, random_state=None, selection='cyclic')
	ld_in=10

def lars(A,y,sn,b,bas_est,c1_est,T):
	if not bas_est:
		b=0
	if not c1_est:
		c1=0

	A,b = lars_regression_noise(y-A*c-b*f, math.pow(sn,2)*T)

def larse_regression_noise(Y,X,positive, noise):
	k=1;
	N = size(X,2)
	maxcomps = N;

	W= np.zeros(N,k)
	active_set = np.zeros((N,k))
	visited_set = np.zeros((N,k))

	lambdas=[]

	Ws = np.zeros((N,K))

	r = X.transpose()* Y(:);
	M = -X*X;

	i=1
	flag=0
	while 1:
		if flag==1:
			W_lam=0
			break

		if i>1 and new and visited_set(new)==0:
			visited_set(new) = 1

		# Compute Gradient
		dQ = r+M*W

		# Compute new W
		if i==1:
			if positive:
				dQa = dQ
			else:
				dQa = abs(dQ)
			(lamb, new) = max(dQa(:))

			if lamb<0:
				print "All negative directions!"
				break
		else:
			# Calculate vector to travel along

			(avec, gamma_plus, gamma_minus) = calcAvec(new, dQ, W, lambda, active_set, M, positive)

			if new==0:
				if dropped_sign==1:
					gamma_plus(dropped) = inf
				else:
					gamma_minus(dropped) = inf

			gamma_plus(active_set == 1 ) = inf
			gamma_plus(gamma_plus <=0 ) = inf
			gamma_plus(gamma_plus > lamb) = inf
			(gp_min, gp_min,ind) = min(gamma_plus[:]);

			if positive:
				gm_min = inf
			else:
				gamma_minus(active_set==1) = inf
				gamma_minus(gamma_minus > lamb) = inf
				gamma_minus(gamma_minus <=0 ) = inf
				(gm_min, gm_min_ind) = min(gamma_minus[:])


			if g_min == inf:
				g_min = lamb

			# LARS check 

			gamma_zero = -W(active_set == 1) /avec
			gamma_zero_full = zeros(N,k)
			gamma_zero_full (active_set==1) = gamma_zero
			gamma_zero_full(gamma_zero_full<=0) = inf
			(gz_min, gz_min_ind) = min(gamma_zero_full[:])

			if gz_min<g_min:
				if verbose:
					print "Dropping active weight:"+gz_min_ind

				active_set(gz_min_ind) = 0
				dropped = gz_min_ind
				dropped_sign = sign(W(dropped));
				W(gz_min_ind) =0

				avec = avec(gamma_zero!=gz_min)
				g_min = gz_min_indnew=0

			elif g_min<lamb:
				if which==1:
					new = gp_min_ind
					if verbose:
						print "New positive component: "+new
					elif:
						new = gm_min_ind
						print "New negative component: "+new

			W(active_set==1) = W(active_set==1) +g_min*avec

			if positive:
				if any (W<0):
					min(W)
					glag=1

			lamb = lamb = g_min

		# Update Wights and lambdas

		lambdas[i] = lamb
		Ws[:,:,i] = W
		res = norm(Y-X*W, 'fro')^2

		if lamb==0 or (new and sum(active_set[:]) == maxcomps) or (res<noise):
			if verbose:
				print "End!"
			break


			active_set=1

		i=i+1
	# End of while loops

	# Final calculation of mus
	if flag==0:
		Ws = np.squeeze(Ws[:,:,1:length(lambdas)])
		w_dir = -Ws(:,i) - Ws(:,i-1)/(lambas[i] - lambdas[i-1])
		Aw = X*w_dir
		y_res = Y-X*(Ws(:,i-1) + w_dr(lambdas[i-1]))
		ld = np.roots([np.linalg.norm(Aw)^2, -2*(Aw*y_res), y_res.transpose()*y_res - noise])
		lam = ld(np.intersect1d(numpy.where(ld>lambdas[i]), numpy.where(ld<lambdas[i-1])))

		if len(lam)==0 or any(lam)<0 or any(np.isreal(lam)):
			lam=lambdas[i]
		W_lam = Ws(:,i-1) + w_dir*[lambdas[i-1] - lam[1]]
	else:
		W_lam=0
		Ws=0
		lambdas=0
		lam=0

# End of function

def calcAvec(new, dQ, W, lambda, active_set, M, positive):

	(r,c) = np.find(active_set)
	Mm = -M(r,r)

	Mm = (Mm + Mm.transpose())/2

	# Verify there is no numerical instability
	eigMm = np.linalg.eigen(Mm)
	if any(eigMm<0):
		min(eigMm)
		flag=1

	b=sign(W)
	if new:
		b(new) = np.sign(dQ(new))

	b= b[active_set==1]

	avec = np.divide(Mm,b)

	if positive:
		if new:
			in1 = sum(active_set[1:new])
			if avec(in1)<0:
				new
				flag=1


	one_vec = np.ones(size(W))
	dQa = np.zeros(size(W))

	for j in range(1, len(r)+1):
		dQa = dQa+avec[j]*M[:,r[j]]

	gamma_plus = np.divide((lamb - dQ),(one_vec + dQa))
	gamma_minus = np.divide((lamb + dQ), (one_vec - dQa))
