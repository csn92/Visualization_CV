import scipy as sp
import numpy as np

def make_G_matrix(T, g):
	# if len(g)==1 and g<0:
	# 	g=0
	print np.shape(g)
	# G = sp.sparse.spdiags(np.dot(np.ones((T,1)), np.flipud(g)).T, len(g),T,T)
	return np.random.rand(T)
	# sl = np.cumsum(varargin.keys()[0][:])
	# for i in range(0, len(sl)-1):
	# 	G[sl[i]+1,sl[i+1]] = 0

def merge_ROIs(Y_res, A_, b, C_,f,P):
	if 'merge_thr' in P:
		thr = 0.85
	else:
		thr = P['merge_thr']

	if 'max_merg' in P:
		mx = 50
	else:
		mx = P['max_merg']

	nr = np.size(A,axis=1)

	(d,T) = np.size(Y_res)
	C_corr = np.corr(np.full(C[0:nr,:]).tranpose())
	FF1 = np.triu(C_corr)>=thr

def plain_foopsi(H,D,I_est,eps):
	ln = len(H)
	step_back_frac = 0.5
	iter = 0
	if nargin == 2:
	    I_est = 1e-3*np.ones(ln,1)
	    eps = 1;
	
	Zin = np.reshape(I_est,[6])

	if nargin == 3:
	    eps = 1;

	while eps>1e-5:
	    n = np.dot(D,Zin)
	    nnd = 10;
	    E = np.linalg.norm(Zin-H)^2 - np.dot(eps,np.sum(log(D*Zin)))
	    grad = 2*(Zin-H) - np.dot(np.dot(eps,D.transpose()),(np.power(n,-1)))
	    Hs = 2*scipy.sparse.speye(ln) + np.dot(np.dot(np.dot(eps,D.transpose()),scipy.sparse.spdiags(np.power(n,-2),0,ln,ln)),D)        
	    while nnd/2>1:
	        iter = iter + 1;
	        Z_dir = np.divide(-Hs,grad)
	        hit = np.divide(-n,(np.dot(D,Z_dir)))
	        if np.all(hit<0):
	            s = 1
	        else:
	            s = np.min(1,.9*np.min(hit(hit>=0)))
	        end
	        E_new = E
	        s = s/step_back_frac
	        x_dg = np.dot(grad.transpose(),Z_dir)
	        while E_new > E + 0.25*s*x_dg:
	            s=s*step_back_frac
	            Z_new = Zin + s*Z_dir
	            n = D*Zin
	            E_new = norm(Z_new-H)^2 - eps*sum(log(D*Z_new));
	        end
	        # E = E_new;
	        Zin = Zin + s*Z_dir
	        nnd = -x_dg
	        E = norm(Zin-H)^2 - eps*sum(log(D*Zin))
	        n = D*Zin
	        grad = 2*(Zin-H) - np.dot(np.dot(eps,D.transpose()),(np.power(n,-1)))
	        Hs = 2*scipy.sparse.speye(ln) + np.dot(np.dot(np.dot(eps,D.transpose()),scipy.sparse.spdiags(np.power(n,-2),0,ln,ln)),D)
	        # %disp(nnd)
	    eps = eps/10

	# %fprintf('Interior point method converged after %i iterations \n',iter);
	ip_it = iter
