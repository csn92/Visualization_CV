import scipy as sp
import numpy as np


def make_G_matrix(T, g):

	G = sp.sparse.spdiags(np.dot(np.ones((1,T)), np.flipud(g).transpose()), len(g),T,T)

	return G

def plain_foopsi(H,D):
	ln = len(H)
	step_back_frac = 0.5
	iter = 0
	I_est = 1e-3*np.ones((ln,1))
	eps = 1
	
	Zin = np.reshape(I_est,(ln,1))

	while eps>1e-5:
		n = D*Zin
		nnd = 10
		E = np.power(np.linalg.norm(Zin-H),2) - np.dot(eps,np.sum(np.log(D*Zin)))

		temp = np.ones(np.shape(n), dtype='float64')

		temp = np.reshape(np.power(n[:,0],-1, dtype='float64'), (T,1))

		grad = 2*(Zin-H) - eps*D.transpose()*temp
		Hs = 2*sp.sparse.eye(ln) + eps*np.dot(np.dot(D.transpose(),sp.sparse.spdiags(np.power(n,-2, dtype='float64').transpose(),0,ln,ln)),D)        
		while nnd/2>1:
			iter = iter + 1
			Z_dir, Zid, Zrank,Zs = np.linalg.lstsq(Hs, grad)
			temp = np.random.rand(200,1)
			hit = np.true_divide(-n,temp)
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
			    E_new = norm(Z_new-H)^2 - eps*np.sum(log(D*Z_new));
		
			Zin = Zin + s*Z_dir
			nnd = -x_dg
			E = norm(Zin-H)^2 - eps*np.sum(log(D*Zin))
			n = D*Zin
			grad = 2*(Zin-H) - np.dot(np.dot(eps,D.transpose()),(np.power(n,-1)))
			Hs = 2*scipy.sparse.speye(ln) + np.dot(np.dot(np.dot(eps,D.transpose()),scipy.sparse.spdiags(np.power(n,-2),0,ln,ln)),D)
			eps = eps/10
			
	ip_it = iter
