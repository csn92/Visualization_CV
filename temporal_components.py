import numpy as np
import scipy as sp
from Utilities import make_G_matrix, plain_foopsi

def update_temporal_components(Y,A,b,Cin,fin,coefs):

	method = 'project'
	restimate_g = 1
	temporal_iter = 1

	(d,T) = np.shape(Y)
	flag_G = True

	g=coefs
	P={}
	print np.shape(g)
	if 'g' not in P:
		flag_G = False
		G = make_G_matrix(T,g)

	nr = 50

	b = np.concatenate((b, np.random.rand(d,1)), axis=1)
	A = np.concatenate((A,b),axis=1)
	Cin = np.concatenate((Cin,fin),axis=0)
	C = Cin
	Iter=1

	if method=='noise_constrained':
		Y_res =Y - np.dot(A,Cin)
		mc = np.minimum(d,15)
		if Ld is None:
			Ld = 10*np.ones(mc,nr)

	else:
		nA = np.sum(np.power(A,2), axis=1)
		YrA = np.dot(Y.transpose(),A) - np.dot(Cin.transpose(), np.dot(A.transpose(),A))
		if method=='constained_foopsi':
			p['gn'] = [0]*nr
			p['b'] = [0]*nr
			p['c1'] = [0]*nr
			p['neuron_sn'] = [0]*nr
			p['bas_nonneg'] = 0

	for it in range(0, Iter):
		perm = np.random.permutation(range(nr))
		for jj in range(0,nr):
			ii = perm[jj]

			if ii<=nr:
				print "sfe", flag_G
				if flag_G:
					G = make_G_matrix(T,g)
					print "G is ", G
				if method=='project':
					print np.shape(nA), np.shape(YrA), np.shape(Cin)
					YrA[:,ii] = YrA[:,ii] + np.dot(nA[ii],Cin[ii,:].transpose())
					cc = plain_foopsi(np.divide(YrA[:,ii],nA[ii]),G)
					C[ii,:] = np.full(cc.transpose())
					YrA[:,ii] = YrA[:,ii] - np.dot(nA[ii],C[ii,:].transpose())

			else:
				YrA[:,ii] = YrA[:,ii] + np.dot(nA[ii],Cin[ii,:]).transpose()
				cc = np.maximum(YrA[:,ii]/nA[ii],0)
				C[ii,:] = np.full(cc.tranpose())
				YrA[:,ii] = YrA[:,ii] - np.dot(nA[ii],C[ii,:]).tranpose()

				if jj%10==0:
					print jj,"out of total ", nr," temporal componenets updated!"

		if np.linalg.norm(Cin - C, 'fro')/np.linalg.norm(C,'fro') <=1e-3:
			print "Temporal compoenents not changing by much!"
			break
		else:
			Cin = C


		Y_res = Y - np.dot(A,C)

	f = C[nr+1:,:]
	C = C[0:nr,:]

	return C, f, Y_res, P, Ld

