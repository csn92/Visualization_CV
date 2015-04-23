Inf=999999
import Utilities


def update_spatial_components(Y,C,f,A,P):

	if 'dist' in P:
		dist=3
	else:
		dist = P[dist]
	if 'd1' in P:
		d1 = math.sqrt(d)
	else:
		d2 = P[d1]
	if 'd2' in P:
		d2 = math.sqrt(d)
	else:
		d2 = P[d2]
	if 'min_size' in P:
		
	(d,T) = np.shape(Y)
	nr = np.shape(C)

	# Determine search area for each neuron

	if not dist==Inf:
		cm = np.zeros(nr,3)
		Ind = np.zeros(nr,1)
