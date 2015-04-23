import numpy as np

wx=1
wy=1
wz=1
x=5
y=5
z=5
A = np.random.rand(5,5,5)
bp = 3,2,3
cx = bp/x
cy = cx/y
cz = math.mod(cy,z)

print cx,cyz,cz

# Return coordinate matrices from two coordinate vectors
X,Y = np.meshgrid([1,2,3])

# points = np.append(A[0].reshape(-wx,wx), A[1].reshape(-wy,wy),A[2].reshape(-wz,wz))
x,y,z = np.meshgrid(A[0],A[1],A[2],indexing='ij')

assert np.all(x[:,0,0] == A[0])
assert np.all(y[0,:,0] == A[1])
assert np.all(z[0,0,:] == A[2])

positions = np.vstack([X.ravel(), Y.ravel()])

print points,A