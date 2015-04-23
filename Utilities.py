
def get_window_coordinates(i,wx,size):
	upper = i+(wx/2)
	lower = i-(wx/2)


	if upper>size or lower<0:
		upper=i
		for x in range(0, wx/2):
			if not upper>size:
				upper = i+x
		lower=i
		for y in range(-wx/2,0):
			if not lower<0:
				lower = i+x

	xv = np.linspace(upper,lower,1)
	return xv

def get_coordinates(d):
	if d<C:
		x=d
		y=0
	else:
		x = (d)%C
		y=(d+1)/C
	# z = y%P

	return x,y
	# return x,y,z

# def find_pixel(x,y,z):
# 	return x*R+y*C+z
def find_pixel(x,y):
	return x*C+y