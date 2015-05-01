import math
import scipy as sp
import numpy as np
from Larse_regression_noise import larse_regression_noise

inf=float("inf")

def update_spatial_components(Y,Cin,f,A,P,C,R,coefs_est):
    (d,T) = np.shape(Y)
    dist=3
    d1 = R
    d2 = C
    Coor={}
    min_size = math.pow(8,2)
    max_size = math.pow(3,2)
    med_filt = np.array((1,2))
    med_filt.fill(3)
    thresh = 0.2

    temp = np.zeros((1,d1))
    temp2 = np.zeros((1,d2))
    for t in range(0,d1):
        temp[0,t] = t
    for t in range(0,d2):
        temp2[0,t] = t
    temp1 = np.zeros((0,d2))
    # print "Shape ofr kron", np.shape(np.kron(np.ones((d2,1)), temp.transpose()))
    Coor.update({'x' : np.kron(np.ones((d2,1)), temp.transpose())})
    Coor.update({'y' : np.kron(temp2.transpose(),np.ones((d2,1)))})
    print "shape",np.shape(Cin), np.shape(f)
    nr = len(Cin[:,0])

    # print "Number of neurons ", nr
    # Determine search area for each neuron
    cm = np.zeros((nr,2))
    vr = [0]*nr
    # Indicator for distance
    Ind = np.zeros((d,nr))
    print A, np.sum(A)
    # print "Shapes,", np.shape(Coor['x'].transpose()), np.shape(A)

    cm[:,0] = np.dot(Coor['x'].transpose(),A)/np.sum(A)
    cm[:,1] = np.dot(Coor['y'].transpose(),A)/np.sum(A)

    # for i in range(0,nr):
    #     cm[:,0] = Coor['x'].transpose()*A[:,1:nr]/np.sum(A[:,1:nr])
    #     cm[i,0] = np.divide(np.dot(Coor['x'].transpose(),A[:,i]), np.sum(A[:,i]))
    #     cm[i,1] = np.divide(np.dot(Coor['y'].transpose(),A[:,i]),np.sum(A[:,i]))
    temp=np.zeros((d,2))
    for i in range(0,nr):
        print "Neuron", i
        # print np.shape(Coor['x']), cm[i,0], np.shape(sp.sparse.spdiags(A[:,i],0,d,d)),nr
        # print Coor['x'][5,0], cm[i,0]
        for j in range(0,d):
            Coor['x'][j,0]-=cm[i,0]
            Coor['y'][j,0] -= cm[i,1]
        # print "Done with subtraction!"
        temp[:,0] = Coor['x'][:,0]
        temp[:,1] = Coor['y'][:,0]
        # print "here", np.shape(sp.sparse.spdiags(A[:,i],0,d,d)), temp
        temp1 = np.dot(sp.sparse.spdiags(A[:,i],0,d,d),temp)
        # print np.shape(temp1), np.shape(temp)
        vr[i] = np.dot(temp.transpose(), temp1)/np.sum(A[:,i])
        # print "here!"
        # print np.shape(vr[i])

        V,D = np.linalg.eig(vr[i])
        d11 = np.minimum(np.power(min_size,2), np.maximum(np.power(max_size,2), D[0,0]))
        d22 = np.minimum(np.power(min_size,2), np.maximum(np.power(max_size,2), D[1,1]))
        # print "here after d11,d22", np.shape(temp.transpose()), np.shape(np.power(np.dot(temp,D[:,0]),2))
        Ind[:,i] = np.sqrt(np.power(np.dot(temp,D[:,0]),2)/d11 + np.power(np.dot(temp,D[:,1]),2)/d22)<=dist
        print Ind

    # print np.shape(Cin), np.shape(f)
    # Cf = np.concatenate((Cin,f),axis=0)
    print np.shape(Cin),np.shape(f)
    Cf = np.concatenate((Cin,f ),axis=0)
    print np.shape(f),len(f[0,:]),np.shape(f[0,:])
    A_new = np.concatenate((np.zeros((d,nr)), np.zeros((d,len(f[0,:])))), axis=1)
    sA = np.zeros((d1,d2))

    # Estimate spatial components
    for px in range(0,d):
        if dist ==inf:
            print "Distance is infinity"
            a = larse_regression_noise(Y[px,:].transpose(), Cf.transpose(), 1, np.dot(np.power(coefs_est[px],2),T))
            A_new[px,:] = a.tranpose()
            sA[px] = np.sum(a)

        else:
            ind = np.nonzero(Ind[px,:])
            print "indexxxx", ind

            if len(ind)!=0:
                # ind = np.array(np.shape(f,1))
                ind = ind[0][0]
                print ind, nr+np.asarray(range(len(f)))
                ind2 = np.array((ind, nr+np.asarray(range(len(f)))))
                print ind,ind2, type(Cf),np.shape(Cf),np.shape(Y[px,:])
                a = larse_regression_noise(Y[px,:],Cf[ind2,:].transpose(), 1, np.dot(np.power(coefs_est[px],2),T))
                print a, px, ind2
                A_new[px,ind2] = a.tranpose()
                sA[px] = np.sum(a)

    A_new[np.isnan(A_new)]=0

    # Perform median filtering

    for i in range(0,nr):
        I_temp = sp.signal.medfilt(np.reshape(A_new[:,i],(d1,d2)),med_filt)
        print "Shapes!",np.shape(np.where(I_temp[:])), np.shape(np.where(A_new[:,i]))

        print "Shaaapes",np.where(I_temp[:]), np.where(A_new[:,i])
        acp = np.intersect1d(np.where(I_temp[:]),np.where(A_new[:,i]))
        print "ACP",np.shape(acp), acp
        print d, 1, acp, np.shape(acp)
        # A[:,i] = A[acp,i]
        # A[:,i]= sp.sparse.coo_matrix(A[acp,i])
        A_new[:,i] = sp.sparse.coo_matrix((A_new[acp,i],(acp,1)),shape=(d,1))

    A_new = sp.sparse.spdiags(A_new)
    Ath=A

    for i in range(0, nr):
        Ath[np.array(filter(lambda x: x>=thresh*np.max(Ath[:,i]), Ath[:,i])),i] = 0
        # Ath = Ath[:,i](filter (lambda x:))
        Bw = sp.ndimage.measurements.label(np.reshape(Ath[:,i],(d1,d2)))
        ml = np.max(Bw[:])
        ln = np.zeros((ml,1))

        for j in range(0, ml):
            ln[j] = np.size(np.where(Bw==j))
        x,ind = np.max(ln)
        Ath[Bw[:]!=ind,i] = 0

    A_new = Ath

    print "Updated spatial components"

    ff =np.where(np.sum(A)==0)
    print ff, np.shape(ff), np.shape(A_new), np.shape(Cin)
    if ff is not None:
        nr = nr - np.size(ff)
        A[:,ff] = np.zeros([])
        Cin[ff,:] = np.zeros([])

    Y_res = Y - np.dot(A[:,0:nr],Cin[0:nr,:])
    A_bas = np.maximum(np.divide(np.dot(Y_res,f.transpose()), np.power(np.linalg.norm(f),2)), 0)
    b = A_bas
    A_new = A_new[:,0:nr]

    return A_new, b


