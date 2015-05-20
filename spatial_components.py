import math
import scipy as sp
import numpy as np
# import spams_python 
import spams
# import cvxopt as cvx
from Larse_regression_noise import larse_regression_noise

inf=float("inf")

def update_spatial_components(Y,Cin,bin,f,A,P,C,R,coefs_est):
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
    Coor.update({'y' : np.kron(temp2.transpose(),np.ones((d1,1)))})
    print "shape",np.shape(Cin), np.shape(f)
    nr = len(Cin[:,0])

    if not dist==inf:
        cm = np.zeros((nr,2))
        vr = [0]*nr

        # Indicator for distance
        Ind = np.zeros((d,nr))

        cm[:,0] = np.dot(Coor['x'].transpose(),A)/np.sum(A)
        cm[:,1] = np.dot(Coor['y'].transpose(),A)/np.sum(A)

        temp=np.zeros((d,2))
        for i in range(0,nr):
            for j in range(0,d):
                Coor['x'][j,:]-=cm[i,0]
                # print j
                Coor['y'][j,:]-= cm[i,1]
            # print "Done with subtraction!"
            temp[:,0] = Coor['x'][:,0]
            temp[:,1] = Coor['y'][:,0]
            temp1= sp.sparse.spdiags(A[:,i],0,d,d)*temp
            vr[i] = np.dot(temp.transpose(),temp1)/np.sum(A[:,i])

        
            V,D = np.linalg.eig(vr[i])
            d11 = np.minimum(np.power(min_size,2), np.maximum(np.power(max_size,2), D[0,0]))
            d22 = np.minimum(np.power(min_size,2), np.maximum(np.power(max_size,2), D[1,1]))
            Ind[:,i] = np.sqrt(np.power(np.dot(temp,D[:,0]),2)/d11 + np.power(np.dot(temp,D[:,1]),2)/d22)<=dist

  
    Cf = np.concatenate((Cin,f ),axis=0)
    A_new = np.concatenate((np.zeros((d,nr)), np.zeros((d,len(f[0,:])))), axis=1)
    sA = np.zeros(d)

    # Estimate spatial components
    for px in range(0,d):
        if dist ==inf:
            print "Distance is infinity"
            a = larse_regression_noise(Y[px,:].transpose(), Cf.transpose(), True, np.dot(np.power(coefs_est[px],2),T), ind)
            A_new[px,:] = a.transpose()
            sA[px] = np.sum(a)

        else:
            ind = np.nonzero(Ind[px,:])
     
            if len(ind[0])!=0:
   
                ind = ind[0]
                s1,s2 = np.shape(f)
                ind2 = np.concatenate((ind, nr+np.arange(0,s1,1)))
                mul = np.dot(bin,f)
                print px
                a = larse_regression_noise(Y[px,:].transpose(),Cf[ind2,:].transpose(), True, np.dot(np.power(coefs_est[px],2),T), ind)
                A_new[px,ind2] = a
                sA[px] = np.sum(a)

    # np.savez('variables1',A_new,sA)
    # nr=50
    # with np.load('variables1.npz') as data:
    #     A_new = data['arr_0']
    #     sA = data['arr_1']

    print "A new ",np.shape(A_new)
    print "sA ", np.shape(sA)
    A_new[np.isnan(A_new)]=0

    # Perform median filtering

    for i in range(0,nr):
        I_temp = sp.signal.medfilt(np.reshape(A_new[:,i],(d1,d2)),med_filt)
        I_temp = np.reshape(I_temp, d1*d2)
        acp = np.intersect1d(np.where(I_temp)[0],np.where(A_new[:,i])[0])
        acp = np.unravel_index(acp,(d1,d2))
        j = np.zeros(np.shape(acp), dtype='int64')
        temp1 = np.reshape(A_new[:,i], (d1,d2))
        temp=[]
        for h in range(0, len(acp[0])):
        
            temp.append(temp1[acp[0][h],acp[1][h]])
  
        t = sp.sparse.coo_matrix((np.asarray(temp),(acp[0],acp[1])),shape=(d1,d2))

        A_new[:,i] = np.reshape(t.todense(), (d,))

    A_new = sp.sparse.coo_matrix(A_new)
    Ath=A

    for i in range(0, nr):
        Ath[Ath[:,i] < thresh*np.max(Ath[:,i]) , i] = 0
        Bw = sp.ndimage.measurements.label(np.reshape(Ath[:,i],(d1,d2)))
        (p,l) = Bw
        ml = np.max(np.reshape(p,d1*d2))
        ln = np.zeros((ml,1))

        for j in range(0, ml):
            ln[j] = np.size(np.where(p==j))
        x= np.max(ln)
        ind = np.argmax(ln)
        Ath[np.reshape(p, d1*d2)!=ind,i] = 0

    A_new = Ath

    print "Updated spatial components"

    ff =np.where(np.sum(A)==0)
    if ff is not None:
        nr = nr - np.size(ff)
        A[:,ff] = np.zeros([])
        Cin[ff,:] = np.zeros([])

    Y_res = Y - np.dot(A[:,0:nr],Cin[0:nr,:])
    A_bas = np.maximum(np.divide(np.dot(Y_res,f.transpose()), np.power(np.linalg.norm(f),2)), 0)
    b = A_bas
    A_new = A_new[:,0:nr]

    return A_new, b


