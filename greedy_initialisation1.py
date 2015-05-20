import numpy as np
import scipy
from PIL import Image

def greedyROI(Y,K,K_size, K_std,C,R,T,P):
    (R,C,T) = np.shape(Y)
    nIter=1
    med = np.median(Y,axis=2)
    # Subtract median from data
    print "Calculating the median"
    for i in range(0,T):
        Y[:,:,i] = Y[:,:,i] - med


    K_size=np.array([1,2])
    K_size.fill(7)

    K_std = np.array([1,2])
    K_std.fill(2)

    basis = np.zeros((R,C,K),dtype='float64')
    trace = np.zeros((T,K),dtype='float64')
    center = np.zeros((K,2),dtype='float64')
    res = np.zeros((R,C,T))

    # Calculate Margin
    K_Half = np.floor(K_size/2)
    K_size = 2*K_Half + 1
    # K_size=1
    rho = np.zeros((R,C,T))
    print "Scanning the whole image", K_size,K_std
    # rho = scipy.ndimage.filters.gaussian_filter(Y,sigma=K_std)
    nan1=[]
    for i in range(0,T):
        print i
        rho[:,:,i] = scipy.ndimage.filters.gaussian_filter(Y[:,:,i],sigma=K_std)
        # rho[:,:,i] = scipy.ndimage.filters.gaussian_filter(Y[:,:,i],sigma=K_std,truncate=K_size)

    print "Calculating the explained variance"
    v = np.sum(np.power(rho,2), axis=2)
    ind1=[]
    ind2=[]
    for k in range(0,K):
    # Finding center with maximum variance
        ind = np.argmax(v)
        i_max,j_max = np.unravel_index(ind,np.shape(Y[:,:,0]))
        center[k,0]=i_max
        center[k,1]=j_max

        iSig = [np.maximum(i_max - K_Half[0],0), np.minimum(i_max + K_Half[0],R-1)]
        iSigLen = iSig[1] - iSig[0]
        jSig = [np.maximum(j_max - K_Half[1],0) , np.minimum(j_max + K_Half[1], C-1)]
        jSigLen = jSig[1] - jSig[0]


        # Fine tuning the shape
        print "Fine tuning the shape", iSig[0],iSig[1],jSig[0],jSig[1]
        dataTemp = Y[iSig[0]:iSig[1], jSig[0]:jSig[1], :]
        traceTemp = rho[i_max, j_max,:]
        coef, score = finetune2d(dataTemp,traceTemp,nIter,T)
        score.flatten()

        (M,N) = np.shape(coef)
        dataSig = np.empty((M,N,T))
        for i in range(0,T):
            for j in range(0,M):
                dataSig[j,:,i] = np.multiply(coef[j,:],score[j])
        basis[iSig[0]:iSig[1], jSig[0]:jSig[1],k] = coef
        trace[:,k] = score.transpose()

        Y[iSig[0]:iSig[1], jSig[0]:jSig[1],:] = Y[iSig[0]: iSig[1], jSig[0]:jSig[1], :] - dataSig
        print "Found",k,"number of neurons out of ",K

        if k<K:
            iMod = [np.maximum(i_max - 2*K_Half[0],0), np.minimum(i_max + 2*K_Half[1], R-1)]
            # Patches to modify
            iModLen = iMod[1] -iMod[0]
            jMod = [np.maximum(j_max - 2*K_Half[0],0), np.minimum(j_max + 2*K_Half[1], C-1)]
            # Patches to modify
            jModLen = jMod[1] -jMod[0]

            iLag = iSig - iMod[0]
            jLag = jSig - jMod[0]
            dataTemp = np.zeros((iModLen,jModLen))
            dataTemp[iLag[0]: iLag[1], jLag[0]:jLag[1]] = np.reshape(coef,[iSigLen, jSigLen])
            dataTemp = scipy.ndimage.filters.gaussian_filter(dataTemp,K_std)x

            (M,N) = np.shape(dataTemp)
            rhoTemp = np.empty((M,N,T))
            for i in range(0,T):
                for j in range(0,M):
                    rhoTemp[j,:,i] = np.multiply(dataTemp[j,:],score[i])
            rhoTemp = rho[iMod[0] : iMod[1], jMod[0] : jMod[1],:] - rhoTemp
            rho[iMod[0]:iMod[1], jMod[0] : jMod[1],:] = rhoTemp
            v[iMod[0] :iMod[1], jMod[0] :jMod[1]]=np.sum(np.power(rhoTemp,2),axis=2)

    return  basis,trace,center,Y


def finetune2d(data, trace, nIter,T):
    for ite in range(0, nIter):

        a = np.sum(np.power(trace,2))

        (M,N,T) = np.shape(data)
        trace = np.reshape(trace,(1,1,T))

        data = np.multiply(data,trace)
        b = np.sum(data,axis=2)
        basis = np.maximum(b/a,0)
        basisNorm = np.linalg.norm(basis.flatten())
        if basisNorm>0:
            basis = basis/basisNorm

        # Updating trace
        temp = np.zeros(np.shape(data))
        for i in range(0,T):
            temp[:,:,i] = np.multiply(data[:,:,i],basis)

        trace = np.squeeze(np.sum(np.sum(trace,axis=0),axis=0))

    return basis, trace



