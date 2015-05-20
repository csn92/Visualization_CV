import numpy as np
import scipy

from cvxopt.modeling import op
import cvxopt

inf=float("inf")

def calcAvec(new, dQ, W, lamb, active_set, M, positive):
    # print "Active set  sv", active_set, np.nonzero(active_set)
    # print "Active ",active_set, np.nonzero(active_set)
    r,c = np.nonzero(active_set)
    # print "r,c ", r,c
  
    # print "print r, c",M, active_set, r,c
    Mm = -M[r,r]
    # print Mm
    Mm = (Mm + Mm.transpose())/2

    # print Mm
    # if Mm:
    # Verify there is no numerical instability
        # eigMm = np.linalg.eig(Mm)
        # if np.any(eigMm<0):
        #     np.min(eigMm)
        #     flag=1

    b=np.sign(W)
    # print "b, dQ, new ", b, dQ[new], dQ, new
    (m,n) = np.shape(dQ)
    dQ = np.reshape(dQ,(m*n))

    if new!=-1:
        b[new] = np.sign(dQ[new])

    dQ = np.reshape(dQ,(m,n))
    b= b[active_set==1]

    avec = np.divide(Mm,b)

    if positive:
        # print "avec positive",new
        if new!=-1:

            in1 = np.sum(active_set[0:new])
            # print "Active set ", in1, new, avec
            if avec[in1]<0:
                flag=1


    one_vec = np.ones(np.shape(W))
    dQa = np.zeros(np.shape(W))
    # print "afccghv jh ", dQa, one_vec
    # print "AVENCCCC",avec, r, len(r),M
    for j in range(0, len(r)):
        # print "svfv"
        (a,) = np.shape(avec[j]*M[:,r[j]])

        dQa = dQa+ np.reshape(avec[j]*M[:,r[j]],(a,1))
        # print "cgfchrtc ",np.shape(avec[j]*M[:,r[j]]), np.shape(avec[j]), np.shape(M[:,r[j]])

    # print "grammar", lamb, one_vec,dQ, dQa
    gamma_plus = np.divide((lamb - dQ),(one_vec + dQa))
    gamma_minus = np.divide((lamb + dQ), (one_vec - dQa))
    # print gamma_plus, gamma_minus, type(gamma_plus), type(gamma_minus), np.shape(gamma_plus), np.shape(gamma_minus)

    return avec, gamma_plus, gamma_minus


def larse_regression_noise(Y,X,positive, noise, ind):

    verbose = False
    k=1
    # print np.shape(Y[:])
    # print "larse ",np.shape(Y)
    T = len(Y)
    # print "Noise ", noise[0], noise[1], noise[3]
    noise =np.sum(noise)/T

    (a,N) = np.shape(X)
    # print "larse regression noie ",N,a
    maxcomps = N

    W = np.zeros((N,k))
    active_set = np.zeros((N,k))
    visited_set = np.zeros((N,k))

    Ws = np.zeros((N,k,maxcomps))
    # print "Shape of Y", np.shape(Y),np.shape(Y[:]),np.shape(X.transpose())
    r = np.dot(X.transpose(),Y[:])
    # print "Shape of rrrrrr", np.shape(r)
    M = np.dot(-X.transpose(),X)
    # print np.shape(M), np.shape(r)
    # new = np.zeros((1,1))
    i=1
    flag=0
    # print np.shape(r), np.shape(M)
    lambdas=[]

    while 1:
        if flag==1:
            W_lam=0
            break

        if i>1 and new!=0:
            # print np.shape(visited_set)
            if len(visited_set.shape)!=1:
                (a,b) = np.shape(visited_set)
                visited_set=np.reshape(visited_set,(a*b))
                if visited_set[new]==0:
                    visited_set[new] = 1
                visited_set=np.reshape(visited_set,(a,b))

            else:
                if visited_set[new]==0:
                    visited_set[new] = 1

        # Compute Gradient
        # print "gradient ", np.shape(M), np.shape(W), np.shape(r)
        if len(r.shape)==1:
            (s,) = np.shape(r)
            r = np.reshape(r,(s,1))
        dQ = r+ np.dot(M,W)

        # print "dQ shaoeeeee ", np.shape(dQ)

        # Compute new W
        if i==1:
            if positive:
                dQa = dQ
            else:
                dQa = np.abs(dQ)
                # print "shape of dQa", np.shape(dQa)
            (m,n) = np.shape(dQa)
            dQa = np.reshape(dQa,(m*n))
            # print "Shape",dQa, np.shape(dQa),dQ

            lamb = np.amax(dQa)
            new = np.argmax(dQa)
            # print "dqa, lamb", dQa, lamb,new
           
            dQa = np.reshape(dQa,(m,n))
           

            # print lamb, new
            if lamb <0:
                print "All negative directions!"
                break
        else:
            # print "visited set",visited_set, new, i
            # Calculate vector to travel along
            # print "else part i ", i, visited_set, active_set
        
            # print type(new), type(dQ), type(W), type(lamb), type(active_set), type(M), type(positive)
            if i==1:
                active_set = visited_set

            avec, gamma_plus, gamma_minus = calcAvec(new, dQ, W, lamb, active_set, M, positive)
            # print "Adcsvwvsvs"
            if new==-1:

                if dropped_sign==1:
                    # print "HSHDH"
                    gamma_plus[dropped] = inf
                else:
                    # print "gbtbb"
                    gamma_minus[dropped] = inf
            # print "I am hreew"
            active_set = visited_set
            # print "Ac tive set", active_set
            gamma_plus[active_set == 1 ] = inf
            gamma_plus[gamma_plus <=0] = inf
            gamma_plus[gamma_plus > lamb] = inf
            gp_min = np.min(gamma_plus)
            gp_min_ind = np.argmin(gamma_plus)
            # print "Gamma plus shape", gp_min,gp_min_ind,gamma_plus 


            if positive:
                # print "positive here"
                gm_min = inf
            else:
                gamma_minus[active_set==1] = inf
                gamma_minus[gamma_minus > lamb] = inf
                gamma_minus[gamma_minus <=0] = inf
                # print "Gamma minus shape", np.shape(gamma_minus), np.shape(np.min(gamma_minus)), np.shape(np.argmin(gamma_minus))
                gm_min = np.min(gamma_minus)
                gm_min_ind = np.argmin(gamma_minus)
                # print "gamma minus ", gm_min_ind,gm_min,gamma_minus


            g_min = np.min([gp_min, gm_min])
            which1 = np.argmin([gp_min,gm_min])
            # print "G min, G_plus, G_min", gp_min, gm_min, g_min, which1

            if g_min == inf:
                g_min = lamb

            # LARS check
            (a,) = np.shape(avec)
            # print "avec,a ",avec,a
            avec = np.reshape(avec, (1,a))
            # print "Avec later", avec, -W[active_set == 1]
            # (b,) = np.shape(-W[active_set == 1])
            # print "av,b ", -W[active_set == 1], b
            # temp = np.reshape(-W[active_set == 1], (b,1))

            gamma_zero = np.divide(-W[active_set == 1],avec)
            # gamma_zero = np.zeros(np.shape(-W[active_set == 1]))
            # print np.shape(gamma_zero), np.shape(-W[active_set == 1])
            # for i in range(0, len(avec)):
            #     gamma_zero[i,:] = -W[active_set == 1][i,:]/avec[i]
            gamma_zero_full = np.zeros((N,k))
            # print W,-W[active_set == 1],np.shape(-W[active_set == 1]), np.shape(avec.transpose()), avec.transpose()
            # print "36 hours ",np.shape(gamma_zero), gamma_zero_full,active_set==1, gamma_zero_full[active_set==1].shape, gamma_zero,-W[active_set == 1], np.shape(-W[active_set == 1]), np.shape(avec), avec
            # for i  in range(0, len(active_set)):
            #     if active_set==1:
            #         gamma_zero_full[i,1] = gamma_zero[i,1]
            # print np.shape(gamma_zero_full[active_set==1]), np.shape(np.ravel(gamma_zero_full))
            # print gamma_zero,  gamma_zero_full[active_set==1], len(gamma_zero_full[active_set==1].shape), np.ravel(gamma_zero)
            # (b,) = np.shape(gamma_zero_full[active_set==1])
            # print b
            # a = np.shape(gamma_zero_full)
            # gamma_zero_full = np.ravel(gamma_zero_full[active_set==1]) 
            # print "gamma zero full ", gamma_zero_full, a
            # gamma_zero_full = gamma_zero
            # gamma_zero_full = np.reshape(gamma_zero_full,a)
            # gamma_zero_full[active_set==1] = np.ravel(gamma_zero)
            if len(gamma_zero_full[active_set==1].shape)==1:
                gamma_zero_full[active_set==1] = np.ravel(gamma_zero)
            else:
                 gamma_zero_full[active_set==1] = gamma_zero
            gamma_zero_full[gamma_zero_full<=0] = inf
            gz_min = np.min(gamma_zero_full)
            gz_min_ind = np.argmin(gamma_zero_full)
            # print "Shape of Gamma Zero Full", gz_min, gz_min_ind, gamma_zero_full

            # print gz_min,g_min,lamb,which1
            if gz_min<g_min:
                if verbose:
                    print "Dropping active weight:",gz_min_ind

                active_set[gz_min_ind] = 0
                dropped = gz_min_ind
                dropped_sign = np.sign(W[dropped]);
                W[gz_min_ind] =0

                avec = avec[gamma_zero!=gz_min]
                g_min = gz_min
                new=-1

            elif g_min<lamb:
                if which1==0:
                    new = gp_min_ind
                    # print "Here g_min ", new
                    if verbose:
                        print "New positive component: ",new
                else:
                    new = gm_min_ind
                    # print "Here gm_min_ind ", new
                    print "New negative component: ",new

            # print np.shape(g_min *avec), np.shape(W[active_set==1])
            if len((g_min *avec))==1:
                # print 'W acrtidb ', W
                W[active_set==1] = W[active_set==1] + np.ravel((g_min *avec))
                # print "gretg after ",W
            else:
                # print 'W acrtidb ', W
                W[active_set==1] = W[active_set==1] + (g_min *avec)
                # print "gretg after ",W

            if positive:
                if np.any(W<0):
                    # print "Shape of W", np.shape(W), np.min(W)
                 
                    flag=1

            lamb = lamb - g_min

        # Update Wights and lambdas
        # print lamb,i
        lambdas.append(lamb)
        Ws[:,:,i-1] = W

        # print "RES Value ",Y,np.shape(W),W, np.shape(X),np.dot(X,W)
        res = np.power(np.linalg.norm(Y-np.dot(X,W), 'fro'),2)

        # print "Prinitng active set ", np.shape(active_set), active_set
        if len(active_set.shape)!=1:
            (a,b) = np.shape(active_set)
            active_set = np.reshape(active_set,a*b)
            # print "New ", new, visited_set
            active_set[new]=1
            active_set = np.reshape(active_set,(a,b))
        else:
            active_set[new]=1
        # print np.shape(lamb), np.shape(res), np.shape(noise),np.shape(new), res

        if (lamb==0 or ((new==-1 and np.sum(active_set) == maxcomps) or (res<noise) )):
            if verbose:
                print "End!"
            break

        i=i+1
        # print "i value ",i, active_set, np.nonzero(active_set)
    # End of while loops

    # Final calculation of mus
    if flag==0:
        if i>1:
            
            # print "Wssssss", Ws, np.shape(Ws), len(lambdas), np.shape(Ws[:,:,0:len(lambdas)])
            Ws = np.squeeze(Ws[:,:,0:len(lambdas)])

            # print "Shape of Ws ",len(lambdas), lambdas, Ws, np.shape(Ws), np.shape(lambdas)

            w_dir = -Ws[:,i-1] - Ws[:,i-2]/(lambdas[i-1] - lambdas[i-2])
            # w_dir = -Ws[i-1] - Ws[i-2]/(lambdas[i-1] - lambdas[i-2])
            Aw = np.dot(X,w_dir)

            # print "shapes of w_dir ", np.shape(w_dir), np.shape(X*Ws[:,i-2]),lambdas[i-2]
            # y_res = Y-np.dot(X,Ws[i-2]+ w_dir[lambdas[i-2]])
            y_res = Y-np.dot(X,Ws[:,i-2] + w_dir*lambdas[i-2])
            ld = np.roots([np.power(np.linalg.norm(Aw),2), -2*(np.dot(Aw.transpose(),y_res)), np.dot(y_res.transpose(),y_res) - noise])
            
            # print "Intesect 1 D  ", np.where(ld>lambdas[i-1])[0], np.where(ld<lambdas[i-2])[0]
            lam = ld[np.intersect1d(np.where(ld>lambdas[i-1])[0], np.where(ld<lambdas[i-2])[0] )]

            # print "lam before",lam
            if len(lam)==0 or (np.any(lam)<0 or np.any(np.iscomplex(lam))):
                lam=lambdas[i-1]

            # print "Shape of lam ",lam, 
            W_lam = Ws[:,i-1] + np.dot(w_dir,(lambdas[i-2]- lam))
            # W_lam = Ws[i-1] + np.dot(w_dir,(lambdas[i-2]- lam[0]))

        else:
            # W_lam = cvxopt.modeling.variable(len(X[0,:]))
            # c1 = (W_lam>=0)
            # c2 = [(Y-X*W_lam)<=np.sqrt(noise)]
            # lp1 = op(np.sum(W_lam),[c1,c2])
            # lp1.solve()
            # print lp1.status
            # lam = 10
            print "ind shape",np.shape(ind)
            (size,) = np.shape(ind)

            W_lam = np.zeros(size+1)
   #          W_lam = [  0.00000000e+00 ,  0.00000000e+00 , -1.33793646e+22 ,  0.00000000e+00,
   # 0.00000000e+00,  -1.31532477e+22,   0.00000000e+00 ,  0.00000000e+00]
    else:
        W_lam=0
        Ws=0
        lambdas=0
        lam=0


    return W_lam
# End of function




