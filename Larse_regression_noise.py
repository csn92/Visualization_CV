import numpy as np
import scipy
inf=float("inf")

def calcAvec(new, dQ, W, lamb, active_set, M, positive):

        (r,c) = np.nonzero(active_set)
        Mm = -M[r-1,r-1]

        Mm = (Mm + Mm.transpose())/2

        print Mm
        if Mm:
        # Verify there is no numerical instability
            eigMm = np.linalg.eig(Mm)
            if np.any(eigMm<0):
                np.min(eigMm)
                flag=1

            b=np.sign(W)
            if len(new)!=0:
                b[new] = np.sign(dQ(new))

            b= b[active_set==1]

            avec = np.divide(Mm,b)

            if positive:
                if new:
                    in1 = np.sum(active_set[0:new])
                    if avec(in1)<0:
                        new
                        flag=1


            one_vec = np.ones(np.shape(W))
            dQa = np.zeros(np.shape(W))

            for j in range(1, len(r)+1):
                dQa = dQa+np.dot(avec[j],M[:,r[j]])

            gamma_plus = np.divide((lamb - dQ),(one_vec + dQa))
            gamma_minus = np.divide((lamb + dQ), (one_vec - dQa))


def larse_regression_noise(Y,X,positive, noise):
        k=1
        print np.shape(Y[:])
        T = len(Y[:])

        a,N = np.shape(X)
        maxcomps = N;

        W = np.zeros((N,k))
        active_set = np.zeros((N,k))
        visited_set = np.zeros((N,k))

        Ws = np.zeros((N,k,maxcomps))

        r = np.dot(X.transpose(),Y[:])
        M = np.dot(-X.transpose(),X)

        i=1
        flag=0
        print np.shape(r), np.shape(M)
        lambdas=[]
        while 1:
            if flag==1:
                W_lam=0
                break

            if i>1 and (len(new)==0 and visited_set(new)==0):
                visited_set[new] = 1

            # Compute Gradient
            dQ = r+np.multiply(M,W)

            # Compute new W
            if i==1:
                if positive:
                    dQa = dQ
                else:
                    dQa = np.abs(dQ)
                    print "shape of dQa", np.shape(dQa)
                (m,n) = np.shape(dQa)
                dQa = np.reshape(dQa,(m*n))
                print "Shape",dQa
                (m,n) = np.shape(dQa)
                dQa = np.reshape(dQa,(m*n))
                print "Shaping",dQa, np.shape(dQa)
                lamb = np.zeros(len(dQa))
                ind = np.zeros(len(dQa))

                lamb = np.amax(dQa,axis=0)
                new = np.argmax(dQa, axis=0)
                # for i in range(0,len(dQa)):
                #     lamb[i] = np.max(dQa[i][:])
                #     ind[i] = np.argmax(dQa[i][:])
                print lamb, new
                if lamb.any()<0:
                    print "All negative directions!"
                    break
            else:
                # Calculate vector to travel along

                avec, gamma_plus, gamma_minus = calcAvec(new, dQ, W, lamb, active_set, M, positive)

                if new==0:
                    if dropped_sign==1:
                        gamma_plus[dropped] = inf
                    else:
                        gamma_minus[dropped] = inf

                gamma_plus[active_set == 1 ] = inf
                gamma_plus[gamma_plus <=0] = inf
                gamma_plus[gamma_plus > lamb] = inf
                (gp_min, gp_min,ind) = np.min(gamma_plus[:]);

                if positive==1:
                    gm_min = inf
                else:
                    gamma_minus[active_set==1] = inf
                    gamma_minus[gamma_minus > lamb] = inf
                    gamma_minus[gamma_minus <=0] = inf
                    (gm_min, gm_min_ind) = np.min(gamma_minus[:])

                [g_min, which1] = np.min([gp_min, gm_min]);
                if g_min == inf:
                    g_min = lamb

                # LARS check

                gamma_zero = -W[active_set == 1] /avec
                gamma_zero_full = np.zeros((N,k))
                gamma_zero_full[active_set==1] = gamma_zero
                gamma_zero_full[gamma_zero_full<=0] = inf
                gz_min, gz_min_ind = np.min(gamma_zero_full[:])

                if gz_min<g_min:
                    if verbose==1:
                        print "Dropping active weight:",gz_min_ind

                    active_set[gz_min_ind] = 0
                    dropped = gz_min_ind
                    dropped_sign = sign(W[dropped]);
                    W[gz_min_ind] =0

                    avec = avec[gamma_zero!=gz_min]
                    g_min = gz_min_ind
                    new=0

                elif g_min<lamb:
                    if which==1:
                        new = gp_min_ind
                        if verbose==1:
                            print "New positive component: ",new
                    else:
                        new = gm_min_ind
                        print "New negative component: ",new

                W[active_set==1] = W[active_set==1] +np.dot(g_min,avec)

                if positive==1:
                    if np.any(W<0):
                        W=np.min(W)
                        flag=1

                lamb = lamb - g_min

            # Update Wights and lambdas
            print lamb,i
            lambdas.append(lamb)
            Ws[:,:,i-1] = W
            res = np.power(np.linalg.norm(Y-np.dot(X,W), 'fro'),2)
            if (lamb==0).any() or ((len(new)==0 and np.sum(active_set[:]) == maxcomps) or (res<noise).any()):
                if verbose==1:
                    print "End!"
                break


                active_set[new]=1

            i=i+1
        # End of while loops

        # Final calculation of mus
        if flag==0:
            if i>1:
                Ws = np.squeeze(Ws[:,:,1:len(lambdas)])
                w_dir = -Ws[:,i] - Ws[:,i-1]/(lambas[i] - lambdas[i-1])
                Aw = np.dot(X,w_dir)
                y_res = Y-np.dot(X,(Ws[:,i-1]) + w_dr(lambdas[i-1]))
                ld = np.roots([np.linalg.norm(Aw)^2, -2*(np.dot(Aw,y_res)), np.dot(y_res.transpose(),y_res) - noise])
                lam = ld(np.intersect1d(np.where(ld>lambdas[i]), np.where(ld<lambdas[i-1])))

                if len(lam)==0 or (np.any(lam)<0 or np.any(np.isreal(lam))):
                    lam=lambdas[i]
                W_lam = Ws[:,i-1] + np.dot(w_dir,(lambdas[i-1]- lam[0]))

        else:
            W_lam=0
            Ws=0
            lambdas=0
            lam=0

    # End of function




