from pulseq_helper import round_up_to_raster

class diff_params():
    """
    class to set and calculate diffusion gradient parameters assuming rectangular shape
    """
    def __init__(self, b_val = None, delta = None, spacing = None):
        self.b_val = b_val
        self.delta = delta
        self.spacing = spacing

    def set_param(self, b_val = None, delta = None, spacing = None):
        if b_val:
            self.b_val=b_val # [s/mm^2]
        if delta:
            self.delta=delta # [s]
        if spacing:
            self.spacing=spacing # [s]

    def calc_params(self, system):
        """ calculate parameter from other two given parameters
            it is assumed that the trapezoid is symmetric
            unit of system parameters is [Hz/m] - contain gamma/2*pi already
        """
        rise_time = round_up_to_raster(system.max_grad/system.max_slew, decimals=5) # maximum gradient strength and slewrate is always used
        if self.spacing is None:
            try:
                self.spacing = (self.b_val*1e6/(2*np.pi*system.max_grad)**2-rise_time**3/30+self.delta*rise_time**2/6)/self.delta**2+self.delta/3 # [s]
                self.spacing = round(self.spacing, 5) # round to gradient raster time
                print('Spacing of diffusion gradients: '+'{:.2f}'.format(self.spacing*1e3)+' ms')
            except:
                raise ValueError('Select 2 of 3 diffusion parameters')
        elif self.delta is None:
            try:
                # calculate roots of polynomial
                p = np.array([-((2*np.pi*system.max_grad)**2)/3, ((2*np.pi*system.max_grad)**2)*self.spacing, -1*(2*np.pi*system.max_grad)**2*rise_time**2/6, -1*self.b_val*1e6+(2*np.pi*system.max_grad)**2*rise_time**3/30])
                roots = np.roots(p)
                # check which root to choose by comparing with targeted b-value
                diff = self.b_val - 1e-6*(2*np.pi*system.max_grad)**2*(roots**2*(self.spacing-roots/3)+rise_time**3/30-roots*rise_time**2/6) 
                self.delta = roots[np.argmin(diff)]
                print('Duration of diffusion gradients: '+'{:.2f}'.format(self.delta*1e3)+' ms')
            except:
                raise ValueError('Select 2 of 3 diffusion parameters')
        elif self.b_val is None:
            try:
                self.b_val = 1e-6*(2*np.pi*system.max_grad)**2*(self.delta**2*(self.spacing-self.delta/3)+rise_time**3/30-self.delta*rise_time**2/6)
                print('b-value: '+'{:.2f}'.format(self.b_val))
            except:
                raise ValueError('Select 2 of 3 diffusion parameters')
        else:
            raise ValueError('Select only 2 of 3 diffusion parameters.')

def calc_bval(grad,spacing):
    """calculates the b-value of a gradient pair
    
    grad: one gradient of the equal gradient pair
    spacing: spacing of the gradients    
    """
    delta = grad.rise_time + grad.flat_time
    return 1e-6*(2*np.pi*grad.amplitude)**2*(delta**2*(spacing-delta/3)+grad.rise_time**3/30-delta*grad.rise_time**2/6)


""" Discoball
Generates a DISCS direction scheme with latitudes starting at the
northpole (DISCSo = "odd") or not at the northpole (DISCSe = 'even).

Input args: 
            N: Requested number of isotropically distributed directions
               in default mode.

            bNisNtheta: (optional) boolean for either
                        - False (default): DISCS scheme with N directions or
                        - True: DISCS scheme with n=N latitudes. The final
                        number of directions corresponds to the optimal
                        configuration with n latitudes.

Output args: s: N x 3 matrix with Cartesian direction coordinates on unit
                sphere where N is the number of directions.


Copyright (C) 2013 Ruediger Stirnberg

Original Matlab code transferred to Python by Marten Veldmann
"""

import numpy as np

def DISCSo(N, bNisNtheta=False):

    # geometry factor
    alpha = np.deg2rad(30)
    epsilon = 1/np.cos(alpha/2) - 1

    if bNisNtheta:
        n0 = N
        k_ = np.zeros(int(np.floor(n0*0.5)))
        n_ = n0
        for i,_ in enumerate(k_):
            if i+1==n0*0.5:
                k_[i] = round(0.5*np.pi/np.arcsin(np.sin(np.pi/n0/2)/np.sin((i+1)*np.pi/n0)))
            else:
                k_[i] = round(np.pi/np.arcsin(np.sin(np.pi/n0/2)/np.sin((i+1)*np.pi/n0)))
        N = int(sum(k_))+1
    elif N==1:
        n_ = 0
    else:
        n0 = np.pi / (2* np.arccos(-np.pi/(4*(1+epsilon)*(N-1)) + np.sqrt((np.pi/(4*(1+epsilon)*(N-1)))**2+1)))
        n0 = [int(np.floor(n0)), int(np.ceil(n0))]
        if 1 in n0:
            n0.remove(1)
    
        # init. final parameters
        m_ = np.inf             # best mean squared distance error    (1 x 1)
        n_ = None               # best number of latitudinal circles  (1 x 1)
        k_ = np.array([])       # best population of lat. circles     (n_ x 1)

        for n in n0:
            # n/2 (max. lat. circle in upper hemisphere)
            nhalf = int(np.floor(n*0.5))

            # mean n.n. distance
            d0 = 2*(1+epsilon)*np.sin(np.pi/2/n)
            # reasonable population range per latitude
            k = np.stack(nhalf * [np.arange(4*n)+1])
            # latitude index matrix of same size
            i = np.stack(k.shape[1] * [np.arange(nhalf)+1]).T
            # corresponding distances
            d = 2*np.sin(i*np.pi/n)*np.sin(np.pi/k)
            if not n%2: # equator            
                d[-1,:] = 2*np.sin(i[-1,:]*np.pi/n)*np.sin(np.pi/k[-1,:]/2)

            # corresponding squared errors
            e = (d0-d)**2
            ind = np.argsort(e, axis=1)
            e = np.sort(e, axis=1)
            for i in range(nhalf):
                k[i,:] = k[i,ind[i,:]]

            # optimal parameters for given n(l):
            E = e[:,0]
            K = k[:,0]
            N0 = 1+sum(K)

            # optimal population for given n -> best population for given N:
            i2=np.array([]); e2=np.array([]); k2=np.array([])
            Ndiff = N0-N
            if Ndiff!=0:
                for i in range(nhalf):
                    if Ndiff>0:
                        # indices available for population reduction
                        tmpj=np.where(k[i,:]<K[i])[0]
                    else:
                        # indices available for population increase
                        tmpj=np.where(k[i,:]>K[i])[0]

                    i2 = np.concatenate((i2,i*np.ones(tmpj.shape)), axis=0)
                    e2 = np.concatenate((e2,e[i,tmpj]), axis=0)
                    k2 = np.concatenate((k2,k[i,tmpj]), axis=0)       

                # change by Ndiff in total on those latitudes where error is least
                ind = np.argsort(e2)
                e2 = np.sort(e2)
                k2 = k2[ind]
                i2 = i2[ind]

                # update best parameters for given n and N
                E[np.int64(i2[:int(abs(Ndiff))])] = e2[:int(abs(Ndiff))]
                K[np.int64(i2[:int(abs(Ndiff))])] = k2[:int(abs(Ndiff))]

            # mean squared error
            m = sum(E*(K-1))/N

            # update final parameters for given N
            if m<m_:
                m_ = m
                n_ = n
                k_ = K.copy()
            
    # compute coordinates
    s = np.zeros([N,3]) # coordinate array
    s[0,:] = [0,0,1]
    o = 0
    for i in range(int(np.floor(n_*0.5))):
        j = (np.arange(k_[i])+1).T

        theta = (i+1)*np.pi/n_
        if i+1<n_*0.5:
            phi = (i+1+j)*2*np.pi/k_[i]
        else: # equator
            phi = (i+1+j)*np.pi/k_[i]

        st = np.sin(theta); ct = np.cos(theta)
        sp = np.sin(phi);   cp = np.cos(phi)
        s[o+1:o+int(k_[i])+1,0] = st*cp
        s[o+1:o+int(k_[i])+1,1] = st*sp
        s[o+1:o+int(k_[i])+1,2] = ct
        o += int(k_[i])

    return s


def DISCSe(N, bNisNtheta=False):

    # geometry factor
    alpha = np.deg2rad(30)
    epsilon = 1/np.cos(alpha/2) - 1

    # ideal number of latitudinal circles over entire sphere
    if bNisNtheta:
        n0 = N
        k_ = np.zeros(int(np.floor(n0*0.5+0.5)))
        n_ = n0
        for i,_ in enumerate(k_):
            if i+0.5<n0*0.5:
                k_[i] = round(np.pi/np.arcsin(np.sin(np.pi/n0/2)/np.sin((i+0.5)*np.pi/n0)))
            else:
                k_[i] = round(0.5*np.pi/np.arcsin(np.sin(np.pi/n0/2)/np.sin((i+0.5)*np.pi/n0)))
        N = int(sum(k_))
    elif N==1:
        n_ = 1
        k_ = np.array([1])
    else:
        n0 = np.pi / (2* np.arccos(-np.pi/(4*(1+epsilon)*(N-1)) + np.sqrt((np.pi/(4*(1+epsilon)*(N-1)))**2+1)))
        n0 = [int(np.floor(n0)), int(np.ceil(n0))]
    
        # init. final parameters
        m_ = np.inf             # best mean squared distance error    (1 x 1)
        n_ = None               # best number of latitudinal circles  (1 x 1)
        k_ = np.array([])       # best population of lat. circles     (n_ x 1)

        for n in n0:
            # n/2 (max. lat. circle in upper hemisphere)
            nhalf = int(np.floor(n*0.5+0.5))

            # mean n.n. distance
            d0 = 2*(1+epsilon)*np.sin(np.pi/2/n)
            # reasonable population range per latitude
            k = np.stack(nhalf * [np.arange(4*n-1)+2])
            # latitude index matrix of same size
            i = np.stack(k.shape[1] * [np.arange(nhalf)+1]).T
            # corresponding distances
            d = 2*np.sin((i-0.5)*np.pi/n)*np.sin(np.pi/k)
            if not (n-1)%2: # equator            
                d[-1,:] = 2*np.sin(i[-1,:]*np.pi/n)*np.sin(np.pi/k[-1,:]/2)

            # corresponding squared errors
            e = (d0-d)**2
            ind = np.argsort(e, axis=1)
            e = np.sort(e, axis=1)
            for i in range(nhalf):
                k[i,:] = k[i,ind[i,:]]

            # optimal parameters for given n(l):
            E = e[:,0]
            K = k[:,0]
            N0 = sum(K)

            # optimal population for given n -> best population for given N:
            i2=np.array([]); e2=np.array([]); k2=np.array([])
            Ndiff = N0-N
            if Ndiff!=0:
                for i in range(nhalf):
                    if Ndiff>0:
                        # indices available for population reduction
                        tmpj=np.where(k[i,:]<K[i])[0]
                    else:
                        # indices available for population increase
                        tmpj=np.where(k[i,:]>K[i])[0]

                    i2 = np.concatenate((i2,i*np.ones(tmpj.shape)), axis=0)
                    e2 = np.concatenate((e2,e[i,tmpj]), axis=0)
                    k2 = np.concatenate((k2,k[i,tmpj]), axis=0)       

                # change by Ndiff in total on those latitudes where error is least
                ind = np.argsort(e2)
                e2 = np.sort(e2)
                k2 = k2[ind]
                i2 = i2[ind]

                # update best parameters for given n and N
                E[np.int64(i2[:int(abs(Ndiff))])] = e2[:int(abs(Ndiff))]
                K[np.int64(i2[:int(abs(Ndiff))])] = k2[:int(abs(Ndiff))]

            # mean squared error
            m = sum(E*(K-1))/N

            # update final parameters for given N
            if m<m_:
                m_ = m
                n_ = n
                k_ = K.copy()
            
    # compute coordinates
    s = np.zeros([N,3]) # coordinate array
    o = -1
    for i in range(int(np.floor(n_*0.5+0.5))):
        j = (np.arange(k_[i])+1).T

        theta = (i+0.5)*np.pi/n_
        if i+1<(n_+1)*0.5:
            phi = (i+0.5+j)*2*np.pi/k_[i]
        else: # equator
            phi = (i+0.5+j)*np.pi/k_[i]

        st = np.sin(theta); ct = np.cos(theta)
        sp = np.sin(phi);   cp = np.cos(phi)
        s[o+1:o+int(k_[i])+1,0] = st*cp
        s[o+1:o+int(k_[i])+1,1] = st*sp
        s[o+1:o+int(k_[i])+1,2] = ct
        o += int(k_[i])

    return s
    