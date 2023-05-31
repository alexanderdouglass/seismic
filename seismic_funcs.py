

import numpy as np
from matplotlib import pyplot as plt

from scipy.sparse import spdiags
from scipy.sparse.linalg import lsqr as splsqr
from spgl1.lsqr import lsqr
from spgl1 import spgl1, spg_lasso, spg_bp, spg_bpdn, spg_mmv
from spgl1.spgl1 import norm_l1nn_primal, norm_l1nn_dual, norm_l1nn_project
from spgl1.spgl1 import norm_l12nn_primal, norm_l12nn_dual, norm_l12nn_project


### Functions ###

# Basis pursuit denoising function
def BPDN(A,b,sigma=0.1,iters=1000, verbs=1):
    # Sigma is desired ||Ax - b||_2

    r,resid,grad,info = spg_bpdn(A, b, sigma, iter_lim=iters, verbosity=verbs)
    # r,resid,grad,info = spg_bp(A, b)

    r = np.array( np.flip( r ) )

    return r

# Reflection inversion with BPDN
def RefInv(data,srcSig,freqs,Fs=500,sigma=0.1,iters=1000):

    N = len(data)
    tMax = int(N/Fs)
    rsave = np.zeros(N)
    t = np.arange(N)/Fs
    dt = 1/Fs
    f1 = freqs[0] # Initial frequency
    f2 = freqs[1] # Final frequency

    if np.mod(N,2) == 1:
        N += 1

    df = Fs/N
    t_bp = np.arange(0,N/Fs,dt)
    f_bp = np.linspace(0,Fs,N)

    if len(srcSig) < N:
        srcSig = np.append(srcSig,np.zeros(N-len(srcSig)))
    else:
        srcSig = srcSig[0:N]
    SRCSIG = np.fft.fft(srcSig,N)

    ARRAY = np.fft.fft(data,N)

    K = int((f2-f1)/df) # freq samples
    if np.mod(K,2) == 1:
        K += 1
    K1 = int(f1/df)
    K2 = int(f2/df) + 1

    # Define I, W
    W = 2*SRCSIG[K1:K2]
    fi = np.vstack(f_bp[K1:K2])
    Tj = t[0:N:2]
    dphi = (N-1)*dt/2
    Wr = np.real(W)
    Wi = np.imag(W)
    mij = np.pi*fi*Tj
    ni = 2*np.pi*fi*dphi

    # Define A
    A11 =  Wr*(np.sin(mij)*np.sin(ni)).T - Wi*(np.sin(mij)*np.cos(ni)).T
    A12 =  Wr*(np.cos(mij)*np.cos(ni)).T + Wi*(np.cos(mij)*np.sin(ni)).T
    A21 =  Wr*(np.sin(mij)*np.cos(ni)).T + Wi*(np.sin(mij)*np.sin(ni)).T
    A22 = -Wr*(np.cos(mij)*np.sin(ni)).T + Wi*(np.cos(mij)*np.cos(ni)).T
    A = np.bmat([[A11.T, A12.T], [A21.T, A22.T]])

    # Define b
    b = np.append(np.real(ARRAY[K1:K2]),np.imag(ARRAY[K1:K2]))

    # Solve for r
    r = BPDN(A,b,sigma*np.linalg.norm(b,2),iters=iters,verbs=0)
    ro = r[0:int(N/2)]
    re = r[int(N/2):]
    r1 = re + ro
    r2 = -np.flip(re - ro)
    rout = np.append(r1,r2)

    return rout


### PLOTS ###

# def waterfall(data,x,y,dy,ampFac=0,xLims=[0,1],yLims[0,1]):
    
#     dataTmp = 0*data.copy()
#     if ampFac == 0:
#         ampFac = np.max(dataTmp[i,:])
        
#     for i in range(len(y)):
#         dataTmp[:,i] = ampFac*data[:,i] + y[i]/1e3
        
#     plt.plot(t,datatmp[:,rec_plot],color)
#     plt.xlabel('Time (s)')
#     plt.ylabel('Range (km)')
#     plt.xlim(xlim)
#     plt.ylim(ylim)
#     plt.gca().invert_yaxis()
#     plt.show()
