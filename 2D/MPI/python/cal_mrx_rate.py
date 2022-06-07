# Calculate MRX rate

import numpy as np
from scipy import integrate

# Magnetic flux integrated along Y
mflux=np.zeros((nt,nx))
for n in range(nt):
    for i in range(nx):
        tmp=np.abs(bx[n,:,i])
        mflux[n,i]=integrate.simps(tmp,y)

# Identify reconnection point at which Magnetic flux is minimum
smx=[0]*nt
for n in range(nt):
    smx[n]=np.where(mflux[n,:] == np.min(mflux[n,:]))[0][-1]
    if (smx[n] == 0):
        smx[n]=1
    if (smx[n] == nx-1):
        smx[n]=nx-2

# Reconnection rate
mrate=np.zeros(nt)
for n in range(1,nt):
    mrate[n]=-0.5*(mflux[n,smx[n]]-mflux[n-1,smx[n-1]])/(t[n]-t[n-1])
    # When MRX in full region, factor 0.5 is needed.

# Time integration of reconnection rate
int_mrate=np.append(0,integrate.cumtrapz(mrate,t))

# Normalization
mrate=mrate/(np.abs(bx[0,ny-1,nx//2])/np.sqrt(ro[0,ny-1,nx//2]))
