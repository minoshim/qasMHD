# Python3 script to load and draw MHD-2D(MPI-merged) data (all period)
# Packages Numpy and Matplotlib are required.

# Call the script in command line:
# > python batch_a.py
# Call the script in Python3 interactive mode:
# >>> exec(open("batch_a.py").read())

import numpy as np
import matplotlib.pyplot as plt
from python import plt2d

#Read independent variables and parameters
while True:
    direc=input("Input data directory (Ctrl-D to exit): ")+"/"
    try:
        x=np.loadtxt(direc+"merge_x.dat",dtype=float)
        y=np.loadtxt(direc+"merge_y.dat",dtype=float)
        t=np.atleast_1d(np.loadtxt(direc+"t.dat",dtype=float))
        para=np.atleast_1d(np.loadtxt(direc+"params.dat",dtype=float))
        break
    except:
        print("Error during file load.")

gam=para[0]

#Number of elements
nx=np.size(x)
ny=np.size(y)
nt=np.size(t)
nd=8 #Number of dependent variables in MHD-2D

#Read MHD data @ all time
data=np.zeros((nt,nd,ny,nx),dtype=np.float32)
sst=0
for sst in range(0,nt):
    tmp=np.fromfile(direc+f"merge_outdat_{sst:05d}.dat",dtype=np.float32).reshape((nd,ny,nx))
    data[sst,:,:,:]=tmp

dx=x[1]-x[0]
dy=y[1]-y[0]

#Primitive variables
ro=data[:,0,:,:]
vx=data[:,1,:,:]/ro
vy=data[:,2,:,:]/ro
vz=data[:,3,:,:]/ro
bx=data[:,4,:,:]                  # @ CT grid (i-1/2,j)
by=data[:,5,:,:]                  # @ CT grid (i,j-1/2)
bz=data[:,6,:,:]
en=data[:,7,:,:]
pr=(gam-1)*(en-0.5*(ro*(vx**2+vy**2+vz**2)+(bx**2+by**2+bz**2)))

# Current
jz=np.zeros((nt,ny,nx))
for j in range(1,ny):
    for i in range(1,nx):
        jz[:,j,i]=(by[:,j,i]-by[:,j,i-1])/dx-(bx[:,j,i]-bx[:,j-1,i])/dy # @ corner (i-1/2,j-1/2)
