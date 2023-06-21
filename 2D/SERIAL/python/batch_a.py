# Python3 script to load and draw MHD-2D(Serial) data (all period)
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
        x=np.loadtxt(direc+"x.dat",dtype=float)
        y=np.loadtxt(direc+"y.dat",dtype=float)
        t=np.atleast_1d(np.loadtxt(direc+"t.dat",dtype=float))
        offs=np.loadtxt(direc+"offsets.dat",dtype=int)
        para=np.atleast_1d(np.loadtxt(direc+"params.dat",dtype=float))
        break
    except:
        print("Error during file load.")

xoff=offs[0]
yoff=offs[1]
gam=para[0]

#Number of elements
nx=np.size(x)
ny=np.size(y)
nt=np.size(t)
nd=8 #Number of dependent variables in MHD-2D

#Read MHD data @ all time
data=np.zeros((nt,nd,ny,nx),dtype=np.float64)
sst=0
for sst in range(0,nt):
    tmp=np.fromfile(direc+f"outdat_{sst:05d}.dat",dtype=np.float64).reshape((nd,ny,nx))
    data[sst,:,:,:]=tmp

# #Slice to remove ghost cells
x=x[xoff:nx-xoff]
y=y[yoff:ny-yoff]
data2=data[:,:,yoff:ny-yoff,xoff:nx-xoff]
dx=x[1]-x[0]
dy=y[1]-y[0]

#Primitive variables
ro=data2[:,0,:,:]
vx=data2[:,1,:,:]/ro
vy=data2[:,2,:,:]/ro
vz=data2[:,3,:,:]/ro
bx=data2[:,4,:,:]
by=data2[:,5,:,:]
bxct=bx
byct=by
#Cell-face to cell-center B
for i in range(0,nx-2*xoff):
    bx[:,:,i]=0.5*(data[:,4,yoff:ny-yoff,i+xoff]+data[:,4,yoff:ny-yoff,i+1+xoff])
for j in range(0,ny-2*yoff):
    by[:,j,:]=0.5*(data[:,5,j+yoff,xoff:nx-xoff]+data[:,5,j+1+yoff,xoff:nx-xoff])
bz=data2[:,6,:,:]
en=data2[:,7,:,:]
pr=(gam-1)*(en-0.5*(ro*(vx**2+vy**2+vz**2)+(bx**2+by**2+bz**2)))
data2=np.array([ro,vx,vy,vz,pr,bx,by,bz])
