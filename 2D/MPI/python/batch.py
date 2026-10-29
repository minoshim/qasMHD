# Python3 script to load and draw MHD-2D(MPI-merged) data
# Packages Numpy and Matplotlib are required.

# Call the script in command line:
# > python batch.py
# Call the script in Python3 interactive mode:
# >>> exec(open("batch.py").read())

import numpy as np
import matplotlib.pyplot as plt
from python import plt2d
from python.mhd_io import read_ct_fields

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

#Optional, time-independent gravitational potential (ghost cells already removed).
potential_path=direc+"merge_g_potential.dat"
try:
    with open(potential_path,"rb") as potential_file:
        potential_bytes=potential_file.read()
except FileNotFoundError:
    phi_g=None
    print("No merged gravitational potential: pressure is not gravity-corrected.")
else:
    expected_bytes=nx*ny*np.dtype(np.float32).itemsize
    if len(potential_bytes) != expected_bytes:
        raise ValueError(f"{potential_path}: expected {expected_bytes} bytes, got {len(potential_bytes)}")
    phi_g=np.frombuffer(potential_bytes,dtype=np.float32).reshape((ny,nx))
    if not np.all(np.isfinite(phi_g)):
        raise ValueError(f"{potential_path}: gravitational potential contains non-finite values")
    print(f"Loaded gravitational potential from {potential_path}")

#Read MHD data @ particular time
sst=-1
while ((sst < 0) or (sst >= nt)):
    sst=int(input(f"Specity time period (0-{nt-1}): "))
    
data=np.fromfile(direc+f"merge_outdat_{sst:05d}.dat",dtype=np.float32).reshape((nd,ny,nx))
bxct,byct=read_ct_fields(direc,sst,nx,ny)

dx=x[1]-x[0]
dy=y[1]-y[0]

#Primitive variables
ro=data[0,:,:]
vx=data[1,:,:]/ro
vy=data[2,:,:]/ro
vz=data[3,:,:]/ro
bx=data[4,:,:]                  # Cell centers; already averaged by merge.out
by=data[5,:,:]                  # Do not interpolate again in Python
bz=data[6,:,:]
en=data[7,:,:]
if phi_g is not None:
    en=en-ro*phi_g #Do not modify the stored total energy in place.
pr=(gam-1)*(en-0.5*(ro*(vx**2+vy**2+vz**2)+(bx**2+by**2+bz**2)))

# CT divergence at all cell centers (second-order diagnostic).
divb=(bxct[:,1:]-bxct[:,:-1])/dx+(byct[1:,:]-byct[:-1,:])/dy
# Current at interior corners only; outer corners need tangential ghost data.
# jz.shape=(ny-1,nx-1); xjz,yjz are its coordinates. bx/by remain cell-centered.
xjz=0.5*(x[1:]+x[:-1])
yjz=0.5*(y[1:]+y[:-1])
jz=(byct[1:-1,1:]-byct[1:-1,:-1])/dx-(bxct[1:,1:-1]-bxct[:-1,1:-1])/dy

# #Plot
val=pr/ro
a=plt2d.image(x=x,y=y,val=val,save=0,title=f"t={t[sst]:.2f}",show=1)
