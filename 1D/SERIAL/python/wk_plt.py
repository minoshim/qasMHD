# FFT plot for MHD wave

from python import plt2d

# Data size for FFT
nx=np.size(x)
nt=np.size(t)
ssx=2**(int(np.log(nx)/np.log(2)))
sst=2**(int(np.log(nt)/np.log(2)))

val=by
# val=by+bz*1j
# val=by-bz*1j
# val=pr
val=val-np.mean(val)

val=val[nt-sst:nt,nx-ssx:nx]
ff=np.fft.fftshift(np.fft.fft2(val)) # call 2D FFT
wnum=np.fft.fftshift(np.fft.fftfreq(ssx,x[1]-x[0]))
freq=np.fft.fftshift(np.fft.fftfreq(sst,t[1]-t[0]))
kk=wnum*2*np.pi
oo=freq*2*np.pi
ff=np.log10(np.abs(ff))

# Wave velocities
ro0=np.mean(ro[:,0])
pr0=np.mean(pr[:,0])
bx0=np.mean(bx[:,0])
bb0=np.sqrt(np.mean(bx[:,0]*bx[:,0]+by[:,0]*by[:,0]+bz[:,0]*bz[:,0]))
ca0=np.sqrt(bb0*bb0/ro0)
tmp=gam*pr0+bb0*bb0
cf0=np.sqrt(0.5*(tmp+np.sqrt(tmp*tmp-4*gam*pr0*bx0*bx0))/ro0)
cs0=np.sqrt(0.5*(tmp-np.sqrt(tmp*tmp-4*gam*pr0*bx0*bx0))/ro0)

# Plot option
xmin=0
xmax=np.max(kk)/2
ymin=0
ymax=np.max(oo)
zmax=np.max(ff)
zmin=zmax-4
xran=(xmin,xmax)
yran=(ymin,ymax)
zran=(zmin,zmax)

# Plot
fig=plt2d.image(kk,oo,ff)
plt.xlim(xran)
plt.ylim(yran)
plt.clim(zran)
plt.xlabel(r'$k_x$')
plt.ylabel(r'$\omega$')

plt.show(block=False)
