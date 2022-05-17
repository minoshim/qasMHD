# FFT plot for MHD wave

from python import plt2d

# Data size for FFT
nx=np.size(x)
nt=np.size(t)
ssx=2**(int(np.log(nx)/np.log(2)))
sst=2**(int(np.log(nt)/np.log(2)))

val=by

val=val[nt-sst:nt,nx-ssx:nx]
ff=np.fft.fftshift(np.fft.fft2(val)) # call 2D FFT
wnum=np.fft.fftshift(np.fft.fftfreq(ssx,x[1]-x[0]))
freq=np.fft.fftshift(np.fft.fftfreq(sst,t[1]-t[0]))
kk=wnum*2*np.pi
oo=freq*2*np.pi
ff=np.log10(np.abs(ff))

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
