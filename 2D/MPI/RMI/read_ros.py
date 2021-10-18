import numpy as np

#Read MHD data (only density) @ all time
ros=np.zeros([nt,ny,nx])
for i in range(nt):
    data=np.fromfile(direc+f"merge_outdat_{i:05d}.dat",dtype=np.float32).reshape((nd,ny,nx))
    ros[i,:,:]=data[0,:,:]

# ssts=[10,15,20,25,30]
ssts=[5,10,15,20]
xran=[np.min(y),0]
xran=[-4,4]

plt.figure()
for i in range(len(ssts)):
    plt.plot(y,ros[ssts[i],:,nx//2],label=f"t={t[ssts[i]]:.1f}")

plt.legend()
plt.xlim(xran)
plt.show(block=False)

