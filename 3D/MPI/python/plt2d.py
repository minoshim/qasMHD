import matplotlib.pyplot as plt

def image(x,y,val,save=0,title="",cmap="jet",xlabel="x",ylabel="y"):
    fig=plt.figure()
    plt.axes().set_aspect("equal")
    plt.pcolormesh(x,y,val,cmap=cmap)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.colorbar()
    if (save):
        fig.savefig("result.eps")
    plt.show(block=False) # console non-blocked
    return 0

def imagex3(x,y,z,val,save=0,title="",cmap="jet"):
    fig=plt.figure(figsize=(9,3))
    nx=len(x)
    ny=len(y)
    nz=len(z)

    plt.subplot(1,3,1)
    tmp=0.5*(val[nz//2-1,:,:]+val[nz//2,:,:])
    plt.pcolormesh(x,y,tmp,cmap=cmap)
    plt.xlabel("x")
    plt.ylabel("y")

    plt.subplot(1,3,2)
    tmp=0.5*(val[:,:,nx//2-1]+val[:,:,nx//2])
    plt.pcolormesh(y,z,tmp,cmap=cmap)
    plt.xlabel("y")
    plt.ylabel("z")
    plt.title(title)
    
    plt.subplot(1,3,3)
    tmp=0.5*(val[:,ny//2-1,:]+val[:,ny//2,:])
    plt.pcolormesh(x,z,tmp,cmap=cmap)
    plt.xlabel("x")
    plt.ylabel("z")

    fig.tight_layout()
    plt.show(block=False) # console non-blocked
    return 0
