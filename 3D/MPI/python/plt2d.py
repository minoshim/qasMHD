import matplotlib.pyplot as plt

def image(x,y,val,
          figsize=(6.4,4.8),
          aspect="auto",
          cmap="jet",
          title="",
          xlabel="",
          ylabel="",
          vmin=None,
          vmax=None,
          xmin=None,
          xmax=None,
          ymin=None,
          ymax=None,
          equal=0,
          show=0,
          save=0,
          colorbar=1,
          filename="result.eps"):
    
    fig=plt.figure(figsize=figsize)

    if (equal != 0):
        aspect="equal"
    plt.axes().set_aspect(aspect)

    if (vmin == None):
        vmin=val.min()
    if (vmax == None):
        vmax=val.max()
        
    # plt.pcolormesh(x,y,val,cmap=cmap,shading='auto')
    axesimage=plt.imshow(val,cmap=cmap,
                         vmin=vmin,vmax=vmax,
                         extent=[x.min(),x.max(),y.min(),y.max()],
                         interpolation='nearest',origin='lower',aspect=aspect) # Fast

    if (xmin == None):
        xmin=x.min()
    if (xmax == None):
        xmax=x.max()
    if (ymin == None):
        ymin=y.min()
    if (ymax == None):
        ymax=y.max()
    plt.xlim((xmin,xmax))
    plt.ylim((ymin,ymax))
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)

    if (colorbar):
        plt.colorbar()
    if (save):
        fig.savefig(filename)
    if (show != 0):
        plt.show(block=False) # console non-blocked

    return fig, axesimage

def imagex3(x,y,z,val,save=0,title="",cmap="jet",xlabel="x",ylabel="y",zlabel="z",figsize=(12.0,4.0),filename="result3.eps",equal=0,show=0,aspect="auto"):
    fig=plt.figure(figsize=figsize)

    if (equal != 0):
        aspect="equal"
    plt.axes().set_aspect(aspect)

    nx=len(x)
    ny=len(y)
    nz=len(z)

    plt.subplot(1,3,1)
    tmp=0.5*(val[nz//2-1,:,:]+val[nz//2,:,:])
    # plt.pcolormesh(x,y,tmp,cmap=cmap,shading='auto')
    plt.imshow(tmp,cmap=cmap,
               extent=[x.min(),x.max(),y.min(),y.max()],
               interpolation='nearest',origin='lower',aspect=aspect) # Fast
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.subplot(1,3,2)
    tmp=0.5*(val[:,:,nx//2-1]+val[:,:,nx//2])
    # plt.pcolormesh(y,z,tmp,cmap=cmap,shading='auto')
    plt.imshow(tmp,cmap=cmap,
               extent=[y.min(),y.max(),z.min(),z.max()],
               interpolation='nearest',origin='lower',aspect=aspect) # Fast
    plt.xlabel(ylabel)
    plt.ylabel(zlabel)
    plt.title(title)
    
    plt.subplot(1,3,3)
    tmp=0.5*(val[:,ny//2-1,:]+val[:,ny//2,:])
    # plt.pcolormesh(x,z,tmp,cmap=cmap,shading='auto')
    plt.imshow(tmp,cmap=cmap,
               extent=[x.min(),x.max(),z.min(),z.max()],
               interpolation='nearest',origin='lower',aspect=aspect) # Fast
    plt.pcolormesh(x,z,tmp,cmap=cmap,shading='auto')
    plt.xlabel(xlabel)
    plt.ylabel(zlabel)

    fig.tight_layout()
    if (save):
        fig.savefig(filename)
    if (show != 0):
        plt.show(block=False) # console non-blocked
    return fig
