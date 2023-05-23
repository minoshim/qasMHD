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
