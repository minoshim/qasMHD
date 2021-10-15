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
