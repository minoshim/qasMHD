import lic
import matplotlib.pyplot as plt
from . import plt2d

def lic_image(x,y,vx,vy,
              length=30,
              figsize=(6.4,4.8),
              aspect="auto",
              cmap="gray",
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
              colorbar=0,
              filename="result.eps"):

    lic_result=lic.lic(vx.transpose(),vy.transpose(),length=length)
    plt2d.image(x,y,lic_result.transpose(),
                figsize=figsize,
                aspect=aspect,
                cmap=cmap,
                title=title,
                xlabel=xlabel,
                ylabel=ylabel,
                vmin=vmin,
                vmax=vmax,
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
                equal=equal,
                show=show,
                save=save,
                colorbar=colorbar,
                filename=filename)

    return lic_result
