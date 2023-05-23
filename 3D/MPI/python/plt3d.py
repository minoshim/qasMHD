from mayavi import mlab
import numpy as np

# 3D volume rendering using mayavi.mlab

# Some default parameters
sizedef=(800,700) #Figure size
opadef=1.0 # Opacity
cmapdef="jet" # Colormap
lwdef=2.0 # Line width
nbldef=2 # Number of labels for axes

def volume(x0,y0,z0,val,
           size=sizedef,
           contours=10,
           transparent=False,
           opacity=opadef,
           colormap=cmapdef,
           vmin=None,
           vmax=None,
           xmin=None,
           xmax=None,
           ymin=None,
           ymax=None,
           zmin=None,
           zmax=None,
           color=None,
           line_width=lwdef,
           reset_zoom=True,
           figure=None,
           nb_labels=nbldef,
           xlabel="",
           ylabel="",
           zlabel=""):

    if (figure == None):
        figure=mlab.figure(size=size)

    dx=x0[1]-x0[0]
    dy=y0[1]-y0[0]
    dz=z0[1]-z0[0]
    X, Y, Z=np.mgrid[np.min(x0):np.max(x0)+dx:dx,
                     np.min(y0):np.max(y0)+dy:dy,
                     np.min(z0):np.max(z0)+dz:dz]
    
    if (xmin == None):
        xmin=np.min(x0)
    if (xmax == None):
        xmax=np.max(x0)
    if (ymin == None):
        ymin=np.min(y0)
    if (ymax == None):
        ymax=np.max(y0)
    if (zmin == None):
        zmin=np.min(z0)
    if (zmax == None):
        zmax=np.max(z0)
        
    obj=mlab.contour3d(X,Y,Z,val.transpose(),
                       contours=contours,
                       transparent=transparent,
                       opacity=opacity,
                       colormap=colormap,
                       vmin=vmin,
                       vmax=vmax,
                       extent=(xmin,xmax,ymin,ymax,zmin,zmax),
                       color=color,
                       line_width=line_width,
                       reset_zoom=reset_zoom,
                       figure=figure)
    mlab.axes(figure=figure,nb_labels=nb_labels,xlabel=xlabel,ylabel=ylabel,zlabel=zlabel)
    mlab.outline(figure=figure)
    return figure, obj

def slice(x0,y0,z0,val,
          size=sizedef,
          plane_orientation='x_axes',
          slice_index=0,
          transparent=False,
          opacity=opadef,
          colormap=cmapdef,
          vmin=None,
          vmax=None,
          xmin=None,
          xmax=None,
          ymin=None,
          ymax=None,
          zmin=None,
          zmax=None,
          color=None,
          line_width=lwdef,
          reset_zoom=True,
          figure=None,
          nb_labels=nbldef,
          xlabel="",
          ylabel="",
          zlabel=""):
    
    if (figure == None):
        figure=mlab.figure(size=size)

    dx=x0[1]-x0[0]
    dy=y0[1]-y0[0]
    dz=z0[1]-z0[0]
    X, Y, Z=np.mgrid[np.min(x0):np.max(x0)+dx:dx,
                     np.min(y0):np.max(y0)+dy:dy,
                     np.min(z0):np.max(z0)+dz:dz]

    if (xmin == None):
        xmin=np.min(x0)
    if (xmax == None):
        xmax=np.max(x0)
    if (ymin == None):
        ymin=np.min(y0)
    if (ymax == None):
        ymax=np.max(y0)
    if (zmin == None):
        zmin=np.min(z0)
    if (zmax == None):
        zmax=np.max(z0)
        
    obj=mlab.volume_slice(X,Y,Z,val.transpose(),
                          plane_orientation=plane_orientation,
                          slice_index=slice_index,
                          transparent=transparent,
                          opacity=opacity,
                          colormap=colormap,
                          vmin=vmin,
                          vmax=vmax,
                          extent=(xmin,xmax,ymin,ymax,zmin,zmax),
                          color=color,
                          line_width=line_width,
                          reset_zoom=reset_zoom,
                          figure=figure)
    mlab.axes(figure=figure,nb_labels=nb_labels,xlabel=xlabel,ylabel=ylabel,zlabel=zlabel)
    mlab.outline(figure=figure)
    return figure, obj 
