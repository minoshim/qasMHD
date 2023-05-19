from mayavi import mlab
import numpy as np

# 3D volume rendering using mayavi.mlab

def volume(x0,y0,z0,val,
           contours=10,
           transparent=False,
           opacity=1.0,
           colormap='jet',
           vmin=None,
           vmax=None,
           xmin=None,
           xmax=None,
           ymin=None,
           ymax=None,
           zmin=None,
           zmax=None,
           color=None,
           line_width=2.0,
           reset_zoom=True,
           figure=None):

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
        
    mlab.contour3d(X,Y,Z,val.transpose(),
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
