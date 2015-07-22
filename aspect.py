# -*- coding: utf-8 -*-
from generic_tools import saveimg, getpoints, gridcellborders
import numpy as np
import sys, time, os.path

def calcAspect(dtm):
    [r,c] = dtm.shape
    print 'Calculate slope angle (%d x %d)...' % (c, r),
    time.sleep(0.1)
    t0 = time.time()
    G = np.zeros((r,c,9))
    # Calculate gradient from center cell to each surrounding cell
    G[1:r-1,1:c-1,0] = (dtm[1:r-1,1:c-1] - dtm[0:r-2,0:c-2])/(2**0.5)
    G[1:r-1,1:c-1,1] = (dtm[1:r-1,1:c-1] - dtm[1:r-1,0:c-2])
    G[1:r-1,1:c-1,2] = (dtm[1:r-1,1:c-1] - dtm[2:r  ,0:c-2])/(2**0.5)
    G[1:r-1,1:c-1,3] = (dtm[1:r-1,1:c-1] - dtm[0:r-2,1:c-1])
    G[1:r-1,1:c-1,4] = 0
    G[1:r-1,1:c-1,5] = (dtm[1:r-1,1:c-1] - dtm[2:r  ,1:c-1])
    G[1:r-1,1:c-1,6] = (dtm[1:r-1,1:c-1] - dtm[0:r-2,2:c  ])/(2**0.5)
    G[1:r-1,1:c-1,7] = (dtm[1:r-1,1:c-1] - dtm[1:r-1,2:c  ])
    G[1:r-1,1:c-1,8] = (dtm[1:r-1,1:c-1] - dtm[2:r  ,2:c  ])/(2**0.5)
    
    Gmax = G.max(axis=2)                # max gradient
    slope = np.arctan(Gmax)/(np.pi/180)   # convert to slope angle
    aspect = np.zeros((r,c))*-1
    for i in np.arange(r):
        for j in np.arange(c):
            if Gmax[i,j] > 0:
                aspect[i,j] = np.argwhere(G[i,j,:]==Gmax[i,j])*45
    # return
    t1 = time.time()
    print 'finished in %1d seconds' % (t1-t0)
    time.sleep(0.1)  
    return aspect

if __name__=='__main__':
    total   = len(sys.argv)
    cmdargs = str(sys.argv)

    if total != 2:      
      print ("The total numbers of args passed to the script should be 1")
      sys.exit()

    filename = sys.argv[1]
    basename = os.path.splitext(filename)[0]
    extension = os.path.splitext(filename)[1]

    In  = filename
    fn  = basename+'_aspect.tif'
    
    # read point cloud
    x,y,z,c = getpoints(In)
    g = c==2    # ground points
    
    # intiate grid system (coordinates of center point grid cells)
    res = 1
    xi, yi = np.arange(x.min(), x.max()+res/2, res), np.arange(y.min(), y.max()+res/2, res)
    extent = [xi.min(), xi.max(), yi.min(), yi.max()] # only used to set extent to plots
    grid = np.zeros((len(yi),len(xi)))*np.NaN  

    # retreive gridcell border coordinates
    bx, by = gridcellborders(xi,yi,res)    
    
    # georeferencing info
    geotransform = (x.min(), res, 0, y.max(), 0, -res)  # used for georeference metadata in geotiff
    proj = '' # mandatory setting used to store projection information in metadata geotiff (not assigned as metadata is not stored in lidar txt)
    
    # create DTM
    from dtm import createDTM
    dtm = createDTM(x[g],y[g],z[g],xi,yi,res,geotransform, proj, method='idw')

    # calculate slope angle
    aspect = calcAspect(dtm) 

    # write tiff
    saveimg(np.flipud(aspect), fn, len(xi), len(yi), geotransform, proj)
    