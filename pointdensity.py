# -*- coding: utf-8 -*-
"""
(C) 2015, Niels Anders, EUBrazilCC
"""
from generic_tools import getpoints, gridcellborders, saveimg
import numpy as np
import sys, os.path, time

def createPD(x,y,pd,bx,by, printYN='no'):
    if printYN=='yes':
        print 'Calculate point density (%d x %d)...' % (len(xi), len(yi)),
        time.sleep(0.1)
        t0 = time.time()
    for i in np.arange(len(bx)-1):
        col = np.argwhere((x > bx[i]) & (x < bx[i+1]))
        for j in np.arange(len(by)-1):
            row = (y[col] >= by[j]) & (y[col] < by[j+1])
            if sum(row) > 0:           
                pd[j,i] = len(y[col[row]])/float(res**2)
    # return
    if printYN=='yes':
        t1 = time.time()
        print 'finished in %1d seconds' % (t1-t0)
        time.sleep(0.1)    
    return pd
    
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
    fn  = basename+'_pointdensity.tif'
    
    # read point cloud
    x,y,z,c = getpoints(In)
    
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
    
    # Calculate point density
    pd = createPD(x,y,grid,bx,by)
    
    # write tiff
    saveimg(np.flipud(pd), fn, len(xi), len(yi), geotransform, proj)
