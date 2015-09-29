# -*- coding: utf-8 -*-
"""
(C) 2015, Niels Anders, EUBrazilCC
"""
from generic_tools import getpoints, gridcellborders, saveimg
import numpy as np
import sys, os.path, time

def createRH(x,y,z,rh,dtm,bx,by,pc):
    ny, nx = rh.shape
    print 'Calculate RH%2d (%d x %d)...' % (pc, nx, ny),
    time.sleep(0.1)
    t0 = time.time()
    for i in np.arange(len(bx)-1):
        col = np.argwhere((x > bx[i]) & (x < bx[i+1]))
        for j in np.arange(len(by)-1):
            row = (y[col] >= by[j]) & (y[col] < by[j+1])
            if sum(row) > 0:           
                rh[j,i] = np.percentile(z[col[row]] - dtm[j,i], pc)
    # return
    t1 = time.time()
    print 'finished in %1d seconds' % (t1-t0)
    time.sleep(0.1)    
    return rh
    
if __name__=='__main__':
    total   = len(sys.argv)
    cmdargs = str(sys.argv)

    if total == 3:      
        pc          = int(sys.argv[1])
        filename    = sys.argv[2]
        
    elif total == 2:
        filename = sys.argv[1]
        pc = 50
    else:
        print 'wrong number of arguments'
        sys.exit()
        
    basename = os.path.splitext(filename)[0]
    extension = os.path.splitext(filename)[1]

    In  = filename
    fn  = basename+'_rh%d.tif' % pc
    
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
    
    # create dtm
    from dtm import createDTM  
    g = c == 2
    dtm = createDTM(x[g],y[g],z[g],xi,yi,res,geotransform, proj, method='idw')    
    
    # Calculate relative height percentiles
    rh = createRH(x,y,z,grid,dtm,bx,by,pc)
    
    # write tiff
    saveimg(np.flipud(rh), fn, len(xi), len(yi), geotransform, proj)
