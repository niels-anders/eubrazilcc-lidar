# -*- coding: utf-8 -*-
"""
(C) 2015, Niels Anders, EUBrazilCC
"""
from generic_tools import saveimg, getpoints, gridcellborders, ETA
import numpy as np
import sys, time, os.path

def createDSM(x,y,z,xi,yi,bx,by,res,geotransform, proj):
    """
    Create DSM based on max point elevation
    """
    from generic_tools import idw
    print 'Create DSM (%d x %d)...' % (len(xi), len(yi))
    time.sleep(0.1)
    t0 = time.time()
    
    # initialize dsm array
    dsm = np.zeros((len(yi), len(xi))) * np.NaN
    
    step = 0
    for i in np.arange(len(xi)):
        col = np.argwhere((x > bx[i]) & (x < bx[i+1]))
        for j in np.arange(len(yi)):
            row = (y[col] >= by[j]) & (y[col] < by[j+1])
            if sum(row) > 0:           
                dsm[j,i] = z[col[row]].max()
        step = ETA(t0,time.time(),step,i,0,len(xi))                
            
    # fill gaps using IDW interpolation ---------------------------
    XI, YI = np.meshgrid(xi,yi, indexing='xy')
    values = np.isnan(dsm)==0
    nodata = np.isnan(dsm)==1    
    dsm[nodata] = idw(XI[values], YI[values], dsm[values], XI[nodata], YI[nodata], 25)
         
    # return
    t1 = time.time()
    print 'finished in %1d seconds' % (t1-t0)
    time.sleep(0.1)
    
    return dsm    
    
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
    fn  = basename+'_dsm.tif'
    
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
    
    # create DSM
    dsm = createDSM(x,y,z,xi,yi,bx,by,res,geotransform, proj)

    # write tiff
    saveimg(np.flipud(dsm), fn, len(xi), len(yi), geotransform, proj)
