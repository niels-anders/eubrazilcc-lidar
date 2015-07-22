# -*- coding: utf-8 -*-
from generic_tools import getpoints, gridcellborders, saveimg
import numpy as np
import sys, os.path, time

def createCHM(dtm,dsm):
    ny, nx = dtm.shape
    print 'Calculate CHM (%d x %d)...' % (nx, ny),
    time.sleep(0.1)
    t0 = time.time()
    chm = dsm - dtm

    # return
    t1 = time.time()
    print 'finished in %1d seconds' % (t1-t0)
    time.sleep(0.1)    
    return chm
    
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
    fn  = '../example_data/output'+basename[basename.rfind('/'):]+'_chm.tif' 
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
    
    # create dsm
    from dsm import createDSM
    dsm = createDSM(x,y,z,xi,yi,bx,by,res,geotransform, proj)
    
    # create chm
    chm = createCHM(dtm,dsm)
    
    # write tiff
    saveimg(np.flipud(chm), fn, len(xi), len(yi), geotransform, proj)
    