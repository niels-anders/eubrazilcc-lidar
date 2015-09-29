# -*- coding: utf-8 -*-
"""
(C) 2015, Niels Anders, EUBrazilCC
"""
from generic_tools import saveimg, getpoints, gridcellborders
import numpy as np
import sys, time, os.path

def createDTM(x, y, z, xi, yi, res, geotransform, proj, method='idw'):    
    print 'Create DTM (%d x %d)...' % (len(xi), len(yi)),
    time.sleep(0.1)
    t0 = time.time()
    
    XI, YI = np.meshgrid(xi, yi, indexing='xy')

    nnear = 25
    if method=='idw':
        from generic_tools import idw
        dtm = idw(x,y,z,XI.flatten(),YI.flatten(), nnear).reshape(XI.shape)
    else:
        from scipy import interpolate
        dtm = interpolate.griddata((x,y),z,(XI,YI), method=method)
    
    dtm[dtm==-999] = np.NaN

    # return
    t1 = time.time()
    print 'finished in %2d seconds' % (t1-t0)
    time.sleep(0.1)
    return dtm  
    
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
    fn  = basename+'_dtm.tif'
    
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
    proj = '+proj=utm +zone=21 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs' # mandatory setting used to store projection information in metadata geotiff (not assigned as metadata is not stored in lidar file)
    
    # create DTM
    dtm = createDTM(x[g],y[g],z[g],xi,yi,res,geotransform, proj, method='idw')

    # write tiff
    saveimg(np.flipud(dtm), fn, len(xi), len(yi), geotransform, proj)
    #string = str('gdal_translate '+ fn+' D:/testdtm.png -of PNG -scale') 
    #print string    
    #exec(string)
