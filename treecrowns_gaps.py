# -*- coding: utf-8 -*-
"""
(C) 2015, Niels Anders, EUBrazilCC
"""
from generic_tools import getpoints, gridcellborders, saveimg, ETA, openimg
import numpy as np
import sys, os.path, time
import os
from osgeo import ogr, gdal

def polygonize(tiff, shp):
    if os.path.isfile(shp) == True:    
        try: 
            basename = os.path.splitext(shp)[0]
            os.remove(basename+'.shp')
            os.remove(basename+'.dbf')
            os.remove(basename+'.shx')
        except: 
            print "Can't delete shapefile"
    src = gdal.Open(tiff)
    srcband = src.GetRasterBand(1)
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dst_ds = driver.CreateDataSource(shp)
    dst_Layer = dst_ds.CreateLayer(shp, srs=None)
    gdal.Polygonize(srcband, None, dst_Layer, -1, [], callback=None)
    print 'File saved as', shp
    
def ndlocalmaxima(data,s=2):
    """
    data    = 2D raster
    s       = search radius (nr of cells)
    """
    windowsize = s*2+1
    centercell = windowsize*s+s
    nx, ny = data.shape
    localmax = np.zeros((nx,ny))
    for i in np.arange(s,nx):
        for j in np.arange(s,ny):
            z = data[i-s:i+s+1,j-s:j+s+1]
            if np.argwhere(z.flatten()==z.max())[0][0]==centercell:
                localmax[i][j] = 1
    return localmax

def treecrowns(chm, seed, crowns, seeds, old_seeds, label, top, s_max, xi, yi, threshold):
    if old_seeds[seed[0], seed[1]] == 0:
        # only continue if seed is not already processed
        crowns[seed[0], seed[1]] = label
        surroundings = chm[seed[0]-1:seed[0]+2, seed[1]-1:seed[1]+2]
        crown = np.argwhere(crowns == label)

        for i in range(surroundings.shape[0]):
            for j in range(surroundings.shape[1]):
                # evaluate merging of object
                global_i = seed[0] + i -1
                global_j = seed[1] + j -1
                temp = np.vstack([crown, np.array([global_i, global_j])]).astype('int')                  
                std = np.std(chm[temp[:,0], temp[:,1]])        
                if  std < threshold:
                    # accept merging
                    crown = np.vstack([crown, np.array([global_i, global_j])])  
                    crowns[global_i, global_j] = label 
        # update seed list
        seeds = seeds[1:,:]
        old_seeds[seed[0], seed[1]] = label
        for row in crown:
            # if seed does not already exist or already processed
            if (old_seeds[row[0], row[1]] == 0) & ((seeds == row).any() == False):
                distance = np.sqrt((top[0] - Xi[row[0], row[1]])**2+(top[1] - Yi[row[0], row[1]])**2)
                # if distance is within treshold tree height / crown width regression                
                if distance < s_max:
                    seeds = np.vstack([seeds, np.array(row)])
    else:
        # update seed list
        seeds = seeds[1:,:]
    return crowns, seeds

def crownheights(chm, crowns):
    """nog niet klaar"""
    chmf = chm.flatten()
    crownsf = crowns.flatten()    
    
    ids = np.unique(crownsf)
    ids = ids[np.isnan(ids)==0]
    ch = np.zeros(crownsf.shape)* np.nan
    for i in ids:
        n = np.argwhere(crowns==i)
        ch[n] = chmf[n].max()
    ch = ch.reshape(crowns.shape)
    return ch
        
    
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
    fn_tree     = basename+'_treecrowns.tif' 
    fn_gaps     = basename+'_forestgaps.tif' 
    fn_heights  = basename+'_treeheights.txt' 
    fn_peaks    = basename+'_localpeaks.txt' 
    fn_shape    = basename+'_treecrowns.shp' 
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
    try:
        dtm, _, _, _, _, _ = openimg(basename+'_dtm.tif')
    except:        
        from dtm import createDTM  
        g = c == 2
        dtm = createDTM(x[g],y[g],z[g],xi,yi,res,geotransform, proj, method='idw')    
    
    # create dsm
    try:
        dtm, _, _, _, _, _ = openimg(basename+'_dsm.tif')
    except:        
        from dsm import createDSM
        dsm = createDSM(x,y,z,xi,yi,bx,by,res,geotransform, proj)
    
    # create chm
    try:
        chm, _, _, _, _, _ = openimg(basename+'_chm.tif')
    except:            
        from chm import createCHM
        chm = createCHM(dtm,dsm)
        
    # create tree crowns map
    localmax = ndlocalmaxima(chm,s=2)   # get local maxima
    coos = np.argwhere(localmax == 1)
    
    Xi, Yi = np.meshgrid(xi, yi)
    peak_x, peak_y = Xi[coos[:,0],coos[:,1]],Yi[coos[:,0],coos[:,1]]
    peak_z = chm[coos[:,0],coos[:,1]]
    f = open(fn_peaks, 'w')
    for i in range(len(peak_x)):
        peak_record = '%.2f, %.2f, %.2f\n' % (peak_x[i], peak_y[i], peak_z[i])
        f.write(peak_record)
    f.close()
        
    # initiate
    crowns = np.zeros(chm.shape)*np.NaN
    old_seeds = np.zeros(chm.shape)
    
    peak = 0
    
    print 'Extract tree crowns ...'
    time.sleep(0.1)
    t0 = time.time()
    threshold = 2.5
    
    step = 0    
    for seed in coos:
        step = ETA(t0,time.time(),0.01,step,peak,0,len(coos))
        top = (peak_x[peak], peak_y[peak], peak_z[peak])
        # s_max is based on regression line with field data of Hunter et al.
        s_max = ((0.36757982565888242*top[2]+-0.87528173222693795)/2)*1.5
        seeds = np.array([seed])
        label = peak+1
        while len(seeds) > 0:  
            seed = seeds[0]
            crowns, seeds = treecrowns(chm, seed, crowns, seeds, old_seeds, label, top, s_max, Xi, Yi, threshold)
        # next tree crown
        peak+=1
        
    nodata = -999
    crowns[np.isnan(crowns)==1] = nodata    

    t1 = time.time()
    print ' --> finished in %2d seconds' % (t1-t0)
    time.sleep(0.1)    

    # peaks
    peaks_vs_trees = crowns[coos[:,0], coos[:,1]]
    
    # Get tree heights and write to file
    f = open(fn_heights, 'w')    
    for treeid in np.unique(peaks_vs_trees):
        if treeid != nodata:
            peaks_loc = np.argwhere(peaks_vs_trees == treeid) 
            coo = coos[peaks_loc,:][0][0]
            if len(coo.shape) == 1: # only one local peak in crown = tree top
                tree_x = Xi[coo[0], coo[1]]
                tree_y = Yi[coo[0], coo[1]]
                tree_z = chm[coo[0], coo[1]]
            else:   # more peaks in one crown, determine highest peak=tree top
                tree_z = chm[coo[:,0], coo[:,1]]
                highest_peak = np.argwhere(tree_z == tree_z.max())
                tree_x = Xi[coo[highest_peak,0],coo[highest_peak,1]]                
                tree_y = Yi[coo[highest_peak,0],coo[highest_peak,1]]
                tree_z = chm[coo[highest_peak,0],coo[highest_peak,1]]
            tree_record = '%.2f, %.2f, %.2f\n' % (tree_x, tree_y, tree_z)
            f.write(tree_record)
    f.close()
    print 'Tree coordinates and heights stored in', fn_heights
    
    
    # write tiff
    saveimg(np.flipud(crowns),      fn_tree, len(xi), len(yi), geotransform, proj)
    saveimg(np.flipud(crowns==-999), fn_gaps, len(xi), len(yi), geotransform, proj)
    # polygonize tree crowns
    polygonize(fn_tree, fn_shape)
    
    # plot
    import pylab as plt
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_subplot(131)
    plt.imshow(chm, extent=[Xi.min(), Xi.max(), Yi.min(), Yi.max()])
    plt.scatter(peak_x, peak_y, c = 'k', marker='+', s=100)
    plt.xlim([Xi.min(),Xi.max()])
    plt.ylim([Yi.min(),Yi.max()])
    title = 'local peaks (n=%d)' % (peak)
    plt.title(title)
    
    crowns[crowns==nodata] = np.NaN
    ax = fig.add_subplot(132)
    plt.imshow(crowns, extent=[Xi.min(), Xi.max(), Yi.min(), Yi.max()], cmap='flag')
    plt.xlim([Xi.min(),Xi.max()])
    plt.ylim([Yi.min(),Yi.max()])
    title = 'tree crowns (n=%d)' % len(np.unique(crowns[np.isnan(crowns)==False]))
    plt.title(title)
    
    ax = fig.add_subplot(133)
    plt.imshow(np.isnan(crowns), extent=[Xi.min(), Xi.max(), Yi.min(), Yi.max()])
    #plt.scatter(tree_x, tree_y, c = 'k', marker='+', s=100)
    plt.xlim([Xi.min(),Xi.max()])
    plt.ylim([Yi.min(),Yi.max()])
    plt.title('gaps')
