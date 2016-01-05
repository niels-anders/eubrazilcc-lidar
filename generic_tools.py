# -*- coding: utf-8 -*-
"""
Generic tools for lidar processing

(c) 2015 EUBrazilCC, Niels Anders
"""
import time, sys
import numpy as np
from osgeo import gdal

def ETA(t0,t1, s, step, i, r, k):
    """
    t0: start start
    t1: current time
    s: step size to print to screen
    step: current step
    i: current row
    r: start row (mostly 0)
    k: total number of rows    
    """
    pc = float(i+1- r)/k
    if pc >= s:
        t1 = time.time()
        time_running = t1-t0
        time_left = time_running/pc-time_running
        ETA = t1+time_left
        sys.stdout.write( '\r%3.1f' % (pc*100) +' % completed --> ETA: '+ time.strftime('%d %b %Y %H:%M:%S',time.localtime(ETA)))
        sys.stdout.flush()
        s = s+ step
    return s

def gridcellborders(xi,yi,res):
    # grid cell borders
    bx = np.arange(xi.min() - res/float(2), xi.max() + res, res)
    by = np.arange(yi.min() - res/float(2), yi.max() + res, res)
    return bx, by
    
def saveimg(Array, filename, cols, rows, geotransform, proj, printYN = 'no'):
    # change NaN to -999 nodata value
    Array[np.isnan(Array)] = -999
    # get driver
    driver = gdal.GetDriverByName('GTiff')
    # create file
    outData = driver.Create(filename, cols, rows, 1, gdal.GDT_Float32)
    # get band to write to
    outBand = outData.GetRasterBand(1)
    # write array
    outBand.WriteArray(Array, 0, 0)
    outBand.FlushCache()
    # set nodata
    outBand.SetNoDataValue(-999)
    # calculate statistics
    outBand.GetStatistics(0,1)
    # georeference
    if geotransform != "":    
        outData.SetGeoTransform(geotransform)
        outData.SetProjection(proj)

    outData = None
    if printYN == 'yes': print 'File saved as '+filename

def openimg(filename):
    inData = gdal.Open(filename, 0)
    if inData is None:
        print 'Could not open ' + filename

    # get metadata
    cols = inData.RasterXSize
    rows = inData.RasterYSize
    # geotransform = (left, cellsize, 0.0, top, 0.0, -cellsize)
    geotransform = inData.GetGeoTransform()
    proj = inData.GetProjection()
    res = geotransform[1]
    # get band
    band = inData.GetRasterBand(1)
    # get data as array
    data = band.ReadAsArray(0,0, cols, rows)
    return data, cols, rows, res, geotransform, proj
    
def idw(x,y,z,xi,yi, nnear=501):
    import invdisttree
    xy = np.append([x],[y],axis=0).transpose() 
    xiyi = np.append([xi],[yi], axis=0).transpose()
    idw = invdisttree.Invdisttree(xy, z, leafsize=10, stat=1)
    weights = np.ones(z.shape)
    zi = idw(xiyi, nnear=nnear, eps=0.5, p=2 , weights=weights)
    #zi = zi.reshape(xi.shape)
    return zi
    
def getpoints(In, printYN='no'):
    if printYN=='yes':
        print 'get points....',
        time.sleep(0.1)
        t0 = time.time()

    if (In[-3:] == 'laz') or (In[-3:] == 'las'):
        import liblas
        m = liblas.file.File(In, mode='r')           
        nrpoints = int(m.header.point_records_count)
        data = np.zeros((nrpoints,4)) 
        x = np.zeros((nrpoints))
        y = np.zeros((nrpoints))
        z = np.zeros((nrpoints))
        c = np.zeros((nrpoints))
        i = 0
        for p in m:
            x[i] = p.x
            y[i] = p.y
            z[i] = p.z
            c[i] = p.classification
            i+=1 
        
    if In[-3:] == 'txt':
        f = open(In)
        lines = f.readlines()
        data = np.zeros((len(lines),6))
        i = 0
        if printYN=='yes':
            t0 = time.time()
        for line in lines:
            line = np.array(filter(None, line.split(' ')))
            data[i] = line.astype('float64')      
            i = i+1      
        x,y,z,c = data[:,0], data[:,1], data[:,2], data[:,5]   
    
    # return
    if printYN=='yes':
        t1 = time.time()
        print 'finished in %2d seconds' % (np.round(t1-t0))
        time.sleep(0.1)  
    return x,y,z,c
