# -*- coding: utf-8 -*-
"""
This script calculates solar radiation for a digital elevation model (DEM)
         over one year for clear sky conditions in W/m2
-------------------------------------------------------------------
USAGE: srad = solarradiation(dem,lat,cs)
where: dem is the digital elevation model used for solar rad. estimates
       lat is the latitude vector for the DEM - same size as size(dem,2)
       cs is the cellsize in meters
       r is the ground reflectance (global value or map, default is 0.2)
      srad is the solar radiation in W/m2 over one year per grid cell

Note: Follows the approach of Kumar et al 1997. Calculates clear sky
      radiation corrected for the incident angle (selfshading) plus
      diffuse and reflected radiation. Insolation is depending on time of year (and day), 
      latitude, elevation, slope and aspect. 
      Relief shading is not considered.
      Script uses simple unweighed gradient of 4 nearest neighbours for slope
      calculation (with "gradient" function).

Reference: Kumar, L, Skidmore AK and Knowles E 1997: Modelling topographic variation in solar radiation in 
           a GIS environment. Int.J.Geogr.Info.Sys. 11(5), 475-497
Translated script of Felix Hebeler, Dept. of Geography, University Zurich,
May 2008. http://www.mathworks.com/matlabcentral/fileexchange/19791-solar-radiation

Niels Anders, EUBrazilCC 2015
"""

import math, sys, os, time
from generic_tools import getpoints, gridcellborders, saveimg
import numpy as np
import pylab as plt
from slope import calcSlope
from aspect import calcAspect

def openimg(filename):
    from osgeo import gdal
    inData = gdal.Open(filename, 0)
    #driver = gdal.GetDriverByName('HFA')
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
    
def calcLatLon(xi, yi):
    lat = np.zeros(xi.shape)
    lon = np.zeros(yi.shape)
    i = 0
    j = 0
    for la in xi: 
        for lo in yi: 
            lat[i], lon[j] = utmToLatLng(20, la, lo, northernHemisphere=False)
            j+=1
        i+=1
        j=0
    return lat, lon
    
def utmToLatLng(zone, easting, northing, northernHemisphere=True):
    if not northernHemisphere:
        northing = 10000000 - northing

    a = 6378137
    e = 0.081819191
    e1sq = 0.006739497
    k0 = 0.9996

    arc = northing / k0
    mu = arc / (a * (1 - math.pow(e, 2) / 4.0 - 3 * math.pow(e, 4) / 64.0 - 5 * math.pow(e, 6) / 256.0))

    ei = (1 - math.pow((1 - e * e), (1 / 2.0))) / (1 + math.pow((1 - e * e), (1 / 2.0)))

    ca = 3 * ei / 2 - 27 * math.pow(ei, 3) / 32.0

    cb = 21 * math.pow(ei, 2) / 16 - 55 * math.pow(ei, 4) / 32
    cc = 151 * math.pow(ei, 3) / 96
    cd = 1097 * math.pow(ei, 4) / 512
    phi1 = mu + ca * math.sin(2 * mu) + cb * math.sin(4 * mu) + cc * math.sin(6 * mu) + cd * math.sin(8 * mu)

    n0 = a / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (1 / 2.0))

    r0 = a * (1 - e * e) / math.pow((1 - math.pow((e * math.sin(phi1)), 2)), (3 / 2.0))
    fact1 = n0 * math.tan(phi1) / r0

    _a1 = 500000 - easting
    dd0 = _a1 / (n0 * k0)
    fact2 = dd0 * dd0 / 2

    t0 = math.pow(math.tan(phi1), 2)
    Q0 = e1sq * math.pow(math.cos(phi1), 2)
    fact3 = (5 + 3 * t0 + 10 * Q0 - 4 * Q0 * Q0 - 9 * e1sq) * math.pow(dd0, 4) / 24

    fact4 = (61 + 90 * t0 + 298 * Q0 + 45 * t0 * t0 - 252 * e1sq - 3 * Q0 * Q0) * math.pow(dd0, 6) / 720

    lof1 = _a1 / (n0 * k0)
    lof2 = (1 + 2 * t0 + Q0) * math.pow(dd0, 3) / 6.0
    lof3 = (5 - 2 * Q0 + 28 * t0 - 3 * math.pow(Q0, 2) + 8 * e1sq + 24 * math.pow(t0, 2)) * math.pow(dd0, 5) / 120
    _a2 = (lof1 - lof2 + lof3) / math.cos(phi1)
    _a3 = _a2 * 180 / math.pi

    latitude = 180 * (phi1 - fact1 * (fact2 + fact3 + fact4)) / math.pi

    if not northernHemisphere:
        latitude = -latitude

    longitude = ((zone > 0) and (6 * zone - 183.0) or 3.0) - _a3

    return (latitude, longitude)

def calcSR(dsm, slope_rad, aspect_rad, xi, yi, reflectance):
    [r,c] = dsm.shape
    print 'Calculate Solar radiation (%d x %d)...' % (c, r),
    time.sleep(0.1)
    t0 = time.time()

    XI,YI = np.meshgrid(xi,yi)
    
    # parameters
    n       = 1         # timestep of calculation over sunshine hours
    tau_a   = 365.      # days per year
    S0      = 1367.      # solar constant W m2
    
    lat, lon = calcLatLon(xi, yi)    
    
    L, _ = np.meshgrid(lat, lon)
    L = L*np.pi/180
    
    srad = 0
    
    term1 = np.sin(L)*np.cos(slope_rad) - np.cos(L)*np.sin(slope_rad)*np.cos(asp_rad)
    term2 = np.cos(L)*np.cos(slope_rad) + np.sin(L)*np.sin(slope_rad)*np.cos(asp_rad)
    term3 = np.sin(slope_rad)*np.sin(asp_rad)

    # loop over days per year
    I=0
    for d in np.arange(1, tau_a+1):   
        # clear sky solar radiation
        IO = S0 * (1+0.0344*np.cos(2*np.pi*d/tau_a)) # extraterr rad per day
        # sun declination
        dS = 23.45 * (np.pi/180) * np.sin(2*np.pi*((284+d)/tau_a))
        # angle at sunrise/sunset
        hsr = np.real(np.arccos(-np.tan(L)*np.tan(dS))) # only works up to 66.5 deg N
        # daylength
        It = np.round(12*(1+np.mean(hsr.flatten())/np.pi) - 12*(1-np.mean(hsr.flatten()/np.pi)))        
        I = 0
        # loop over sunshine hours per day
        for t in np.arange(1,It+1, n):

            hs = hsr - (np.pi*t/It)
            sinALpha = np.sin(L)*np.sin(dS)+np.cos(L)*np.cos(dS)*np.cos(hs)
            M = np.sqrt(1229+((614*sinALpha))**2)-614*sinALpha
            tau_b = 0.56*(np.exp(-0.65*M)+np.exp(-0.095*M))
            tau_d = 0.271-0.294*tau_b # radiation diffusion coefficient for diffuse insolation
            tau_r = 0.271+0.706*tau_b # reflectance transmitivity
            # correct for local incident angle
            cos_i = np.sin(dS)*term1 + np.cos(dS)*np.cos(hs)*term2 + np.cos(dS)*np.sin(hs)*term3
            Is = IO * tau_b # potential incoming shortwave radiation at surface normal (equator)
            # R = potential clear sky solar radiation W m2            
            R = Is * cos_i
            R[R<0] = 0      # correct negative values
            Id = IO * tau_d * (np.cos(slope_rad)*np.cos(slope_rad))/ 2*sinALpha # diffuse radiation
            Ir = IO * reflectance * tau_r * (np.sin(slope_rad)*np.sin(slope_rad))/2*sinALpha # reflectance
            R = R+ Id + Ir
            R[R<0] = 0      # correct negative values
            I=I+R           # solar radiation per ay (sunshine hours)
        # add up radiation
        srad = srad+I
    
    # return
    t1 = time.time()
    print 'finished in %1d seconds' % (t1-t0)
    time.sleep(0.1)     
    
    return srad
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
    fn  = basename+'_SolarRadiation.tif' 
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
    
    # create dsm
    from dsm import createDSM
    dsm = createDSM(x,y,z,xi,yi,bx,by,res,geotransform, proj)
    saveimg(np.flipud(dsm), basename+"_dsm.tif", len(xi), len(yi), geotransform, proj)
   
    # calculate slope and aspect in radians
    slope_rad = np.core.umath.deg2rad(calcSlope(dsm))
    asp_rad = np.core.umath.deg2rad(calcAspect(dsm))
    
    slope_rad = np.flipud(slope_rad)    
    asp_rad = np.flipud(asp_rad)
    
    # calculate solar radiation
    SR = calcSR(dsm, slope_rad, asp_rad, xi, yi, 0.2)
    
    # write tiff
    saveimg(SR, fn, len(xi), len(yi), geotransform, proj)

    plt.figure(figsize=(15,5))
    plt.imshow(asp_rad)
    plt.colorbar(shrink=0.5)
    plt.title('aspect [rad]')
        
    plt.figure(figsize=(15,5))
    plt.imshow(slope_rad)
    plt.colorbar(shrink=0.5)    
    plt.title('slope angle [rad]')

    plt.figure(figsize=(15,5))
    SR[SR==-999] = np.NaN
    plt.imshow(SR)
    plt.colorbar(shrink=0.5)
    plt.title('solar radiation')
    

