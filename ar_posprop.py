# ;Edited by S Murray to get boundaries of blobs for Hexa (2015-07-30)
# ;Input a processed (pre-run mag cos-corr) data map and mask
# ;The tot. area, pos. area, neg. area, and
# ;	total, signed, fractional signed, negative, and positive flux
# ;	are determined
# ;The returned structure contains position info for the whole
# ;detection.
# ;OUTPOSSTR has position information for the positive pixels in the
# ;detection
# ;OUTNEGSTR has position '' for negative ''
# ;
# ;Status:
# ;	0 = initialised value
# ;	7 = All desired positions should have been calculated
# ;	1 = No positive pixels found; positive positions not calulated/output
# ;	2 = No negative pixels found; negative positions not calulated/output
# ;	3 = No positive or nagative pixels found; Shouldn't happen!! (expect crash in this case)
# ;	4 = No detection 'where' values found for detection mask; Shouldn't happen!!
# ;------------------------------------------------------------------------------>
#
# ;Determine the AR positions, given mask indices
# ;Input/Output:
# ;   inoutstr = a blank structure that will be filled
# ;Input:
# ;   WAR = an array of the AR pixel indices
# ;   DATA = the masked abs(data) array (cosine-area and -mag corrected; in flux units)
# ;

import numpy as np
from ar_detect import xyrcoord, ar_pxscale
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import pandas as pd

def ar_posprop(thismap, thismask, cosmap):
    """
    """
#    pxmmsq =  ar_pxscale(thismap, cmsqr=False, mmppx=False, cmppx=False)
#    pxcmsq = ar_pxscale(thismap, cmsqr=True, mmppx=False, cmppx=False)
    nmask = np.max(thismask)
    # Create dataframes
    ardf = pd.DataFrame(columns = ['arid',
                                   'xminbnd', 'yminbnd', 'xmaxbnd', 'ymaxbnd', 'xcenbnd', 'ycenbnd',
                                   'xcenflx', 'ycenflx', 'xcenarea', 'ycenarea',
                                   'hcxbnd', 'hcybnd', 'hcxflx', 'hcyflx', 'hcxarea', 'hcyarea',
                                   'hglonbnd', 'hglatbnd', 'hglonflx', 'hglatflx', 'hglonarea', 'hglatarea',
                                   'carlonbnd', 'carlonflx', 'carlonarea'])
    posdf = ardf.copy()
    negdf = ardf.copy()
    #For each AR...
    for i in range(1, np.int(nmask)+1):
        # Zero pixels outside detection boundary
        #tmpmask = np.copy(thismask)
        #tmpmask[np.where(thismask != i)] = 0.
        tmpdat = np.copy(thismap.data)
        tmpdat[np.where(thismask != i)] = 0.
        tmpabs = np.abs(tmpdat)
        tmpflx = tmpabs*cosmap
        # Where are values within the detection boundary
        #tmpmask[np.where(thismask == i)] = 1.
#        nothresh=0 & nopos=0 & noneg=0 & noposbnd=0 & nonegbnd=0
        # Where are the signed values
#        wneg = np.where(tmpdat < 0.)
#        wpos = np.where(tmpdat > 0.)
        # Fill the position structure for whole detection
        arstr = ar_posprop_findpos(i, thismap, np.where(thismask == i), tmpflx)
        # Same for positive and negative regions
        posstr = ar_posprop_findpos(i, thismap, np.where(tmpdat > 0), tmpflx)
        negstr = ar_posprop_findpos(i, thismap, np.where(tmpdat < 0), tmpflx)
        # Fill position dataframes
        ardf = ardf.append([arstr])
        posdf = posdf.append([posstr])
        negdf = negdf.append([negstr])
    #TO DO: ?? check differences in lists between posprop and findpos??
    return ardf#, posdf, negdf

def ar_posprop_findpos(arid, thismap, war, data):
    """
    """
    data = np.abs(data)
    date = thismap.meta['t_obs']
    dx = thismap.meta['cdelt1']
    dy = thismap.meta['cdelt2']
    rsun = thismap.meta["rsun_obs"]
    xc = thismap.meta["crval1"] + thismap.meta["cdelt1"]*(((thismap.meta["naxis1"] + 1)/2) - thismap.meta["crpix1"])
    #xc = thismap.meta["crval1"] + thismap.meta["cdelt1"]*np.cos(thismap.meta["crota2"])*(((thismap.meta["naxis1"] + 1)/2) - thismap.meta["crpix1"]) - thismap.meta["cdelt2"]*np.sin(thismap.meta["crota2"])*(((thismap.meta["naxis2"] + 1)/2) - thismap.meta["crpix2"])
    yc = thismap.meta["crval2"] + thismap.meta["cdelt2"] * (((thismap.meta["naxis2"] + 1) / 2) - thismap.meta["crpix2"])
    #yc = thismap.meta["crval2"] + thismap.meta["cdelt1"]*np.sin(thismap.meta["crota2"])*(((thismap.meta["naxis1"] + 1)/2) - thismap.meta["crpix1"]) + thismap.meta["cdelt2"]*np.cos(thismap.meta["crota2"])*(((thismap.meta["naxis2"] + 1)/2) - thismap.meta["crpix2"])
    # get x and y indices
    sz = thismap.data.shape
    xx, yy, rr = xyrcoord(sz)
    xwar = war[1]
    ywar = war[0]
    # Determine the bounding box center positions
    xminmax = np.array((np.min(xwar), np.max(xwar)))
    yminmax = np.array((np.min(ywar), np.max(ywar)))
    xcenbnd = np.mean(xminmax) #in pixels from lower left corner of the image FOV
    ycenbnd = np.mean(yminmax) #in pixels
    xminbnd = xminmax[0] #boundary coords
    yminbnd = yminmax[0]
    xmaxbnd = xminmax[1]
    ymaxbnd = yminmax[1]
    hcxbnd, hcybnd = px2hc(xcenbnd, ycenbnd, dx, dy, xc, yc, sz)
    #hgxbnd, hgybnd, carxbnd = hc2hg(hcxbnd, hcybnd, date, rsun, carx=True)
    hgxbnd, hgybnd, carxbnd = hc2hg(thismap, hcxbnd, hcybnd)
    # Determine the area weighted centroids
    xcenarea = np.sum(xwar * xx[war]) / np.sum(xx[war])
    ycenarea = np.sum(ywar * yy[war]) / np.sum(yy[war])
    hcxarea, hcyarea = px2hc(xcenarea, ycenarea, dx, dy, xc, yc, sz)
    hgxarea, hgyarea, carxarea = hc2hg(thismap, hcxarea, hcyarea)
    # Determine the Flux weighted centroids
    xcenflx = np.sum(xwar * data[war]) / np.sum(data[war])
    ycenflx = np.sum(ywar * data[war]) / np.sum(data[war])
    hcxflx, hcyflx = px2hc(xcenflx, ycenflx, dx, dy, xc, yc, sz)
    hgxflx, hgyflx, carxflx = hc2hg(thismap, hcxflx, hcyflx)
    # Fill structure
    posstr = {'arid': arid,
              'xminbnd': xminbnd, 'yminbnd': yminbnd,
              'xmaxbnd': xmaxbnd, 'ymaxbnd': ymaxbnd,
              'xcenbnd': xcenbnd, 'ycenbnd': ycenbnd,
              'xcenflx': xcenflx, 'ycenflx': ycenflx,
              'xcenarea': xcenarea, 'ycenarea': ycenarea,
              'hcxbnd': hcxbnd, 'hcybnd': hcybnd,
              'hcxflx': hcxflx, 'hcyflx': hcyflx,
              'hcxarea': hcxarea, 'hcyarea': hcyarea,
              'hglonbnd': hgxbnd, 'hglatbnd': hgybnd,
              'hglonflx': hgxflx, 'hglatflx': hgyflx,
              'hglonarea': hgxarea, 'hglatarea': hgyarea,
              'carlonbnd': carxbnd, 'carlonflx': carxflx, 'carlonarea': carxarea}
    return posstr


def px2hc(xpx, ypx, dx, dy, xc, yc, sz):
    """Given positions in pixel coordinates (x = [0->nx-1]; y = [0->ny-1]),
    determine the HC x and y coordinates in arcseconds.
    xpx, ypx = pixel coordinates of AR positions.
    dx, dy = arcsec/px
    xc, yc = FOV center in arcsec from solar disk center
    xs, ys = x and y px size of the image
    """
    xs = sz[0]
    ys = sz[1]
    hcx = ((xpx - (xs/2.))*dx) + xc
    hcy = ((ypx - (ys/2.))*dy) + yc
    return hcx, hcy

def hc2hg(thismap, hcxbnd, hcybnd):
    """Convert heliocentric to heliographic coordinates using sunpy coordinates
    hc coords are in arcsecs
    hg coords are in degrees
    Also outputs Carrington longitude in degrees
    """
    heliocentric = SkyCoord(hcxbnd, hcybnd, unit='arcsec', frame=thismap.coordinate_frame)
    heliographic = heliocentric.transform_to(frames.HeliographicStonyhurst)
    carrington = heliographic.transform_to(frames.HeliographicCarrington)
    return heliographic.lon.value, heliographic.lat.value, carrington.lon.value


def old_hc2hg(hcxbnd, hcybnd, date, rsun, carx):
    """Convert heliocentric to heliographic coordinates using the usual SSW routines.
    hc coords are in arcsecs
    hg coords are in degrees
    rsun is in arcsec
    Optionally outputs the Carrington longitude
    """
    return


