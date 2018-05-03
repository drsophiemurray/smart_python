'''
    SMART position property code
    ============================
    Written by Sophie A. Murray, code originally developed by Paul Higgins (ar_posprop.pro).

    Developed under Python 3 and Sunpy 0.8.3
    - Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)

    Provides position information of detected SMART regions in various coordinate forms.

    Inputs:
    - inmap: Processed magnetogram
    - inmask: Output SMART mask from ar_detect_core

'''

import numpy as np
from detect import xyrcoord, ar_pxscale
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import pandas as pd

def main(inmap, inmask, cosmap):
    """
    Determine the AR positions, given mask indices.
    """
#    pxmmsq = ar_pxscale(inmap, cmsqr=False, mmppx=False, cmppx=False)
#    pxcmsq = ar_pxscale(inmap, cmsqr=True, mmppx=False, cmppx=False)
    nmask = np.max(inmask)
    ## Create dataframes to output (total, positive, negative)
    ardf = pd.DataFrame(columns = ['arid',
                                   'xminbnd', 'yminbnd', 'xmaxbnd', 'ymaxbnd', 'xcenbnd', 'ycenbnd',
                                   'xcenflx', 'ycenflx', 'xcenarea', 'ycenarea',
                                   'hcxbnd', 'hcybnd', 'hcxflx', 'hcyflx', 'hcxarea', 'hcyarea',
                                   'hglonbnd', 'hglatbnd', 'hglonflx', 'hglatflx', 'hglonarea', 'hglatarea',
                                   'carlonbnd', 'carlonflx', 'carlonarea'])
    posdf = ardf.copy()
    negdf = ardf.copy()
    ## For each AR get position information
    for i in range(1, np.int(nmask)+1):
        ## Zero pixels outside detection boundary
#        tmpmask = np.copy(inmask)
#        tmpmask[np.where(inmask != i)] = 0.
        tmpdat = np.copy(inmap.data)
        tmpdat[np.where(inmask != i)] = 0.
        tmpabs = np.abs(tmpdat)
        tmpflx = tmpabs*cosmap
        ## Where are values within the detection boundary
#        tmpmask[np.where(inmask == i)] = 1.
        ## Where are the signed values
#        wneg = np.where(tmpdat < 0.)
#        wpos = np.where(tmpdat > 0.)
        # Fill the position structure for whole detection
        arstr = ar_posprop_findpos(i, inmap, np.where(inmask == i), tmpflx)
        # Same for positive and negative regions
        posstr = ar_posprop_findpos(i, inmap, np.where(tmpdat > 0), tmpflx)
        negstr = ar_posprop_findpos(i, inmap, np.where(tmpdat < 0), tmpflx)
        # Fill position dataframes
        ardf = ardf.append([arstr], ignore_index=True)
        posdf = posdf.append([posstr], ignore_index=True)
        negdf = negdf.append([negstr], ignore_index=True)
    #TO DO: check differences in lists between posprop and findpos??
    ## Currently only outputting total region positional information...
    return ardf#, posdf, negdf

def ar_posprop_findpos(arid, inmap, war, data):
    """
    Find positional information in diffferent coordinate frames.
    """
    data = np.abs(data)
    date = inmap.meta['t_obs']
    dx = inmap.meta['cdelt1']
    dy = inmap.meta['cdelt2']
    rsun = inmap.meta["rsun_obs"]
    xc = inmap.meta["crval1"] + inmap.meta["cdelt1"]*(((inmap.meta["naxis1"] + 1)/2) - inmap.meta["crpix1"])
#    xc = inmap.meta["crval1"] + inmap.meta["cdelt1"]*np.cos(inmap.meta["crota2"])*(((inmap.meta["naxis1"] + 1)/2) - inmap.meta["crpix1"]) - inmap.meta["cdelt2"]*np.sin(inmap.meta["crota2"])*(((inmap.meta["naxis2"] + 1)/2) - inmap.meta["crpix2"])
    yc = inmap.meta["crval2"] + inmap.meta["cdelt2"] * (((inmap.meta["naxis2"] + 1) / 2) - inmap.meta["crpix2"])
#    yc = inmap.meta["crval2"] + inmap.meta["cdelt1"]*np.sin(inmap.meta["crota2"])*(((inmap.meta["naxis1"] + 1)/2) - inmap.meta["crpix1"]) + inmap.meta["cdelt2"]*np.cos(inmap.meta["crota2"])*(((inmap.meta["naxis2"] + 1)/2) - inmap.meta["crpix2"])
    ## Get x and y indices
    sz = inmap.data.shape
    xx, yy, rr = xyrcoord(sz)
    xwar = war[1]
    ywar = war[0]
    ## Determine the bounding box center positions
    xminmax = np.array((np.min(xwar), np.max(xwar)))
    yminmax = np.array((np.min(ywar), np.max(ywar)))
    xcenbnd = np.mean(xminmax) #in pixels from lower left corner of the image FOV
    ycenbnd = np.mean(yminmax) #in pixels
    xminbnd = xminmax[0] #boundary coords
    yminbnd = yminmax[0]
    xmaxbnd = xminmax[1]
    ymaxbnd = yminmax[1]
    hcxbnd, hcybnd = px2hc(xcenbnd, ycenbnd, dx, dy, xc, yc, sz)
    hgxbnd, hgybnd, carxbnd = hc2hg(inmap, hcxbnd, hcybnd) #old version used cxbnd, hcybnd, date, rsun, carx=True
    ## Determine the area weighted centroids
    xcenarea = np.sum(xwar * xx[war]) / np.sum(xx[war])
    ycenarea = np.sum(ywar * yy[war]) / np.sum(yy[war])
    hcxarea, hcyarea = px2hc(xcenarea, ycenarea, dx, dy, xc, yc, sz)
    hgxarea, hgyarea, carxarea = hc2hg(inmap, hcxarea, hcyarea)
    ## Determine the Flux weighted centroids
    xcenflx = np.sum(xwar * data[war]) / np.sum(data[war])
    ycenflx = np.sum(ywar * data[war]) / np.sum(data[war])
    hcxflx, hcyflx = px2hc(xcenflx, ycenflx, dx, dy, xc, yc, sz)
    hgxflx, hgyflx, carxflx = hc2hg(inmap, hcxflx, hcyflx)
    ## Fill structure
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
    """
    Given positions in pixel coordinates (x = [0->nx-1]; y = [0->ny-1]),
    determine the HC x and y coordinates in arcseconds.
    - xpx, ypx = pixel coordinates of AR positions
    - dx, dy = arcsec/px
    - xc, yc = FOV center in arcsec from solar disk center
    - xs, ys = x and y px size of the image
    """
    xs = sz[0]
    ys = sz[1]
    hcx = ((xpx - (xs/2.))*dx) + xc
    hcy = ((ypx - (ys/2.))*dy) + yc
    return hcx, hcy

def hc2hg(inmap, hcxbnd, hcybnd):
    """
    Convert heliocentric to heliographic coordinates using sunpy coordinates.
    hc coords are in arcsecs
    hg coords are in degrees
    Also outputs Carrington longitude in degrees
    """
    heliocentric = SkyCoord(hcxbnd, hcybnd, unit='arcsec', frame=inmap.coordinate_frame)
    heliographic = heliocentric.transform_to(frames.HeliographicStonyhurst)
    carrington = heliographic.transform_to(frames.HeliographicCarrington)
    return heliographic.lon.value, heliographic.lat.value, carrington.lon.value

if __name__ == '__main__':
    main()