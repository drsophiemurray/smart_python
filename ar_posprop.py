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


def ar_posprop(map, mask, cosmap):
    """	blankstr={datafile:indatafile[0],arid:0,xcenbnd:0d, ycenbnd:0d, xcenflx:0d, ycenflx:0d, xcenarea:0d, ycenarea:0d, $
          hcxbnd:0d, hcybnd:0d, hcxflx:0d, hcyflx:0d, hcxarea:0d, hcyarea:0d, $
          hglonbnd:0d, hglatbnd:0d, hglonflx:0d, hglatflx:0d, hglonarea:0d, hglatarea:0d, $
          carlonbnd:0d,  carlonflx:0d,  carlonarea:0d, $
          xminbnd:0d, yminbnd:0d, xmaxbnd:0d, ymaxbnd:0d}
    """
#    pxmmsq =  ar_pxscale(map, cmsqr=False, mmppx=False, cmppx=False)
#    pxcmsq = ar_pxscale(map, cmsqr=True, mmppx=False, cmppx=False)
    nmask = np.max(mask)
    #For each AR...
    for i in range(1, np.int(nmask)+1):
        # Zero pixels outside detection boundary
        thismask = np.copy(mask)
        thismask[np.where(mask != i)] = 0.
        thisdat = np.copy(map.data)
        thisdat[np.where(mask != i)] = 0.
        thisabs = np.abs(thisdat)
        thisflx = thisabs*cosmap
        # Where are values within the detection boundary
        thismask[np.where(mask == i)] = 1.
#        nothresh=0 & nopos=0 & noneg=0 & noposbnd=0 & nonegbnd=0
        # Where are the signed values
#        wneg = np.where(thisdat < 0.)
#        wpos = np.where(thisdat > 0.)
        # Same for positive and negative regions
        posstr? = ar_posprop_findpos(map, np.where(thisdat > 0), thisflx)
        negstr? = ar_posprop_findpos(map, np.where(thisdat < 0), thisflx)
        # Fill the position structure for - whole detection
        arstr? = ar_posprop_findpos(map, np.where(mask == i), thisflx)
    # Fill position structure
    ?? check differences in lists between posprop and findpos??
    return arstr, posstr, negstr

def ar_posprop_findpos(map, war, data):
    """
    """
    data = np.abs(data)
    date = map.meta['t_obs']
    dx = map.meta['cdelt1']
    dy = map.meta['cdelt2']
    rsun = map.meta["rsun_obs"]
    xc = map.meta["crval1"] + map.meta["cdelt1"]*(((map.meta["naxis1"] + 1)/2) - map.meta["crpix1"])
    #xc = map.meta["crval1"] + map.meta["cdelt1"]*np.cos(map.meta["crota2"])*(((map.meta["naxis1"] + 1)/2) - map.meta["crpix1"]) - map.meta["cdelt2"]*np.sin(map.meta["crota2"])*(((map.meta["naxis2"] + 1)/2) - map.meta["crpix2"])
    yc = map.meta["crval2"] + map.meta["cdelt2"] * (((map.meta["naxis2"] + 1) / 2) - map.meta["crpix2"])
    #yc = map.meta["crval2"] + map.meta["cdelt1"]*np.sin(map.meta["crota2"])*(((map.meta["naxis1"] + 1)/2) - map.meta["crpix1"]) + map.meta["cdelt2"]*np.cos(map.meta["crota2"])*(((map.meta["naxis2"] + 1)/2) - map.meta["crpix2"])
    # get x and y indices
    sz = map.data.shape
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
    hgxbnd, hgybnd, carxbnd = hc2hg(map, hcxbnd, hcybnd)
    # Determine the area weighted centroids
    xcenarea = np.sum(xwar * xx[war]) / np.sum(xx[war])
    ycenarea = np.sum(ywar * yy[war]) / np.sum(yy[war])
    hcxarea, hcyarea = px2hc(xcenarea, ycenarea, dx, dy, xc, yc, sz)
    hgxarea, hgyarea, carxarea = hc2hg(map, hcxarea, hcyarea)
    # Determine the Flux weighted centroids
    xcenflx = np.sum(xwar * data[war]) / np.sum(data[war])
    ycenflx = np.sum(ywar * data[war]) / np.sum(data[war])
    hcxflx, hcyflx = px2hc(xcenflx, ycenflx, dx, dy, xc, yc, sz)
    hgxflx, hgyflx, carxflx = hc2hg(map, hcxflx, hcyflx)
# ;Fill structure
# inoutstr.xminbnd=xminbnd
# inoutstr.yminbnd=yminbnd
# inoutstr.xmaxbnd=xmaxbnd
# inoutstr.ymaxbnd=ymaxbnd
# inoutstr.xcenbnd=xcenbnd
# inoutstr.ycenbnd=ycenbnd
# inoutstr.xcenflx=xcenflx
# inoutstr.ycenflx=ycenflx
# inoutstr.xcenarea=xcenarea
# inoutstr.ycenarea=ycenarea
# inoutstr.hcxbnd=hcxbnd
# inoutstr.hcybnd=hcybnd
# inoutstr.hcxflx=hcxflx
# inoutstr.hcyflx=hcyflx
# inoutstr.hcxarea=hcxarea
# inoutstr.hcyarea=hcyarea
# inoutstr.hglonbnd=hgxbnd
# inoutstr.hglatbnd=hgybnd
# inoutstr.hglonflx=hgxflx
# inoutstr.hglatflx=hgyflx
# inoutstr.hglonarea=hgxarea
# inoutstr.hglatarea=hgyarea
# inoutstr.carlonbnd=carxbnd
# inoutstr.carlonflx=carxflx
# inoutstr.carlonarea=carxarea
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

def hc2hg(map, hcxbnd, hcybnd):
    """Convert heliocentric to heliographic coordinates using sunpy coordinates
    hc coords are in arcsecs
    hg coords are in degrees
    Also outputs Carrington longitude in degrees
    """
    heliocentric = SkyCoord(hcxbnd, hcybnd, unit='arcsec', frame=map.coordinate_frame)
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


