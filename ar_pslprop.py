# ;Determine many PSL properties of an AR
# ;
# ;INPUTS:
# ;	inmap = data map of magnetogram
# ;	inmask = indexed image mask of detected features
# ;
# ;OPT. INPUTS:
# ;	refimg = an image to be projected in the same manner as the data and mask
# ;
# ;KEYWORDS:
# ;	param = The SMART2 parameter structure
# ;	fparam = Filename of desired parameter structure to use (PARAM takes precedence)
# ;	doproj = set to do a stereographic deprojection when determining PSL props.
# ;	projscl = DEFUNCT choose the factor to increase the projection image dimensions by, as compared to the original image
# ;	dobppxscl = set to use the bipole separation length to determine the deprojected pixel scaling (= bp_sep/bp_sep_proj)
# ;	outproj = ouput the STG projected magnetogram
# ;	outscl = output the 'true' pixel scaling for each projected pixel using the gradient of Rdeg map
# ;	outbpscl = output the conversion between great-circle and STG projected bipole separation length to use as conversion between scaling
# ;	outmaskproj = STG projected mask
# ;	outrefproj = STG projected version of input REFIMG
# ;	projlimbxy = x,y positions of the limb in STG projected space
# ;	outpslmask = mask of PSL in STG space
# ;	outgradpsl = gradient image of PSL in STG space
# ;
# ;NOTES:
# ;	1. If /DOPROJ is NOT set, then all outputs labeled 'STG projected' will just be in LOS (HCP) space, corresponding to the input map.

projscl = 1.5 # F; the fractional increase/decrease in image size for the projected image
mdi_noisethresh= 70.0 #; F; a noise threshold used in processing MDI magnetograms


import numpy as np
import sunpy.map
from ar_detect import xyrcoord, ar_grow
from ar_posprop import px2hc, hc2hg
from sunpy.sun import constants
import pandas as pd
import scipy.interpolate
import scipy.ndimage
import astropy.units as u

def ar_pslprop(map, thismask, dproj, projmaxscale):
    """.arpslstr={arid:0,psllength:0.,pslsglength:0.,pslcurvature:0d,rvalue:0.,wlsg:0.,bipolesep_mm:0.,bipolesep_px:0.,bipolesep_proj:0.}
    """
    sz = map.data.shape
    xscale = map.meta['cdelt1']
    nmask = np.max(mask)
    for i in range(1, np.int(nmask)+1):
        # Zero pixels outside detection boundary
        thismask = np.copy(mask)
        thismask[np.where(mask != i)] = 0.
        thismask[np.where(mask == i)] = 1.
        thisdat = map.data*thismask
        # Take a sub-map around the AR
        thisdatmap = sunpy.map.Map(thisdat, map.meta)
        maskmap = sunpy.map.Map(thismask, map.meta)
        #xrange = ((np.min(np.where(thismask == 1)[1])-1, np.max(np.where(thismask == 1)[1])+1))
        #yrange = ((np.min(np.where(thismask == 1)[0])-1, np.max(np.where(thismask == 1)[0])+1))
        bottom_left_pixels = ((np.min(np.where(thismask == 1)[1])-1, np.min(np.where(thismask == 1)[0])-1))
        #bl = pix_to_arc(map, bottom_left_pixels[0], bottom_left_pixels[1])
        top_right_pixels = ((np.max(np.where(thismask == 1)[1])+1, np.max(np.where(thismask == 1)[0])+1))
        #tr = pix_to_arc(map, top_right_pixels[0], top_right_pixels[1])
        submask = maskmap.submap(bottom_left_pixels * u.pixel, top_right_pixels * u.pixel)
        submag = thisdatmap.submap(bottom_left_pixels * u.pixel, top_right_pixels * u.pixel)
        #Convert to wcs structure??
        # Determine the bipole separation properties
        bipsepstr = ar_bipolesep(submag)
        # Stereographic deprojection --- TO DO!!!!! DOPROJ=1
        if doproj is False:
            projpxscl = np.ones(sz)
            projmag = np.copy(submag.data)
            rim = np.copy(submag.data)
            projmask = np.copy(submask.data)
            bisepstrproj = bipsepstr
            projpxscl_bpsep = 1.
        projsz = projmag.shape
        # Choose whether to use the Rdeg gradient or the bipole separation conversion to determine the projected pixel scaling
        if dobppxscl is True:
            projmmscl = ar_pxscale(submag, cmsqr=False, mmppx=True, cmppx=False) * projpxscl_bpsep
        else:
            projmmscl = ar_pxscale(submag, cmsqr=False, mmppx=True, cmppx=False) * projpxscl
        kernpsl = [[0, 0, 1, 1, 1, 0, 0], [0, 1, 1, 1, 1, 1, 0], [1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1],
                   [1, 1, 1, 1, 1, 1, 1], [0, 1, 1, 1, 1, 1, 0], [0, 0, 1, 1, 1, 0, 0]]
        kernpsl = np.array(kernpsl)
        kernsz = kernpsl.shape
        # Resize the kernel based on the scale conversion
        if (np.min(kernsz[0]/projpxscl_bpsep) < 1) or (np.isnan(np.min(kernsz[0]/projpxscl_bpsep)) is True)):
            kernpsl = rebin(kernpsl, (kernsz[0], kernsz[1]))
        else:
            factor = (kernsz[0]/projpxscl_bpsep)/kernsz[0]
            kernpsl = scipy.ndimage.zoom(kernpsl, factor) #need to get congrid working but this will do for now
        projmagg = ar_grow(projmag, 1, gauss = True)
        psz = projmagg.shape
        nmask = np.zeros(projmagg.shape)
        pmask = np.copy(nmask)
        nmask[np.where(projmagg < (-mdi_noisethresh*2))] = 1.
        pmask[np.where(projmagg > mdi_noisethresh*2)] = 1.
        pmaskg = ar_grow(pmask, 1./2., gauss=False, kern=kernpsl)
        nmaskg = ar_grow(nmask, 1./2., gauss=False, kern=kernpsl)
        pslmask = np.zeros(projmagg.shape)
        pslmask[np.where(pmaskg + nmaskg == 2)] = 1.


def ar_bipolesep(map):
    """Determine the flux-weighted bipole separation distance between the pos and neg centroids
    ;in degrees if a map is input, or in Px if only an image is input
    """
    image = np.copy(map.data)
    imgsz = image.shape
    yy, xx, rr = xyrcoord(imgsz) #need to check this shit theres definitely something wrong - array flipped
    imageneg = np.copy(image)
    imageneg[np.where(image > 0.)] = 0.
    imagepos = np.copy(image)
    imagepos[np.where(image < 0.)] = 0.
    pxpxloc = np.sum(xx * imagepos) / np.sum(imagepos)
    nxpxloc = np.sum(xx * np.abs(imageneg)) / np.sum(np.abs(imageneg))
    pypxloc = np.sum(yy * imagepos) / np.sum(imagepos)
    nypxloc = np.sum(yy * np.abs(imageneg)) / np.sum(np.abs(imageneg))
    pxsep = np.sqrt((pxpxloc-nxpxloc)**2.+(pypxloc-nypxloc)**2.)
    # Now the map outputs
    xc = map.meta["crval1"] + map.meta["cdelt1"]*(((map.meta["naxis1"] + 1) / 2) - map.meta["crpix1"])
    yc = map.meta["crval2"] + map.meta["cdelt2"] * (((map.meta["naxis2"] + 1) / 2) - map.meta["crpix2"])
    phcxflx, phcyflx = px2hc(pxpxloc, pypxloc, map.meta['cdelt1'], map.meta['cdelt2'], xc, yc, imgsz[::-1]) #again flipped wtf
    phgxflx, phgyflx, carpxflx = hc2hg(map, phcxflx, phcyflx)
    nhcxflx, nhcyflx = px2hc(nxpxloc, nypxloc, map.meta['cdelt1'], map.meta['cdelt2'], xc, yc, imgsz[::-1]) #again flipped wtf
    nhgxflx, nhgyflx, carnxflx = hc2hg(map, nhcxflx, nhcyflx)
    gcdist_deg, outeqnode,  outadist = gc_dist(np.array((phgxflx, phgyflx)), np.array((nhgxflx, nhgyflx)), nonan=False)
    gcdist_mm = (gcdist_deg / 360.) * 2. * np.pi * (constants.radius.value / 1e6)
    gcdist_px = (gcdist_deg / 360.) * 2. * np.pi * (submag.meta['rsun_obs'] / submag.meta['cdelt1'])
    sepstr = {'pxcen': pxpxloc, 'pycen': pypxloc, 'nxcen': nxpxloc, 'nycen': nypxloc,
              'plon': phgxflx, 'plat': phgyflx, 'nlon': nhgxflx, 'nlat': nhgyflx, 'pxsep': pxsep,
              'gdist_deg': gcdist_deg, 'gcdist_mm': gcdist_mm, 'gcdist_px': gcdist_px}
    return sepstr

def gc_dist(alonlat, blonlat, nonan):
    """
    Return the distance along the great circle between two reference points from Leonard (1953)
    Coordinates are in geographic longitude latitude
    alonlat: first reference point on great circle (GC), [longitude, latitude] in degrees
    blonlat: second reference point on great circle (GC), [longitude, latitude] in degrees
    OutEqNode: output the location of the nearest equatorial node to the reference point in degrees
    OutADist: distance between reference point A and the equatorial node of the GC. in degrees
    NoNaN: where ALONLAT is equal to BLONLAT GC_DIST will return a NaN. To replace the NaNs with 0, set /NoNaN
    """
    # Bit of a fudge apparently - looks like he doesnt want zeros
    alonlat[np.where(alonlat == 0.)] = 0.001
    blonlat[np.where(blonlat == 0.)] = 0.001
    # Inclination between GC and equator
    inc, nlonlat = gc_inc(alonlat, blonlat)
    # Convert to radians
    alonlat = np.radians(alonlat)
    blonlat = np.radians(blonlat)
    nlonlat = np.radians(nlonlat)
    # Distance of equatorial node to point A along GC
    da = np.arccos(np.cos(alonlat[0] - nlonlat) * np.cos(alonlat[1]))
    outda = np.degrees(da)
    db = np.arccos(np.cos(blonlat[0] - nlonlat) * np.cos(blonlat[1]))
    # Distance between the two points
    alatsign = alonlat[1] / np.abs(alonlat[1])
    blatsign = blonlat[1] / np.abs(blonlat[1])
    diffhemisign = alatsign * blatsign
    dd = np.abs(da - db * diffhemisign)
    # Find the mid point position between the two reference points
    dmid = -(da + db)/2.
    mlon = np.arctan(np.tan(dmid) * np.cos(inc)) + (nlonlat)
    mlat = np.arctan(np.tan(inc) * np.sin(mlon - nlonlat))
    mlonlat = np.array((mlon, mlat))
    outmid = np.degrees(mlonlat)
    if nonan is True:
        if np.isnan(dd) is True:
            dd = 0.
    # Convert the distance to degrees
    dist = np.degrees(dd)
    # When distances go above 180deg they are wrong! Subtract from 360.
    if (dist > 180.) is True:
        dist = 360. - dist
    #return various things
    return dist, outda, outmid

def gc_inc(alonlat, blonlat):
    """"
    Returns the inclination between the arc connecting two reference points along a GC and the equator from Leonard 1953
    Coordinates are in geographic longitude latitude
    alonlat: first reference point on great circle (GC), [longitude, latitude] in degrees
    blonlat: second reference point on great circle (GC), [longitude, latitude] in degrees
    """
    # Convert to radians
    alonlat = np.radians(alonlat)
    blonlat = np.radians(blonlat)
    # Try to correct for problem when determining angles when A and B are on either side of the equator
    lonshift = alonlat[0]
    # Get position of equatorial node nearest to reference mid point between A and B
    nlonlat = np.arctan( (np.sin(alonlat[0] - lonshift) * np.tan(blonlat[1]) - np.sin(blonlat[0] - lonshift) * np.tan(alonlat[1]))
                      / (np.cos(alonlat[0] - lonshift) * np.tan(blonlat[1]) - np.cos(blonlat[0] - lonshift) * np.tan(alonlat[1])) ) + lonshift
    #nlonlat = np.array((nlon, nlon.size-1)) #dont see why this is necessary so commenting out
    nlonlat = np.degrees(nlonlat)
    # Get inclination between GC and equator
    inc = np.arctan(np.tan(alonlat[1]) / np.sin(alonlat[0] - nlonlat))
    return inc, nlonlat

def ar_losgrad():

def pix_to_arc(map, x, y):
    """Convert pixel location to arcsecond location in map
    """
    arc_x = (x - (map.meta['crpix1'] - 1)) * map.meta['cdelt1'] + map.meta['crval1']
    arc_y = (y - (map.meta['crpix2'] - 1)) * map.meta['cdelt2'] + map.meta['crval2']
    return ((arc_x, arc_y))