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
mdi_noisethresh = 70.0 #; F; a noise threshold used in processing MDI magnetograms
psl_grad =  250. # F; a gradient threshold for determining the location of strong-gradient PSL
r_kernsz = 10 # I; the FWHM of the smoothing kernal for calculating Schrijver R value


import numpy as np
import sunpy.map
from ar_detect import xyrcoord, ar_grow, ar_pxscale
from ar_posprop import px2hc, hc2hg
from sunpy.sun import constants
import pandas as pd
import scipy.interpolate
import scipy.ndimage
import astropy.units as u
from ar_processmag import remove_nans
from skimage.morphology import skeletonize
from skimage import measure
from astropy.convolution import convolve, Box2DKernel


def ar_pslprop(map, mask, doproj, projmaxscale):
    """.arpslstr={arid:0,psllength:0.,pslsglength:0.,pslcurvature:0d,rvalue:0.,wlsg:0.,bipolesep_mm:0.,bipolesep_px:0.,bipolesep_proj:0.}
    """
    sz = map.data.shape
    xscale = map.meta['cdelt1']
    nmask = np.max(mask)
    psldf = pd.DataFrame(columns = ['arid',
                                    'psllength', 'pslsglength', 'pslcurvature',
                                    'rvalue', 'wlsg',
                                    'bipolesep_mm', 'bipolesep_px'])#TO DO ADD FOLLOWING:, 'bipolesep_proj'])
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
        #commented out below as dont get the point of it (seems to divide by itself later to equal zero)
        #if dobppxscl is True:
        #    projmmscl = ar_pxscale(submag, cmsqr=False, mmppx=True, cmppx=False) * projpxscl_bpsep
        #else:
        #    projmmscl = ar_pxscale(submag, cmsqr=False, mmppx=True, cmppx=False) * projpxscl
        projmmscl = ar_pxscale(submag, cmsqr=False, mmppx=True, cmppx=False) * projpxscl_bpsep
        kernpsl = [[0, 0, 1, 1, 1, 0, 0], [0, 1, 1, 1, 1, 1, 0], [1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1],
                   [1, 1, 1, 1, 1, 1, 1], [0, 1, 1, 1, 1, 1, 0], [0, 0, 1, 1, 1, 0, 0]]
        kernpsl = np.array(kernpsl)
        kernsz = kernpsl.shape
        # Resize the kernel based on the scale conversion
        if ((np.min(kernsz[0]/projpxscl_bpsep) < 1) is True) or (np.isnan(np.min(kernsz[0]/projpxscl_bpsep)) is True):
            kernpsl = rebin(kernpsl, (kernsz[0], kernsz[1]))
        else:
            factor = (kernsz[0]/projpxscl_bpsep)/kernsz[0]
            kernpsl = scipy.ndimage.zoom(kernpsl, factor) #need to get congrid working but this will do for now
        projmagg = ar_grow(projmag, 1, gauss = True, kern=None)
        psz = projmagg.shape
        nmask = np.zeros(projmagg.shape)
        pmask = np.copy(nmask)
        nmask[np.where(projmagg < (-mdi_noisethresh*2))] = 1.
        pmask[np.where(projmagg > mdi_noisethresh*2)] = 1.
        pmaskg = ar_grow(pmask, 1./2., gauss=False, kern=kernpsl)
        nmaskg = ar_grow(nmask, 1./2., gauss=False, kern=kernpsl)
        pslmask = np.zeros(projmagg.shape)
        pslmask[np.where(pmaskg + nmaskg == 2)] = 1.
        gradmag = ar_losgrad(projmagg)
        mapscl = ar_pxscale(map, cmsqr=False, mmppx=True, cmppx=False)
        gradpsl = pslmask*gradmag*projmmscl/mapscl
        pslmaskthresh = np.copy(pslmask)
        pslmaskthresh[np.where(gradpsl < psl_grad)] = 0.
        pslmaskt = skeletonize(pslmask)
        pslmaskt_thresh = skeletonize(pslmaskthresh)
        # Find the largest segment of PSL and indicate terminals
        pslmaslt_skel = skeletonize(ar_largest_blob(pslmask, gradpsl))
        # Large commented out section skipped
        # Determine the longest PSLs Skeleton length and curvature
        pslcurvature = 0.
        meanmmscl = np.mean(projmmscl)
        psllength = np.sum(pslmaskt * projmmscl)
        psllengtht = np.sum(pslmaskt_thresh * projmmscl) #strong
        # Determine  R
        # compute pos and neg polarity maps, with product defining polarity inversion line:
        prim = np.copy(rim)
        prim[np.where(rim < 150.)] = 0.
        p1p = convolve(prim, Box2DKernel(3))
        p1p[np.where(p1p > 0)] = 1.
        nrim = np.copy(rim)
        nrim[np.where(rim > -150.)] = 0.
        p1n = convolve(nrim, Box2DKernel(3))
        p1n[np.where(p1n < 0)] = 1.
        pmap = ar_r_smear((p1p * p1n), r_kernsz)
        rmap = pmap*np.abs(rim)
        rmap[np.where(rmap < 0.)] = 1.
        rmasked = rmap * projmask
        rmasked[np.where(np.isnan(rmasked))] = 0.
        thisr = np.sum(rmasked)
        # DETERMINE SUMMED GRADIENT (WLsg)
        wlsgdat = gradpsl * pslmask #thresh
        wlsgdat[np.where(np.isnan(wlsgdat))] = 0.
        thiswlsg = np.sum(wlsgdat)
        # Fill structure
        psldf = psldf.append([{'arid': i,
                               'psllength': psllength, 'pslsglength': psllengtht, 'pslcurvature': pslcurvature,
                               'rvalue': thisr, 'wlsg': thiswlsg,
                               'bipolesep_mm': bipsepstr['gcdist_mm'],
                               'bipolesep_px': bipsepstr['gcdist_px']}])
                               #'bipolesep_proj': bipsepstr['gcdist_proj']}])
        # TO DO - SHUOLD BE bisepstrproj for the last valeu:
        # thispslstr.bipolesep_proj = bipsepstrproj.pxsep
    return psldf

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
    gcdist_px = (gcdist_deg / 360.) * 2. * np.pi * (map.meta['rsun_obs'] / map.meta['cdelt1'])
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

def ar_losgrad(data):
    """
    Take the gradient in the horizontal plane of the LOS magnetic field.
    """
    # Buffer the image to reduce edge effects
    imgsz = data.shape
    dataint = np.zeros([imgsz[0]+10,imgsz[1]+10])
    dataint[:] = np.nan
    dataint[5:imgsz[0] + 5, 5:imgsz[1] + 5] = data
    dataint = remove_nans(dataint)
    xgrad = np.gradient(dataint)[1]
    ygrad = np.gradient(dataint)[0] # np.rot90(np.gradient(np.rot90(dataint,3))[1])
    gradmag = np.sqrt(xgrad**2. + ygrad**2.)
    return gradmag[5:imgsz[0] + 5, 5:imgsz[1] + 5]

def ar_largest_blob(mask, data):
    """
    ;Input MASK with 1=feature, 0=quiet
    ;Returns MASK with all features zeroed except for the largest
    ;set flux to take the flux weighted largest blob
    """
    masksep = measure.label(mask, background=0)
    ncont = np.max(masksep)
    narr = np.zeros(ncont)
    for i in range(1, np.int(ncont)+1):
        narr[i - 1] = np.where(masksep == i)[0].size
    wnbest = np.int((np.where(narr == np.max(narr)))[0] + 1.)
    wbig = np.where(masksep == wnbest)
    w0 = np.where(masksep != wnbest)
    mask[w0] = 0.
    return mask

def ar_r_smear(image, szkernel):
    """
;+
;  NAME:
;    ar_r_smear
;  PURPOSE:
;    convolve an image with [default] a gaussian profile of FWHM width
;    n and boxwidth 4n, or alternatively with a specified kernel
;
;  CALLING SEQUENCE:
;    image=ar_r_smear(image,n,kernel=kernel)
;
;  INPUTS:
;    image  image to be processed
;    n      fwhm value of the gaussian smearing that is applied
;
;  OUTPUTS:
;           returns convolved image
;
;  OPTIONAL INPUT KEYWORDS:
;    kernal        convolution kernel, 2d real array
;   [edge_truncate as in IDL routine convol - obsolete 2006/04/17, see below]
;
;  OPTIONAL OUTPUT KEYWORDS:
;
;  METHOD:
;    uses standard IDL routine convol
;
;  MODIFICATION HISTORY:
;	C.J. Shcrijver - 11-Feb-2014 - written
;	P.A. Higgins - 12-Feb-2014 - modified using M. Bobra's suggestion (changed
;		kernal width to 4*n+1 rather than 4n) and standardised code to fit
;		within the SMART_LIBRARY repository:
;		http://github.com/pohuigin/smart_library/
;
;-
    """
    n = szkernel
    sigma = n / (2. * np.sqrt(2. * np.log(2.)))
    kernel = np.zeros(np.int(4 * n + 1.))
    for i in range(0, len(kernel)):
        kernel[i] = np.exp(-(i-(2 * n - 0.5)) ** 2 / (2 * sigma ** 2))
    kernel = np.outer(kernel, kernel)
    kernel = kernel/np.sum(kernel)
    return convolve(image, kernel)


def pix_to_arc(map, x, y):
    """Convert pixel location to arcsecond location in map
    """
    arc_x = (x - (map.meta['crpix1'] - 1)) * map.meta['cdelt1'] + map.meta['crval1']
    arc_y = (y - (map.meta['crpix2'] - 1)) * map.meta['cdelt2'] + map.meta['crval2']
    return ((arc_x, arc_y))