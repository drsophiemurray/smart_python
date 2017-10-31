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

import numpy as np
import sunpy.map
from ar_detect import xyrcoord
from ar_posprop import px2hc, hc2hg

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


def ar_bipolesep(map):
    """Determine the flux-weighted bipole separation distance between the pos and neg centroids
    ;in degrees if a map is input, or in Px if only an image is input
    """
    image = map.data[::-1]
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
    pxsep = np.sum((pxpxloc-nxpxloc)^2.+(pypxloc-nypxloc)^2.)
    # Now the map outputs
    xc = map.meta["crval1"] + map.meta["cdelt1"]*(((map.meta["naxis1"] + 1)/2) - map.meta["crpix1"])
    yc = map.meta["crval2"] + map.meta["cdelt2"] * (((map.meta["naxis2"] + 1) / 2) - map.meta["crpix2"])
    phcxflx, phcyflx = px2hc(pxpxloc, pypxloc, map.meta['cdelt1'], map.meta['cdelt2'], xc, yc, imgsz)
    phgxflx, phgyflx, carpxflx = hc2hg(map, phcxflx, phcyflx)
    nhcxflx, nhcyflx = px2hc(nxpxloc, nypxloc, map.meta['cdelt1'], map.meta['cdelt2'], xc, yc, imgsz)
    nhgxflx, nhgyflx, carnxflx = hc2hg(map, nhcxflx, nhcyflx)
# sepstr.plon = phgxflx
# sepstr.plat = phgyflx
# sepstr.nlon = nhgxflx
# sepstr.nlat = nhgyflx
    gc_dist(((phgxflx, phgyflx)), ((nhgxflx, nhgyflx)))

def gc_dist():


def pix_to_arc(map, x, y):
    """Convert pixel location to arcsecond location in map
    """
    arc_x = (x - (map.meta['crpix1'] - 1)) * map.meta['cdelt1'] + map.meta['crval1']
    arc_y = (y - (map.meta['crpix2'] - 1)) * map.meta['cdelt2'] + map.meta['crval2']
    return ((arc_x, arc_y))


