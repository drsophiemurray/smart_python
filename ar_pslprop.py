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
        xrange = ((np.min(np.where(thismask == 1)[1])-1, np.max(np.where(thismask == 1)[1])+1))
        yrange = ((np.min(np.where(thismask == 1)[0])-1, np.max(np.where(thismask == 1)[0])+1))
        bottom_left_pixels = ((np.min(np.where(thismask == 1)[1])-1, np.min(np.where(thismask == 1)[0])-1))
        top_right_pixels = ((np.max(np.where(thismask == 1)[1])+1, np.max(np.where(thismask == 1)[0])+1))


def pix_to_arc(xrange, yrange):
    xmap =
    ymap =


def get_map_xrange(map):
    nx = map.data.shape[0]
    xc = map.meta["crval1"] + map.meta["cdelt1"] * (((map.meta["naxis1"] + 1) / 2) - map.meta["crpix1"])
    dx = map.meta['cdelt1']
    xmin = np.min(xc-dx*(nx-1.)/2.)
    xmax = np.max(xc+dx*(nx-1.)/2.)
    return ((xmin, xmax))

def get_map_yrange(map):
    ny = map.data.shape[0]
    yc = map.meta["crval2"] + map.meta["cdelt2"] * (((map.meta["naxis2"] + 1) / 2) - map.meta["crpix2"])
    dy = map.meta['cdelt2']
    ymin = np.min(yc-dy*(ny-1.)/2.)
    ymax = np.max(yc+dy*(ny-1.)/2.)
    return ((ymin, ymax))

def
    dxs = (max(xrange) - min(xrange)) / 2.
    dys = (max(yrange) - min(yrange)) / 2.
    xrange = [xrange[0] - dxs, xrange[0] + dxs]
    yrange = [yrange[1] - dys, yrange[1] + dys]


