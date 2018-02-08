'''
    SMART HMI magnetgram .fits processing code
    =========================================
    Written by Sophie A. Murray, code originally developed by Paul Higgins (ar_processmag.pro).

    Developed under Python 3 and Sunpy 0.8.3
    - Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)

    Inputs:
    - inmap: processed magnetogram
    - cosmicthresh: a hard threshold for detecting cosmic rays
    - medfiltwidth: the width of a box to use to perform median filter on magnetogram

    Optional keywords:
    - medianfilter: if TRUE, 3x3 median filter noisy values (default FALSE in original code)

    Steps:
    - Cosmic ray removal
    - Non-finite removal
    - Zero off-limb pixels
    - Rotate solar north = up
    - Cosine correction
    - Median filter

    Notes:
    - sunpy.wcs has been deprecated and needs to be replaced
    - there has to be a better way for median filtering a.l.a filer_image.pro
'''

import numpy as np
import scipy
from scipy import interpolate
import sunpy.wcs.wcs as wcs
from configparser import ConfigParser
import sunpy.map
import astropy.units as u

def ar_processmag(inmap, medianfilter):
    """Load input paramters, remove cosmic rays and NaNs,
    then make all off-limb pixels zero, and clean up limb,
    rotate map and do a cosine correction.
    """
    ## Load configuration file
    config = ConfigParser()
    config.read("config.ini")
    ## Rotate
#    inmap = inmap.rotate(angle=int(inmap.meta['crota2'])*u.deg)
#    inmap = inmap.resample(u.Quantity([1024, 1024], u.pixel))
#    inmap.meta['crota2'] = 0.
    # I commented out above as it was just adding way too much limb noise, so instead doing a hacky way
    if (inmap.meta['crota2'] >= 100.):
        data = np.flip(inmap.data, 1)[::-1]
        inmap = sunpy.map.Map(data, inmap.meta)
        inmap.meta['crota2'] = 0.
    # Already floats in numpy data array so skipped first line converting double to float
    imgsz = len(inmap.data) #not 1024 x 1024 like idl
    # Didnt load parameters - need to add to config file
    # Didnt bother with indextag
    ## Cosmic ray removal
    data = cosmicthresh_remove(inmap.data, np.float(config.get('processing', 'cosmicthresh')))
    ## Clean NaNs
    # Higgo used bilinear interpoaltion
    data = remove_nans(data)
    ## Zero off-limb pixels
    # Clean edge - make all pixels off limb equal to 0. TO DO - can be commented out as done during nan removal above!
    data = edge_remove(data)
    ## Create cosine map
    inmap = sunpy.map.Map(data, inmap.meta)
    cosmap, rrdeg, limbmask = ar_cosmap(inmap)
    ## Fix remaining limb issues
    data, limbmask = fix_limb(inmap.data, rrdeg, limbmask)
    ## Median filter noisy values
    # TO DO: do it properly like Higgo - check sunpy stuff - also commented as not used in main program
    if medianfilter is True:
        data = median_filter(data, np.float(config.get('processing', 'medfiltwidth')))
    inmap = sunpy.map.Map(data, inmap.meta)
    ## Magnetic field cosine correction
    inmap = sunpy.map.Map(data, inmap.meta)
    data, cosmap = cosine_correction(inmap, cosmap)
    return inmap, cosmap, limbmask

def cosmicthresh_remove(data, cosmicthresh):
    """
    Search for cosmic rays using hard threshold defined in config file.
    Remove if greater than 3-sigma detection than neighbouring pixels.
    """
    wcosmic = np.where(data>cosmicthresh)
    ncosmic = len(wcosmic[0])
    print('Cosmic Ray Candidates Found: ', ncosmic)
    if (ncosmic == 0.):
        return data
    else:
        wcx = wcosmic[0]
        wcy = wcosmic[1]
        neighbours = np.int_([-1, 0, 1, 1, 1, 0, -1, -1])
        for i in range(0, ncosmic):
            wcx_neighbours = wcx[i] + neighbours
            wcy_neighbours = wcy[i] + neighbours
            wc_logic = (3*np.std(data[wcx_neighbours, wcy_neighbours]))+np.mean(data[wcx_neighbours, wcy_neighbours])
            if (data[wcx[i], wcy[i]] >= wc_logic):
                data[wcx[i], wcy[i]] = np.mean(data[wcx_neighbours, wcy_neighbours])
        return data

def remove_nans(array):
    """
    Clean NaNs.
    Includes zero-value as 'missing'.
    """
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    ## Mask invalid values
    array = np.ma.masked_invalid(array)
    xx, yy = np.meshgrid(x, y)
    ## Get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
                               (xx, yy),
                               method='cubic', fill_value=0.)
    return GD1

def edge_remove(data):
    """
    Get rid of crazy arbitrary values at the edge of the image.
    Basically set everything beyond the limb equal to zero.
    """
    edgepix = data[0, 0]
    wblankpx = np.where(data == edgepix)
    data[wblankpx] = 0.
    return data

def ar_cosmap(inmap):
    """
    Get the cosine map and off-limb pixel map using WCS.
    Generate a map of the solar disk that is 1 at disk center and goes radially outward as the cos(angle to LOS), which
    is = 2 at 60 degrees from LOS.
    Other outputs:
    - rrdeg: gives degrees from disk center
    - offlimb: map of 1=on-disk and 0=off-disk
    """
    ## Take off an extra half percent from the disk to get rid of limb effects
    fudge=0.999
    #
    ## Get helioprojective_coordinates
    xx, yy = wcs.convert_pixel_to_data(inmap.data.shape,
                                       [inmap.meta["CDELT1"], inmap.meta["CDELT2"]],
                                       [inmap.meta["CRPIX1"], inmap.meta["CRPIX2"]],
                                       [inmap.meta["CRVAL1"], inmap.meta["CRVAL2"]])
    rr = ((xx**2.) + (yy**2.))**(0.5)
    #
    coscor = np.copy(rr)
    rrdeg = np.arcsin(coscor / inmap.meta["RSUN_OBS"])
    coscor = 1. / np.cos(rrdeg)
    wgt = np.where(rr > (inmap.meta["RSUN_OBS"]*fudge))
    coscor[wgt] = 1.
    #
    offlimb = np.copy(rr)
    wgtrr = np.where(rr >= (inmap.meta["RSUN_OBS"]*fudge))
    offlimb[wgtrr] = 0.
    wltrr = np.where(rr < (inmap.meta["RSUN_OBS"]*fudge))
    offlimb[wltrr] = 1.
    #
    return coscor, rrdeg, offlimb

def fix_limb(data, rrdeg, limbmask):
    """
    Zero off-limb pixels (zero from 80 degrees to LOS).
    This is making the edge a bit smaller.
    """
    maxlimb = 80.
    wofflimb = np.where(((rrdeg/(2.*np.pi))*360.) > maxlimb)
    data[wofflimb] = 0.
    limbmask[wofflimb] = 0.
    return data*limbmask, limbmask

def median_filter(data, medfiltwidth):
    """
    Median filter noisy values. See here for inspiration:
    http://docs.sunpy.org/en/stable/generated/gallery/image_bright_regions_gallery_example.html
    """
    return scipy.ndimage.gaussian_filter(data, medfiltwidth)

def cosine_correction(inmap, cosmap):
    """
    Do magnetic field cosine correction.
    Limit correction to having 1 pixel at edge of the Sun. This is the maximum factor of pixel area
    covered by a single pixel at the solar limb as compared with at disk centre.
    """
    thetalim = np.arcsin(1. - inmap.meta["CDELT1"] / inmap.meta["RSUN_OBS"])
    coscorlim = 1. / np.cos(thetalim)
    cosmaplim = np.where((cosmap) > coscorlim)
    cosmap[cosmaplim] = coscorlim
    return inmap.data*cosmap, cosmap

def myround(x, base=5):
    """
    Round a number to nearest '5'.
    """
    return int(base * round(float(x)/base))

if __name__ == '__main__':
    ar_processmag()
