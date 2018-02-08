'''
    SMART detection code
    ====================
    Written by Sophie A. Murray, code originally developed by Paul Higgins (ar_detect.pro).

    Developed under Python 3 and Sunpy 0.8.3
    - Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)

    Assumes input map has been (via ar_processmag):
    - rotated to solar north up
    - magnetic field values cosine corrected
    - offlimb pixels zeroed

    Inputs:
    - thismap: Processed magnetogram
    - limbmask: Mask with offlimb pixels
    - cmpmm: Number centimeters in a Mm.
    - smoothphys: Physical gaussian smoothing HWHM in Mm. The radius of a characteristic supergranule
    - smooththresh: Segmentation threshold used in processing magnetograms
    (for 1kx1k gaussian smoothed; determined from plot_hmi_crosscal_resid.pro by comparing HMI and MDI detection masks)
    - magthresh: Secondary segmentation threshold used to find all flux fragments in an image.
'''

from configparser import ConfigParser
import numpy as np
from sunpy.sun import constants
import sunpy.map
import scipy
from skimage import measure
from scipy import ndimage as nd
from skimage.morphology import watershed
from skimage import filters
import cv2


def ar_detect(thismap, limbmask):
    """
    Make detection mask of ARs with mask values ordered.
    """
    ## Load configuration file
    config = ConfigParser()
    config.read("config.ini")
    ## Initialise blank mask map
    sz = thismap.data.shape
    mask = np.zeros(sz)
    ## Gaussian smoothing
    # TO DO: same thing as ar_grow(thismap.data,smoothhwhm,gauss=True,kern=None)
    datasm, smoothhwhm = gauss_smooth(thismap,
                                      np.float(config.get('detection', 'smoothphys')),
                                      np.float(config.get('constants', 'cmpmm')))
    ## Make a mask of detections
    wmask = np.where(np.abs(datasm) > np.float(config.get('detection', 'smooththresh')))
    mask[wmask] = 1.
    ## Segment the non-smoothed magnetogram to grab near by fragments and connect adjacent blobs
    fragmask = np.zeros(sz)
    wfrag = np.where(np.abs(thismap.data) > np.float(config.get('detection', 'magthresh')))
    fragmask[wfrag] = 1.
    smfragmask = ar_grow(fragmask, smoothhwhm/2., gauss=False, kern=None)
    ## Region grow the smooth detections
    poismask = np.where(mask == 1.)
    grmask = watershed(filters.sobel(smfragmask), mask, mask=smfragmask)
#    grmask = cv2.dilate(grmask, np.ones((2,2), np.uint8), iterations=1) #make sure areas connect TO DO: TEST
    ## Mask offlimb pixels
    grmask = grmask*limbmask
    ## Return to 1s and 0s
    grmask[np.where(grmask < 0.5)] = 0. # really dont think this is needed in python
    grmask[np.where(grmask >= 0.5)] = 1.
    ## Separate the detections by assigning numbers
    maskfull, num = nd.label(grmask)
#   measure.label(grmask, background=0, return_num=True) #found other way above as this connected regions that werent connected in idl
    ## Order the detections by number of pixels
#    maskorder = ar_order_mask(maskfull, sz) # seems to drop a region, and dont think I need to order them for the next part as just goes binary again
#    maskmap = sunpy.map.Map(maskorder, thismap.meta)
    maskmap = sunpy.map.Map(maskfull, thismap.meta)
    return maskmap

def gauss_smooth(thismap, rsgrad, cmpmm):
    """
    Gaussian smooth the magnetogram
    """
    ## Get smoothing Gaussian kernel HWHM
    smoothhwhm = (rsgrad*cmpmm)/ar_pxscale(thismap, cmsqr=False, mmppx=False, cmppx=True)
    datasm = ar_grow(thismap.data, smoothhwhm, gauss=True, kern=None)
    return datasm, smoothhwhm

def ar_pxscale(thismap, cmsqr, mmppx, cmppx):
    """
    Calculate the area of an magnetogram pixel at disk-centre on the solar surface.
    """
    ## Area of pixel in Mm^2
    rsunmm = constants.get('radius').value/1e6
    mmperarcsec = rsunmm/thismap.meta["RSUN_OBS"] # Mm/arcsec
    pixarea = ((thismap.meta["CDELT1"] * mmperarcsec) * (thismap.meta["CDELT2"] * mmperarcsec)) # Mm^2
    if cmsqr is True:
        pixarea = pixarea*1e16
    ## Length of a side of a pixel
    retmmppx = (thismap.meta["CDELT1"]/thismap.meta["RSUN_OBS"])*rsunmm # Mm/px
    if cmppx is True:
        return retmmppx*1e16
    elif mmppx is True:
        return retmmppx
    else:
        return pixarea

def ar_grow(data, fwhm, gauss, kern):
    """
    Returns a dilated mask.
    Provide fwhm in pixels.
    If a NL mask is provided and gauss is TRUE, then the result will be a Shrijver R-mask.
    The convolution stucture will fall off as a gaussian.
    Notes:
    - For nice circular kernel binary kernel, width of kernel mask will be 2*fwhm+1,
    with a 1px boundary of 0 around the outside.
    2. Setting radius to 1 will result in a 3x3 structuring element for binary kernels,
    with a total array size of 5x5
    """
    gsig = fwhm / (np.sqrt(2. * np.log(2.)))
    imgsz=(int(4. * fwhm), int(4. * fwhm))
    struc = np.zeros(imgsz)
    ## Generate coordinate maps
    xcoord, ycoord, rcoord = xyrcoord(imgsz)
    struc[np.where(rcoord <= fwhm)] = 1.
    ## Crop to the edges of kernel with 1px boundary
    wxbound = ((np.min(np.where(struc.sum(axis=1))), np.max(np.where(struc.sum(axis=1)))))
    wybound = ((np.min(np.where(struc.sum(axis=0))), np.max(np.where(struc.sum(axis=0)))))
    struc = struc[(wxbound[0]-1):(wxbound[1] + 2),(wybound[0]-1):(wybound[1]+2)]
    if kern is None:
        struc = struc
    else:
        struc = kern
    struc[np.where(np.isnan(struc))] = 0.
    outkernal = struc
    ## Get Gaussian
    if gauss is True:
        mu = 0. # mean
#        max = 1.
        gstruc = gaussian(rcoord, mu, gsig)
        ## Normalize gstruc so that the volume is 1
        gstruc = gstruc/gstruc.sum()
        if kern is None:
            gstruc = gstruc
        else:
            gstruc = kern
        gstruc[np.where(np.isnan(gstruc))] = 0.
        outkernal = gstruc
        ## Convolve
        if (np.min(data.shape) > np.min(gstruc.shape)):
            return scipy.signal.convolve2d(data, outkernal, mode='same')
        else:
            print("ar_grow: kernel is too big compared to image!")
            return data
    else:
        if (np.min(data.shape) > np.min(struc.shape)):
            return scipy.ndimage.morphology.binary_dilation(data, structure=outkernal).astype(data.dtype)
        else:
            print("ar_grow: kernel is too big compared to image!")
            return data

def xyrcoord(imgsz):
    """
    Generate coordinate arrays.
    """
    rows = np.array(range(imgsz[0]))
    columns = np.array(range(imgsz[1]))
    ycoord,xcoord = np.meshgrid(columns, rows) #should i switch back? so confused!
    rcoord = np.sqrt((xcoord - imgsz[0] / 2.)**2. + (ycoord - imgsz[1] / 2.)**2)
    return xcoord, ycoord, rcoord

def gaussian(x, mu, sig):
    """
    Compute the 1-d Gaussian function.
    - x: independent variable of Gaussian function.
    - mu: mean value (center) of Gaussian.
    - sig: standard deviation (sigma) of Gaussian.
    Note: maximum value (factor) of Gaussian not included here like in IDL (= 1.)
    """
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def ar_order_mask(maskfull, sz):
    """
    Order the detections by number of pixels.
    NOTE: theres a bug in this with missing the last region. To fix someday if ever actually needed.
    """
    nar = np.max(maskfull)
    arnpix = np.histogram(maskfull, bins=range(nar+1))
    maskorder = np.zeros(sz)
    rank = np.argsort(-arnpix[0])
    for i in range(nar):
        maskorder[np.where(maskfull == rank[i])] = i
    return maskorder

if __name__ == '__main__':
    ar_detect()