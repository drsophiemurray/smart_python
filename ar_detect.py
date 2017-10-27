#Routine to make a detection mask (map structure) of ARs with mask values ordered from largest to smallest
#Assumes input mask has been:
#	-rotated to solar north up
#	-magnetic field values cosine corrected
#	-offlimb pixels zeroed
#DOPROCESS = set to do the above processing
#
#	If running MDI data, it is suggested to read in raw magnetograms and
#	set /DOPROCESS because the smoothed and unsmoothed masks will use
# 	different processing due to the MDI noise problem.
#
#MAPPROC = Pull out the processed map
#REBIN4k21k = Do the detections on a magnetogram rebinned to 1kx1k
#STATUS = output keyword indicating whether detections were found or
#not
#		   -1 - The initialised value
#			0 - Detections were found
#			1 - No detections found in gaussian mask (aborted)
#			2 - No detections found and fragment mask had no detections (souldn't occur- in this case status should be 0 if there were detections and 1 if no detections)
#			3 - Region-grown mask had no detections (aborted; shouldn't occur... might mean error in code)
#			4 - Final indexed mask had no detections (aborted; shouldn't occur... might mean error in code)

import numpy as np
from sunpy.sun import constants
import sunpy.map
import scipy
from skimage import measure
from scipy import ndimage as nd
from skimage.morphology import watershed

status=-1

#params
cmpmm =  1e16 #Number centimeters in a mega meter
smoothphys = 16. #Physical gaussian smoothing HWHM in Mm. The radius of a characteristic supergranule
smooththresh = 15.0 # a segmentation threshold used in processing magnetograms (for 1kx1k gaussian smoothed; determined from plot_hmi_crosscal_resid.pro by comparing HMI and MDI detection masks)
magthresh =  350.0 #a secondary segmentation threshold (for MDI) used to find all flux fragments in an image


def ar_detect(map, limbmask):
    """
    """
    sz = map.data.shape
    ## Initialise blank mask map
    mask = np.zeros(sz)
    ## Gaussian smoothing
    datasm, smoothhwhm = gauss_smooth(map, smoothphys, cmpmm)
    ## Make a mask of detections
    wmask = np.where(np.abs(datasm) > smooththresh)
    mask[wmask] = 1.
    ## Segment the non-smoothed magnetogram to grab near by fragments and connect adjacent blobs
    fragmask = np.zeros(sz)
    wfrag = np.where(np.abs(map.data) >  magthresh)
    fragmask[wfrag] = 1.
    smfragmask = ar_grow(fragmask, smoothhwhm/2., gauss=False)
    ## Region grow the smooth detections
    poismask = np.where(mask == 1.)
    # wgrow=region_grow(smfragmask,poismask,thresh=[0.5,1.5]) #array, where to grow, threshold between which new region should fall
    #grmask = mask
    #grmask[wgrow] = 1
    # https://stackoverflow.com/questions/31848309/classifying-python-array-by-nearest-seed-region
    #distance = nd.distance_transform_edt(smfragmask)
    #grmask = watershed(distance, mask, mask=smfragmask)
    grmask = watershed(sobel(smfragmask), mask, mask=smfragmask)
    ## Mask offlimb pixels
    grmask = grmask*limbmask
    ## Return to 1s and 0s ndi.binary_fill_holes(segmentation - 1)
    grmask[np.where(grmask < 0.5)] = 0.
    grmask[np.where(grmask >= 0.5)] = 1.
    ## Separate the detections by assigning numbers
    maskfull = measure.label(grmask, background=0)
    ## Order the detections by number of pixels
    nar = np.max(maskfull)
    arnpix = np.histogram(maskfull, bins=range(nar+1))
    maskorder = np.zeros(sz)
    rank = np.argsort(-arnpix[0])
    for i in range(nar):
        maskorder[np.where(maskfull == rank[i])] = i
    maskmap = sunpy.map.Map(maskorder, map.meta)
    return maskmap


def gauss_smooth(map, rsgrad, cmpmm):
    """Gaussian smooth the magnetogram
    """
    ## Get smoothing Gaussian kernel HWHM
    smoothhwhm = (rsgrad*cmpmm)/ar_pxscale(map)
    datasm = ar_grow(map.data, smoothhwhm, gauss=True) #to do: make gaussian an option
    return datasm, smoothhwhm

def ar_pxscale(map):
    """Calculate the area of an magnetogram pixel at disk-centre on the solar surface.
    """
    ## Area of pixel in Mm^2
    rsunmm = constants.get('radius').value/1e6
    mmperarcsec = rsunmm/map.meta["RSUN_OBS"] # Mm/arcsec
    pixarea = ((map.meta["CDELT1"] * mmperarcsec) * (map.meta["CDELT2"] * mmperarcsec)) # Mm^2
    ## Length of a side of a pixel
    retmmppx = (map.meta["CDELT1"]/map.meta["RSUN_OBS"])*rsunmm # Mm/px
    return retmmppx*1e16

def ar_grow(data, fwhm, gauss):
    """
    Returns a dilated mask. If a NL mask is provided and GAUSSIAN is set, then the result will be a Shrijver R-mask.
    Provide RADIUS or FWHM in pixels. FWHM is actually half width at half max!!!
    The convolution stucture will fall off as a gaussian.
    Notes:
    1. For nice circular kernel binary kernel, width of kernel mask will be 2*radius+1, with a 1px boundary of 0 around the outside.
    2. Setting radius to 1 will result in a 3x3 structuring element for binary kernels, with a total array size of 5x5
    # TO DO: MAKE GAUSSIAN BELOW AN OPTION
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
    struc[np.where(np.isnan(struc))] = 0.
    outkernal = struc
    ## Get Gaussian
    if gauss is True:
        mu = 0. #mean
        # max = 1.
        gstruc = gaussian(rcoord, mu, gsig)
        ## Normalize gstruc so that the volume is 1
        gstruc = gstruc/gstruc.sum()
        gstruc[np.where(np.isnan(gstruc))] = 0.
        outkernal = gstruc
        ## Convolve
        if (np.min(data.shape) > np.min(gstruc.shape)):
            return convolution2d(data, outkernal)
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
    """
    rows = np.array(range(imgsz[0]))
    columns = np.array(range(imgsz[1]))
    ycoord,xcoord = np.meshgrid(columns,rows)
    rcoord = np.sqrt((xcoord - imgsz[0] / 2.)**2. + (ycoord - imgsz[1] / 2.)**2)
    return xcoord, ycoord, rcoord

def gaussian(x, mu, sig):
    """
    """
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def convolution2d(image, kernel):
    """TO DO: scipy.signal.convolve2d
    """
    m, n = kernel.shape
    if (m == n):
        y, x = image.shape
        y = y - m + 1
        x = x - m + 1
        new_image = np.zeros((y,x))
        for i in range(y):
            for j in range(x):
                new_image[i][j] = np.sum(image[i:i+m, j:j+m]*kernel)
    return new_image

