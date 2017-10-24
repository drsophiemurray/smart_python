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

status=-1

#Physics
cmpmm =  1e16 #Number centimeters in a mega meter
smoothphys = 16. #Physical gaussian smoothing HWHM in Mm. The radius of a characteristic supergranule


def ar_detect(map):
    """
    """
    szorig = map.data.shape
    # initialise blank mask map
    maskmap = sunpy.map.Map(np.zeros(szorig), map.meta)


def gauss_smooth(map, rsgrad, cmpmm):
    """Gaussian smooth the magnetogram
    """
    # Get smoothing Gaussian kernel HWHM
    ssmoothwhm = (rsgrad*cmpmm)/ar_pxscale(map)
    datasm = ar_grow(map.data, smoothhwhm)

def ar_pxscale(map):
    """Calculate the area of an magnetogram pixel at disk-centre on the solar surface.
    """
    #Area of pixel in Mm^2
    rsunmm = constants.get('radius').value/1e6
    mmperarcsec = rsunmm/map.meta["RSUN_OBS"] # Mm/arcsec
    pixarea = ((map.meta["CDELT1"] * mmperarcsec) * (map.meta["CDELT2"] * mmperarcsec)) # Mm^2
    #Length of a side of a pixel
    retmmppx = (map.meta["CDELT1"]/map.meta["RSUN_OBS"])*rsunmm # Mm/px
    return retmmppx*1e16

def ar_grow(arr):
    """
    Returns a dilated mask. If a NL mask is provided and GAUSSIAN is set, then the result will be a Shrijver R-mask.
    Provide RADIUS or FWHM in pixels. FWHM is actually half width at half max!!!
    The convolution stucture will fall off as a gaussian.
    Notes:
    1. For nice circular kernel binary kernel, width of kernel mask will be 2*radius+1, with a 1px boundary of 0 around the outside.
    2. Setting radius to 1 will result in a 3x3 structuring element for binary kernels, with a total array size of 5x5
    """
    arr0 = arr
    fwhm = radius0 = 5.
    gsig = fwhm / (np.sqrt(2. * np.log(2.)))
    imgsz=(int(4. * radius0), int(4. * radius0))
    struc = np.zeros(imgsz)
    # Generate coordinate maps
    xcoord, ycoord, rcoord = xyrcoord(imgsz)
    struc[np.where(rcoord <= radius0)] = 1.
    # Crop to the edges of kernel with 1px boundary
    wxbound = ((np.min(np.where(struc.sum(axis=1))), np.max(np.where(struc.sum(axis=1)))))
    wybound = ((np.min(np.where(struc.sum(axis=0))), np.max(np.where(struc.sum(axis=0)))))
    test=struc[(wxbound[0]-1):(wxbound[1] + 2),(wybound[0]-1):(wybound[1]+2)]

    # column total
     # row total

def xyrcoord(imgsz):
    """
    """
    rows = np.array(range(imgsz[0]))
    columns = np.array(range(imgsz[1]))
    ycoord,xcoord = np.meshgrid(columns,rows)
    rcoord = np.sqrt((xcoord - imgsz[0] / 2.)**2. + (ycoord - imgsz[1] / 2.)**2)
    return xcoord, ycoord, rcoord

