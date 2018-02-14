'''
    SMART core detection code
    ====================
    Written by Sophie A. Murray, code originally developed by Paul Higgins (ar_detect_core.pro).

    Developed under Python 3 and Sunpy 0.8.3
    - Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)

    Provides detections of more limited extent after the initial ar_detect

    Inputs:
    - inmap: Processed magnetogram
    - inmask: Output mask from ar_detect
    - cmpmm: Number centimeters in a Mm.
    - smoothphys: Physical gaussian smoothing HWHM in Mm. The radius of a characteristic supergranule
    - smooththresh: Segmentation threshold used in processing magnetograms
    (for 1kx1k gaussian smoothed; determined from plot_hmi_crosscal_resid.pro by comparing HMI and MDI detection masks)
    - strongthresh: Segmentation threshold to detect strong field fragments which are then dilated to form the 'cores' of ARs.
    - polethresh: Imbalanced fraction of total area covered by a given (+ or -) polarity (=abs(npos-nneg)/abs(npos+nneg)),
    below which, the detection will be considered multipolar (core passes the flux balance test).
'''

from configparser import ConfigParser
from ar_detect import ar_pxscale, ar_grow
from skimage.morphology import watershed
import numpy as np
from skimage.morphology import skeletonize
from skimage import measure
from skimage import filters
import sunpy.map
from scipy import ndimage as nd


def ar_detect_core(inmap, inmask):
    """
    Make detection mask of ARs with mask values ordered
    to a more limited extent that ar_detect
    """
    ## Load configuration file
    config = ConfigParser()
    config.read("config.ini")
    ## Array size
    sz = inmap.data.shape
    #
    ## Get smoothing gaussian kernal HWHM
    smoothhwhm = np.float(config.get('detection', 'smoothphys'))*np.float(config.get('constants', 'cmpmm'))/ar_pxscale(inmap, cmsqr=False, mmppx=False, cmppx=True)
    ## Smooth the data (used for finding the PSL and PSL mask)
    datasm = ar_grow(inmap.data, smoothhwhm, gauss=True, kern=None)
    ## Get ridge skeleton
    ridgemask = ar_ridgemask(datasm, thresh=np.float(config.get('detection', 'smooththresh')))
    ## Get PSL map
    pslmask = ar_pslmask(datasm, smoothhwhm, np.float(config.get('detection', 'smooththresh')), skeleton=False)
    wpsl = np.where((ridgemask + pslmask) == 2)
    psltrace = np.zeros(sz)
    psltrace[wpsl] = 1.
    ## Dilate PSL trace
    pslblobmask = ar_grow(psltrace, smoothhwhm, gauss=False, kern=None)
    #
    ## Make strong field masks
    wstrong = np.where(np.abs(inmap.data) > np.float(config.get('detection', 'strongthresh')))
    strongmask = np.zeros(sz)
    strongmask[wstrong]=1.
    datastrongsm = ar_grow(np.abs(inmap.data*strongmask), smoothhwhm, gauss=True, kern=None)
    ## Determine strong blob pixels
    wstrongblob = np.where(np.abs(datastrongsm) > np.float(config.get('detection', 'smooththresh')))
    strongblobmask = np.zeros(sz)
    strongblobmask[wstrongblob] = 1.
    #
    ## Trim the PSL blob array to where it overlaps with the strong blob array
    pslblobmasktrim = strongblobmask * pslblobmask
    ## Region grow with PSL map as a base and strong field map as a kernel
    distance = nd.distance_transform_edt(pslblobmasktrim)
    arcoremask = watershed(-distance, pslblobmasktrim, mask=strongblobmask)
    ## Filter to only take multi polar detections
    # If >95% of pixels (above strong field threshold) are of one polarity then disregard the detection.
    arcoremaskid = measure.label(arcoremask, background=0)
    nid = np.max(arcoremaskid)
    polfracarr = np.zeros(nid)
    for j in range(1, nid+1):
        wthiscore = np.where(arcoremaskid == j)
        thiscoremask = np.zeros(sz)
        thiscoremask[wthiscore] = 1.
        wpos = np.where((strongblobmask*thiscoremask*inmap.data) > np.float(config.get('detection', 'strongthresh')))
        wneg = np.where((strongblobmask*thiscoremask*inmap.data) < -np.float(config.get('detection', 'strongthresh')))
        npos = len(wpos[0])
        nneg = len(wneg[0])
        if (npos == 0.) or (nneg == 0.):
            polfrac = 1.
        else:
            polfrac = np.abs(npos - nneg) / np.float(npos + nneg)
        polfracarr[j - 1] = polfrac
        if (polfrac > np.float(config.get('detection', 'polethresh'))):
            arcoremaskid[wthiscore] = 0.
    #
    ## Reset all values to 1.
    # Final detection mask is original detections + core detections.
    arcoremaskmpole = arcoremaskid >= 1.
    ## Re-index the final core mask
    arcoremaskmpoleid = measure.label(arcoremaskmpole, background=0)
    ## Make combined core and normal AR detection using region grow
    # Identifies regions that are attached to AR cores
    inmask = inmask >= 1.
    distance = nd.distance_transform_edt(arcoremaskmpole)
    arcoresmartcomb = watershed(-distance, arcoremaskmpole, mask=inmask)
    ## Create resulting map outputs
    arcoremap = sunpy.map.Map(arcoremaskmpoleid, inmap.meta)
#    pslmaskmap = sunpy.map.Map(arcoremaskmpole*psltrace, inmap.meta)
    pslmaskmap = sunpy.map.Map(arcoremaskmpole*ar_grow(psltrace, 1., gauss=False, kern=None), inmap.meta)
    ## Optional outputs, that for now have not been included..
    smartmaskmap = sunpy.map.Map(inmask*psltrace, inmap.meta) #a SMART mask of 0's and 1's
    smartconnmap = sunpy.map.Map(arcoresmartcomb*psltrace, inmap.meta) #a SMART mask with only the blobs that are touching cores
    return arcoremap, pslmaskmap

def ar_ridgemask(data, thresh):
    """Pull out a ridge skeleton and insert into a mask to create a 1 pixel-wide trace for
    determining the PSL mask, using a watershed transform.
    Alternative:
    - watershed(data, markers=8, connectivity=8)
    """
#    sz = data.shape
#    data[np.where(data<thresh)] = thresh
#    data = -data
#    from skimage.morphology import medial_axis
#    return medial_axis(np.invert(abs(datasm)>thresh))
    return skeletonize(np.invert(abs(data)>thresh))

def ar_pslmask(data, radius, thresh, skeleton):
    """Create a PSL mask, given some LOS B input data.
    - radius: dilation radius in pixels
    - thresh: magnetic threshold in G
    If skeleton TRUE, make PSL thinner by returning the skeleton (could also use thin or skeletonize_3d).
    """
    sz = data.shape
    blank = np.zeros(sz)
    ## Make a negative and a positive mask
    negmask = np.zeros(sz)
    negmask[np.where(data < -thresh)] = 1.
    posmask = np.zeros(sz)
    posmask[np.where(data > thresh)] = 1.
    ## Dilate the polarity masks
    negmaskgrow = ar_grow(negmask, radius, gauss=False, kern=None)
    posmaskgrow = ar_grow(posmask, radius, gauss=False, kern=None)
    ## Determine the overlap of the two masks
    outmask = np.zeros(sz)
    outmask[np.where((posmaskgrow + negmaskgrow)  == 2)] = 1.
    if skeleton is True:
        return skeletonize(outmask)
    else:
        return outmask

def ar_core2mask(data):
    """
    Separates out the three-layered core mask.
    NO LONGER NEEDED...
    """
    core_mask = np.ceil( ( (data - 3.) > 0 ) / 100. )
    smartmask = (data-core_mask*100.)<1.
    smartconn = (data-core_mask*100.)>1.<2.
    return coremask, np.ceil(smartmask), np.ceil(smartconn)

if __name__ == '__main__':
    ar_detect_core()