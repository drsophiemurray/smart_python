# ;provide detections of more limited extent.
# ;
# ;STATUS =  -1 - the initialised value in the routine
# ;			0 - skiped file because no detections were found by normal SMART run
# ;			1 - skipped file because no PSL detections were found
# ;			2 - no ridge/PSL mask overlay pixels found
# ;			3 - no strong pixels found
# ;			4 - no strong blobby pixels found
# ;			5 - no stong-psl mask overlap found (using region grow)
# ;			6 - no final core detections (failed multi-polarity test)
# ;			7 - Core detections were found (lucky number seven!)
# ;
# ;NOTES:
# ;	1. 	For each AR, id values go like 100,200,300
# ;	2.	All positions with SMART mask pixels have 1. added to the mask value allowing one to subtract 1 and remove all SMART masks
# ;	3. 	If the SMART mask pixels are touching a core detection the 2. is added, giving a value of 2. or 3.
# ;	4. 	Thus to isolate only the core detections and have them increment by 1 (instead of 100), one would do: core_mask = ceil( ( ( combined_mask - 3. ) > 0 ) / 100. )
# ;

#params

from ar_detect import ar_pxscale, ar_grow
from skimage.morphology import watershed
import numpy as np
from skimage.morphology import skeletonize
from skimage import measure
from skimage import filters
import sunpy.map
from scipy import ndimage as nd

cmpmm =  1e16 #Number centimeters in a mega meter
smoothphys = 16. #Physical gaussian smoothing HWHM in Mm. The radius of a characteristic supergranule
smooththresh = 15.0 # a segmentation threshold used in processing magnetograms (for 1kx1k gaussian smoothed; determined from plot_hmi_crosscal_resid.pro by comparing HMI and MDI detection masks)
strongthresh = 230.0 # a segmentation threshold to detect strong field fragments which are then dilated to form the 'cores' of ARs.
polethresh = .8 # imbalanced fraction of total area covered by a given (+ or -) polarity (=abs(npos-nneg)/abs(npos+nneg)), below which, the detection will be considered multipolar (core passes the flux balance test)


def ar_detect_core(thismap, smartmask):
    """
    """
    sz = thismap.data.shape
   # maporig=mapmsk=thismap
    # the smoothing gaussian kernal HWHM
    smoothhwhm = smoothphys*cmpmm/ar_pxscale(thismap, cmsqr=False, mmppx=False, cmppx=True)
    # Smooth the data (used for finding the PSL and PSL mask)
    datasm = ar_grow(thismap.data, smoothhwhm, gauss=True, kern=None)
    # Get ridge skeleton
    # use datasm, smooththresh
    ridgemask = ar_ridgemask(datasm,thresh=smooththresh)
    # Get PSL map
    # use datasm, smoothwhm, smooththresh
    pslmask = ar_pslmask(datasm, smoothhwhm, smooththresh, skeleton=False)
    wpsl = np.where((test + pslmask) == 2)
    psltrace = np.zeros(sz)
    psltrace[wpsl] = 1.
    # Dilate PSL trace
    pslblobmask = ar_grow(psltrace, smoothhwhm, gauss=False, kern=None)
    # Make strong field masks
    wstrong = np.where(np.abs(thismap.data) > strongthresh)
    strongmask = np.zeros(sz)
    strongmask[wstrong]=1.
    datastrongsm = ar_grow(np.abs(thismap.data*strongmask), smoothhwhm, gauss=True, kern=None)
    # Determine strong blob pixels
    wstrongblob = np.where(np.abs(datastrongsm) > smooththresh)
    strongblobmask = np.zeros(sz)
    strongblobmask[wstrongblob] = 1.
    # Trim the PSL blob array to where it overlaps with the strong blob array
    pslblobmasktrim = strongblobmask * pslblobmask
    # Region grow with PSL map as a base and strong field map as a kernel
    distance = nd.distance_transform_edt(pslblobmasktrim)
    arcoremask = watershed(-distance, pslblobmasktrim, mask=strongblobmask)
    # Filter to only take multi polar detections -> if >95% of pixels (above strong field threshold) are of one polarity then disregard the detection.
    arcoremaskid = measure.label(arcoremask, background=0)
    nid = np.max(arcoremaskid)
    polfracarr = np.zeros(nid)
    for j in range(1, nid+1):
        wthiscore = np.where(arcoremaskid == j)
        thiscoremask = np.zeros(sz)
        thiscoremask[wthiscore] = 1.
        wpos = np.where((strongblobmask*thiscoremask*thismap.data) > strongthresh)
        wneg = np.where((strongblobmask*thiscoremask*thismap.data) < -strongthresh)
        npos = len(wpos[0])
        nneg = len(wneg[0])
        if (npos == 0.) or (nneg == 0.):
            polfrac = 1.
        else:
            polfrac = np.abs(npos - nneg) / np.float(npos + nneg)
        polfracarr[j - 1] = polfrac
        if (polfrac > polethresh):
            arcoremaskid[wthiscore] = 0.
    # Final detection mask is original detections + core detections.
    # Reset all values to 1.
    arcoremaskmpole = arcoremaskid >= 1.
    # Re-index the final core mask
    arcoremaskmpoleid = measure.label(arcoremaskmpole, background=0)
    # Make combined core and normal AR detection using region grow
    # identifies regions that are attached to AR cores
    smartmask = smartmask >= 1.
    distance = nd.distance_transform_edt(arcoremaskmpole)
    arcoresmartcomb = watershed(-distance, arcoremaskmpole, mask=smartmask)
    # Output result
    arcoremap = sunpy.map.Map(arcoremaskmpoleid, thismap.meta)
    pslmaskmap = sunpy.map.Map(strongblobmask*psltrace, thismap.meta)
    smartmaskmap = sunpy.map.Map(smartmask*psltrace, thismap.meta) #a SMART mask of 0's and 1's
    smartconnmap = sunpy.map.Map(arcoresmartcomb*psltrace, thismap.meta) #a SMART mask with only the blobs that are touching cores
    return arcoremap, pslmaskmap

def ar_ridgemask(data, thresh):
    """Pull out a ridge skeleton and insert into a mask to create a 1 pixel-wide trace for
    determining the PSL mask, using a watershed transform
    """
    #sz=data.shape
    #data[np.where(data<thresh)]=thresh
    #data=-data
    # OR from skimage.morphology import medial_axis
    # medial_axis(np.invert(abs(datasm)>thresh))
    return skeletonize(np.invert(abs(data)>thresh)) #watershed(data, markers=8, connectivity=8)

def ar_pslmask(data, radius, thresh, skeleton):
    """Create a PSL mask, given some LOS B input data.
    Note, radius = dilation radius in pixels
          thresh = magnetic threshold in G
    If skeleton set, make psl thinner by returning the skeleton
    (could also use thin or skeletonize_3d).
    """
    sz = data.shape
    blank = np.zeros(sz)
    # Make a negative and a positive mask
    negmask = np.zeros(sz)
    negmask[np.where(data < -thresh)] = 1.
    posmask = np.zeros(sz)
    posmask[np.where(data > thresh)] = 1.
    # Dilate the polarity masks
    negmaskgrow = ar_grow(negmask, radius, gauss=False, kern=None)
    posmaskgrow = ar_grow(posmask, radius, gauss=False, kern=None)
    # Determine the overlap of the two masks
    outmask = np.zeros(sz)
    outmask[np.where((posmaskgrow + negmaskgrow)  == 2)] = 1.
    if skeleton is True:
        return skeletonize(outmask)
    else:
        return outmask

def ar_core2mask(data):
    """NOT NEEDED"""
    core_mask = np.ceil( ( (data - 3.) > 0 ) / 100. )
    smartmask = (data-core_mask*100.)<1.
    smartconn = (data-core_mask*100.)>1.<2.
    return coremask, np.ceil(smartmask), np.ceil(smartconn)

if __name__ == '__main__':
    ar_detect_core()