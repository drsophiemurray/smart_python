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

cmpmm =  1e16 #Number centimeters in a mega meter
smoothphys = 16. #Physical gaussian smoothing HWHM in Mm. The radius of a characteristic supergranule
smooththresh = 15.0 # a segmentation threshold used in processing magnetograms (for 1kx1k gaussian smoothed; determined from plot_hmi_crosscal_resid.pro by comparing HMI and MDI detection masks)
strongthresh = 230.0 # a segmentation threshold to detect strong field fragments which are then dilated to form the 'cores' of ARs.
polethresh = .8 # imbalanced fraction of total area covered by a given (+ or -) polarity (=abs(npos-nneg)/abs(npos+nneg)), below which, the detection will be considered multipolar (core passes the flux balance test)


def ar_detect_core(map, smartmask):
    """
    """
    sz = map.data.shape
    maporig=mapmsk=map
    # the smoothing gaussian kernal HWHM
    smoothhwhm=smoothphys*cmpmm/ar_pxscale(map)
    # Smooth the data (used for finding the PSL and PSL mask)
    datasm = ar_grow(map.data, smoothhwhm, gauss=True)
    # Get ridge skeleton
    # use datasm, smooththresh
    ridgemask = ?ar_ridgemask(abs(datasm),thresh=params.smooththresh)
    # Get PSL map
    # use datasm, smoothwhm, smooththresh
    pslmask = ar_pslmask(datasm, smoothhwhm, smoothresh, skeleton=False)
    wpsl = np.where((ridgemask + pslmask) == 2)
    psltrace = np.zeros(sz)
    psltrace[wpsl] = 1.
    # Dilate PSL trace
    pslblobmask = ar_grow(psltrace, smoothhwhm, gauss=False)
    # Make strong field masks
    wstrong = np.where(np.abs(map.data) > strongthresh)
    strongmask = np.zeros(sz)
    strongmask[wstrong]=1.
    datastrongsm = ar_grow(np.abs(map.data*strongmask), smoothhwhm, gauss=True)
    # Determine strong blob pixels
    wstrongblob = np.where(np.abs(datastrongsm) > smooththresh)
    strongblobmask = np.zeros(sz)
    strongblobmask[wstrongblob] = 1.
    # Trim the PSL blob array to where it overlaps with the strong blob array
    pslblobmasktrim = strongblobmask * pslblobmask
    # Region grow with PSL map as a base and strong field map as a kernel
    arcoremask = watershed(sobel(pslblobmasktrim),strongblobmask, mask=pslblobmasktrim)
    # Filter to only take multi polar detections -> if >95% of pixels (above strong field threshold) are of one polarity then disregard the detection.
    arcoremaskid = measure.label(arcoremask, background=0)
    nid = np.max(arcoremaskid)
    polfracarr = np.zeros(nid)
    for j in range(1, nid+1):
        wthiscore = np.where(arcoremaskid == j)
        thiscoremask = np.zeros(sz)
        thiscoremask[wthiscore] = 1.
        wpos = np.where((strongblobmask*thiscoremask*map.data) > strongthresh)
        wneg = np.where((strongblobmask*thiscoremask*map.data) < -strongthresh)
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
    arcoremaskmpole = arcoremaskid < 1.
    # Re-index the final core mask
    arcoremaskmpoleid = measure.label(arcoremaskmpole, background=0)
    # Make combined core and normal AR detection using region grow
    # identifies regions that are attached to AR cores
    smartmask = smartmask < 1.
    arcoresmartcomb = watershed(sobel(smartmask), arcoremaskpole, mask=smartmask)
    # Change the values to uniquely identify the ARs and blobs so that only 1 output mask needs to be saved
    arcoremaskfinal = smartmask + (arcoresmartcomb * 2.) + (arcoremaskmpoleid * 100.)
    # Output result
    arcoremap = sunpy.map.Map(arcoremaskfinal, map.meta)
    pslmaskmap = sunpy.map.Map(strongblobmask*psltrace, map.meta)
    return arcoremap, pslmaskmap


def ar_ridgemask(data, thresh):
    """Pull out a ridge skeleton and insert into a mask to create a 1 pixel-wide trace for
    determining the PSL mask, using a watershed transform
    """
    sz=data.shape
    data[np.where(data<thresh)]=thresh


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
    negmaskgrow = ar_grow(negmask, radius, gauss=False)
    posmaskgrow = ar_grow(posmask, radius, gauss=False)
    # Determine the overlap of the two masks
    outmask = np.zeros(sz)
    outmask[np.where((posmaskgrow + negmaskgrow)  == 2)] = 1.
    if skeleton is True:
        return skeletonize(outmask)
    else:
        return outmask