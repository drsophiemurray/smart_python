# ;Input a processed data map and mask
# ;cosine correction for magnetic field values should already be done
# ;The tot. area, pos. area, neg. area, and
# ;	total, signed, fractional signed, negative, and positive flux
# ;	are determined
# ;
# ;STATUS = output status array keyword (one status for each AR)
# ;		0: initialised value
# ;		7: Everything went swimmingly, and there should be valid magnetic properties for each AR
# ;		1: No ARs were present in the mask! Should only happen if a blank array was read in

import numpy as np
from ar_detect import ar_pxscale
import pandas as pd

magthresh= 350.0 # F; a secondary segmentation threshold (for MDI) used to find all flux fragments in an image

def ar_magprop(map, mask, cosmap):
    """.	blankstr={datafile:indatafile[0],arid:0,areabnd:0., posareabnd:0., negareabnd:0., posarea:0., negarea:0., totarea:0., $
          bmax:0d, bmin:0d, bmean:0d, $
          totflx:0., imbflx:0., frcflx:0d, negflx:0., posflx:0.}
    """
    pxmmsq =  ar_pxscale(map, cmsqr=False, mmppx=False, cmppx=False)
    pxcmsq = ar_pxscale(map, cmsqr=True, mmppx=False, cmppx=False)
    nmask = np.max(mask)
    magdf = pd.DataFrame(columns = ['arid',
                                    'areabnd', 'posareanbnd', 'negareabnd',
                                    'posarea', 'negarea', 'totarea',
                                    'bmax', 'bmin', 'bmean',
                                    'posflx', 'negflx', 'totflx',
                                    'imbflx', 'frcflx'])

    #For each AR...
    for i in range(1, np.int(nmask)+1):
        # Zero pixels outside of detection boundary
        thismask = np.copy(mask)
        thismask[np.where(mask != i)] = 0.
        thisdat = np.copy(map.data)
        thisdat[np.where(mask != i)] = 0.
        thisabs = np.abs(thisdat)
        # Where are values within the detection boundary
        thismask[np.where(mask == i)] = 1.
        # Where values above mag threshold
        wthresh = np.where(thisabs >= magthresh)
        # Where negative values above thresh
        wneg = np.where((thisdat < 0) & (thisabs >= magthresh))
        # Where positive values above thresh
        wpos = np.where((thisdat > 0) & (thisabs >= magthresh))
        # Where negative values in boundary
        wnegbnd = np.where(thisdat < 0)
        # Where positive values in boundary
        wposbnd = np.where(thisdat > 0)
        # Magnetic moments calculated for values within boundary [G]
        bmax = np.max(thisdat[wval])
        bmin = np.min(thisdat[wval])
        bmean = np.mean(thisdat[wval])
        # Area of detection boundary [Mm^2]
        areabnd = np.sum(cosmap * thismask * pxmmsq)
        posareabnd = np.sum(cosmap[wposbnd] * thismask[wposbnd] * pxmmsq)
        negareabnd = np.sum(cosmap[wnegbnd] * thismask[wnegbnd] * pxmmsq)
        posarea = np.sum(cosmap[wpos] * thismask[wpos] * pxmmsq)
        negarea = np.sum(cosmap[wneg] * thismask[wneg] * pxmmsq)
        totarea = np.sum(cosmap[wthresh] * thismask[wthresh] * pxmmsq)
        # Flux measurements [Mx = G cm^2]
        totflx = np.sum(cosmap * thismask * pxcmsq * thisabs)
        imbflx = np.sum(cosmap * thismask * pxcmsq * thisdat)
        negflx = np.sum(cosmap[wneg] * thismask[wneg] * pxcmsq * thisabs[wneg])
        posflx = np.sum(cosmap[wpos] * thismask[wpos] * pxcmsq * thisabs[wpos])
        frcflx = (posflx - negflx) / totflx
        # Add to dataframe
        magdf = magdf.append([{'arid': i,
                             'areabnd': areabnd, 'posareanbnd': posareabnd, 'negareabnd': negareabnd,
                             'posarea': posarea, 'negarea': negarea, 'totarea': totarea,
                             'bmax': bmax, 'bmin': bmin, 'bmean': bmean,
                             'posflx': posflx, 'negflx': negflx, 'totflx': totflx,
                             'imbflx': imbflx, 'frcflx': frcflx}])
    return magdf


