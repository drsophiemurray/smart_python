'''
    SMART magnetic property code
    ============================
    Written by Sophie A. Murray, code originally developed by Paul Higgins (ar_magprop.pro).

    Developed under Python 3 and Sunpy 0.8.3
    - Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)

    Provides magnetic information of detected SMART regions, i.e.
    - Total, positive, negative, area
    - Total, postive, negative, signed, fractional signed, flux
    - Max, min, mean, magnetic field strength

    Inputs:
    - inmap: Processed magnetogram
    - inmask: Output SMART mask from ar_detect_core
    - cosmap: Cosine map from processing
    - magthresh: Secondary segmentation threshold (for MDI) used to find all flux fragments in an image
'''

from configparser import ConfigParser
import numpy as np
from ar_detect import ar_pxscale
import pandas as pd

def ar_magprop(inmap, inmask, cosmap):
    """
    Determine simple AR magnetic properties.
    """
    ## Load configuration file
    config = ConfigParser()
    config.read("config.ini")
    ## Set up parameters and output dataframe
    pxmmsq =  ar_pxscale(inmap, cmsqr=False, mmppx=False, cmppx=False)
    pxcmsq = ar_pxscale(inmap, cmsqr=True, mmppx=False, cmppx=False)
    nmask = np.max(inmask)
    magdf = pd.DataFrame(columns = ['arid',
                                    'areabnd', 'posareanbnd', 'negareabnd',
                                    'posarea', 'negarea', 'totarea',
                                    'bmax', 'bmin', 'bmean',
                                    'posflx', 'negflx', 'totflx',
                                    'imbflx', 'frcflx'])
    ## For each AR...
    for i in range(1, np.int(nmask)+1):
        ## Zero pixels outside of detection boundary
        thismask = np.copy(inmask)
        thismask[np.where(inmask != i)] = 0.
        thisdat = np.copy(inmap.data)
        thisdat[np.where(inmask != i)] = 0.
        thisabs = np.abs(thisdat)
        ## Where are values within the detection boundary
        thismask[np.where(inmask == i)] = 1.
        ## Where values above mag threshold
        wthresh = np.where(thisabs >= np.float(config.get('detection', 'magthresh')))
        ## Where negative values above thresh
        wneg = np.where((thisdat < 0) & (thisabs >= np.float(config.get('detection', 'magthresh'))))
        ## Where positive values above thresh
        wpos = np.where((thisdat > 0) & (thisabs >= np.float(config.get('detection', 'magthresh'))))
        ## Where negative values in boundary
        wnegbnd = np.where(thisdat < 0)
        ## Where positive values in boundary
        wposbnd = np.where(thisdat > 0)
        ## Magnetic moments calculated for values within boundary [G]
        bmax = np.max(thisdat[np.where(inmask == i)])
        bmin = np.min(thisdat[np.where(inmask == i)])
        bmean = np.mean(thisdat[np.where(inmask == i)])
        ## Area of detection boundary [Mm^2]
        areabnd = np.sum(cosmap * thismask * pxmmsq)
        posareabnd = np.sum(cosmap[wposbnd] * thismask[wposbnd] * pxmmsq)
        negareabnd = np.sum(cosmap[wnegbnd] * thismask[wnegbnd] * pxmmsq)
        posarea = np.sum(cosmap[wpos] * thismask[wpos] * pxmmsq)
        negarea = np.sum(cosmap[wneg] * thismask[wneg] * pxmmsq)
        totarea = np.sum(cosmap[wthresh] * thismask[wthresh] * pxmmsq)
        ## Flux measurements [Mx = G cm^2]
        totflx = np.sum(cosmap * thismask * pxcmsq * thisabs)
        imbflx = np.sum(cosmap * thismask * pxcmsq * thisdat)
        negflx = np.sum(cosmap[wneg] * thismask[wneg] * pxcmsq * thisabs[wneg])
        posflx = np.sum(cosmap[wpos] * thismask[wpos] * pxcmsq * thisabs[wpos])
        frcflx = (posflx - negflx) / totflx
        ## Add to dataframe
        magdf = magdf.append([{'arid': i,
                             'areabnd': areabnd, 'posareanbnd': posareabnd, 'negareabnd': negareabnd,
                             'posarea': posarea, 'negarea': negarea, 'totarea': totarea,
                             'bmax': bmax, 'bmin': bmin, 'bmean': bmean,
                             'posflx': posflx, 'negflx': negflx, 'totflx': totflx,
                             'imbflx': imbflx, 'frcflx': frcflx}], ignore_index=True)
    return magdf