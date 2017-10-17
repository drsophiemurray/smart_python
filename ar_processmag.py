
# params needed for config file
#Magnetogram Processing
cosmicthresh = 5000.  # A hard threshold for detecting cosmic rays

import numpy as np
from scipy import interpolate

def ar_processmag(map):
    """...
    """
    # already floats in numpy data array so skipped first line converting double to float
    imgsz = len(map.data) #not 1024 x 1024 like idl
    # didnt load parameters - need to add to config file
    # indextag didnt bother with
    data = cosmicthresh_remove(map.data)
    # Clean NaNs - Higgo used bilinear interpoaltion
    data = nan_remove(data)
    # Clean edge - make all pixels off limb equal to 0. -- commented out as done during nan removal above!
    # map.data = edge_remove(map.data)
    # Cosine map
    cosmap, limbmask = ar_cosmap(map)
    # Add property to map object
    new_map = sunpy.map.Map(data, map.meta)

def cosmicthresh_remove(data):
    """
    Search for cosmic rays using hard threshold defined in config file.
    Remove if greater than 3sigma detection than neighbooring pixels
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

def nan_remove(array):
    """
    Clean NaNs
    Includes zero-value as 'missing'
    """
    x = np.arange(0, array.shape[1])
    y = np.arange(0, array.shape[0])
    # mask invalid values
    array = np.ma.masked_invalid(array)
    xx, yy = np.meshgrid(x, y)
    # get only the valid values
    x1 = xx[~array.mask]
    y1 = yy[~array.mask]
    newarr = array[~array.mask]
    GD1 = interpolate.griddata((x1, y1), newarr.ravel(),
                               (xx, yy),
                               method='cubic', fill_value=0.)
    return GD1

def edge_remove(data):
    """
    Get rid of crazy arbitrary values at the edge of the image
    Basically set everything beyond the limb equal to zero
    """
    edgepix = data[0, 0]
    wblankpx = np.where(data == edgepix)
    data[wblankpx] = 0.
    return data

def ar_cosmap(map):
    """
    Get the cosine map and off-limb pixel map using WCS.
    ;Generate a map of the solar disk that is 1 at disk center and goes radially outward as the cos(angle to LOS) 
    ;(= 2 at 60 degrees from LOS)
    ;optionally output:	rrdeg = gives degrees from disk center
    ;					wcs = wcs structure from input map file
    ;					offlimb = map of 1=on-disk and 0=off-disk
    ;					edgefudge = take off an extra half percent from the disk to get rid of limb effects
    """
    return cosmap, limbmask

def offlimb(data, limbmask):
    """
    zero off-limb pixels
    zero from 80 degrees to LOS
    """

def median_filter():
    """
    Median filter noisy values
    Do a 3x3 median filter
    """

def cosine_correction():
    """
    Do magnetic field cosine correction
    Limit correction to having 1 pixel at edge of the Sun
    """

def solar_rot():
    """
    ;Rotate solar north = up
    """

