
# params needed for config file
#Magnetogram Processing
cosmicthresh = 5000.  # A hard threshold for detecting cosmic rays

import numpy as np

def ar_processmag(map):
    """...
    """
    #already floats in numpy data array so skipped first line converting double to float
    imgsz = len(map.data) #not 1024 x 1024 like idl
    #didnt load parameters - need to add to config file
    #indextag didnt bother with
    map.data = cosmicthresh_remove(map.data)
    #Clean NaNs
    map.data = nan_remove(map.data)
    #Clean edge
    map.data = edge_remove(map.data)
    #Cosine map
    cosmap, limbmask = ar_cosmap(map)
    #add property to map object

def cosmicthresh_remove(data):
    """
    Search for cosmic rays using hard threshold defined in config file.
    Remove if greater than 3sigma detection than neighbooring pixels
    """
    wcosmic = np.where(data>cosmicthresh)
    ncosmic = len(wcosmic[0])
    print('Cosmic Ray Candidates Found: ', ncosmic)
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


def nan_remove(data):
    """
    Clean NaNs
    Includes zero-value as 'missing'
    """
    wnan = np.where(np.isnan(data))
    data[wnan] = -9999.
    fill_missing, data, -9999., 1. #IDL CODE - TO CHECK ON
    return data

def edge_remove(data):
    """
    Get rid of crazy arbitrary values at the edge of the image

    """
    edgepix = data[0, 0]
    wblankpx = np.where(data == edgepix)
    data[wblankpx] = 0.
    return data

def ar_cosmap(map):
    """
    Get the cosine map and off-limb pixel map using WCS
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

