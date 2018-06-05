'''
    Solar Monitor Active Region Tracker
    ===================================
    Written by Sophie A. Murray, code originally developed by Paul Higgins (IDL SMART wrapper).
    - Code translated line-by-line and operates as close to IDL code as possible.

    Developed under Python 3 and Sunpy 0.8.3
    - Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)

    Steps:
    - Read HMI magnetogram and process it
    - Find 'magnetically interesting regions' and calculate properties of interest
    - Track regions from previous six hours
    - Visualise output for SolarMonitor.org

    TODO:
    - Merge Sean's tracking code
    - Finish de-projection code in ar_pslprop.py (merge from Sean's work)
    - Get rid of 'ar_' naming scheme and main function naming a.l.a Higgins version
    - Fix numbers in json files (all those unnecessary significant figures...)
    - Make better way to plot lat/lon grid rather than hacky version
'''


from configparser import ConfigParser
import json
import os
import time
import subprocess
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
import sunpy.map
from sunpy.visualization import wcsaxes_compat
import astropy.units as u
from astropy.coordinates import SkyCoord

import input_data
import process_magnetogram
import detect
import detect_core
import position_properties
import magnetic_properties
import psl_properties
from plot_detections import grid_overlay


def main(**fits_file):
    """

    :param file:
    :return:
    """
    # Load configuration file
    config = ConfigParser()
    config.read("config.ini")
    # Get directory
    data_dir = "".join([os.path.expanduser('~'), config.get('paths', 'data_dir')])
    # Load time to see how long this will take
    start_time = time.time()
    # First load the latest HMI data file
    if not fits_file:
        inmap = input_data.main(data_dir)
    else:
        inmap = sunpy.map.Map(data_dir+fits_file['fits_file'])

    # Downsample if 4096x4096 (generally shouldnt be if using near-real-time JSOC data)
    if inmap.dimensions[0].value != 1024:
        inmap = inmap.resample(u.Quantity([1024, 1024], u.pixel))
        inmap.meta['naxis1'] = 1024
        inmap.meta['naxis2'] = 1024

    # Now process magnetogram
    print('Processing data')
    processedmap, cosmap, limbmask = process_magnetogram.main(inmap, medianfilter=False)

    # Create AR masks
    print('Making core detections')
    detectionmap = detect.main(processedmap, limbmask)
    coredetectionmap, pslmap = detect_core.main(processedmap, detectionmap.data)
    print('SMART detections found: ', np.max(np.unique(coredetectionmap.data)))

    # Get properties
    print('Calculating properties')
    posprop = position_properties.main(processedmap, coredetectionmap.data, cosmap)
    magprop = magnetic_properties.main(processedmap, coredetectionmap.data, cosmap)
    pslprop = psl_properties.main(processedmap, coredetectionmap.data,
                                  doproj=False, projmaxscale=1024)

    # Output to json
    smartdate = inmap.date.strftime('%Y%m%d_%H%M')
    out = {'meta': {'dateobs': smartdate,
                    'dimension': inmap.dimensions[0].value,
                    'instrument': inmap.instrument}}
    out['posprop'] = posprop
    out['magprop'] = magprop
    out['pslprop'] = pslprop
    out = pd.io.json.dumps(out, data_dir+smartdate+'_properties.json')
    with open(data_dir+smartdate+'_properties.json', 'w') as outfile:
        outfile.write(out)
    # Beautifying it
    with open(data_dir+smartdate+'_properties.json') as infile:
        obj = json.load(infile)
    with open(data_dir+smartdate+'_properties.json', 'w') as outfile:
        json.dump(obj, outfile,
                  indent=4, separators=(',', ': '))

    # Visualise
    ## Just something simple for my testing
    ## - to be replaced by proper SolarMonitor stuff eventually...
    figure = plt.figure()
    # Get same axes
    bottom_left = SkyCoord(-1000*u.arcsec, -1000*u.arcsec, frame=processedmap.coordinate_frame)
    top_right = SkyCoord(1000*u.arcsec, 1000*u.arcsec, frame=processedmap.coordinate_frame)
    submap = processedmap.submap(bottom_left, top_right)
    axes = wcsaxes_compat.gca_wcs(processedmap.wcs)
    image = processedmap.plot(vmin=-500, vmax=500, axes=axes)
    axes.coords.grid(False)
    # Draw solar lat/lon grid
    overlay = grid_overlay(axes, grid_spacing=10 * u.deg)
    plt.colorbar(label='B [G]')
    # Overlay PILs and SMART detections
    plt.contour(pslmap.data, origin='lower',
                colors='lightblue', linewidths=1,
                vmin=0., vmax=np.max(np.unique(coredetectionmap.data))+1)
    plt.contour(coredetectionmap.data > 0., origin='lower',
                colors='blue', linewidths=1.0,
                vmin=0., vmax=np.max(np.unique(coredetectionmap.data))+1)
    plt.savefig(data_dir+smartdate+'_detections.eps')

    # Delete fits files for now for Met Office
    sys_call = "".join(['rm -r {}'.format(data_dir+'*.fits')])
    subprocess.call(sys_call, shell=True)

    # How long did that take?
    print('Runtime:', round(time.time() - start_time), 'seconds')


if __name__ == "__main__":
    main()
