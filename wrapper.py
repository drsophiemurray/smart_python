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
    - Finish de-projection code in ar_pslprop.py (merge from Sean's work)
    - Fix numbers in json files (all those unnecessary significant figures...)
    - Make better way to plot lat/lon grid rather than hacky version
    - Add meta data to the properties so you know what they are
'''


from configparser import ConfigParser
import json
import os
import time
import subprocess
import numpy as np
import pandas as pd
import sunpy.map
import astropy.units as u

import input_data
import process_magnetogram
import detect
import detect_core
import position_properties
import magnetic_properties
import psl_properties
import plot_detections


def main(*fits_file):
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
        inmap = sunpy.map.Map(data_dir+fits_file[0])

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
    plot_detections.main(processedmap, coredetectionmap,  pslmap,
                         data_dir, smartdate)

    # Delete fits files
    sys_call = "".join(['rm -r {}'.format(data_dir+'*.fits')])
    subprocess.call(sys_call, shell=True)

    # How long did that take?
    print('Runtime:', round(time.time() - start_time), 'seconds')


if __name__ == "__main__":
    main()
