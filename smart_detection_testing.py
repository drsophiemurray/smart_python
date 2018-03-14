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
    - Finish de-projection code in ar_pslprop.py
    - Get rid of 'ar_' naming scheme
    - Put stuff in sub codes under main and not their own function a.l.a Higgo
    - Fix numbers in json files (all those unnecessary significant figures...)
    - Make better way to plot lat/lon grid rather than hacky version
'''


from ar_readmag import ar_readmag
from ar_processmag import ar_processmag
from ar_detect import ar_detect
from ar_detect_core import ar_detect_core
from ar_posprop import ar_posprop
from ar_magprop import ar_magprop
from ar_pslprop import ar_pslprop
import astropy.units as u
import sunpy.map
import matplotlib.pylab as plt
import pandas as pd
import datetime
import numpy as np
import json
import time
from ar_plot import grid_overlay
from sunpy.visualization import wcsaxes_compat
from astropy.coordinates import SkyCoord

if __name__ == "__main__":
    # First load the latest HMI data file
    start_time = time.time()
    import os
    data_dir = "/Users/sophie/data/smart/python_test/"
    for file in os.listdir(data_dir):
        thismap = sunpy.map.Map(data_dir+file)

    # Need to downsample if 4096x4096
    # (generally shouldnt be if using near-real-time JSOC data)
        if thismap.meta['naxis1'] != 1024:
            thismap = thismap.resample(u.Quantity([1024, 1024], u.pixel))
            thismap.meta['naxis1'] = 1024
            thismap.meta['naxis2'] = 1024

    # Now process magnetogram
        print('Processing data')
        magproc, cosmap, limbmask = ar_processmag(thismap, medianfilter=False)

    # Create AR masks
        print('Making core detections')
        thissm = ar_detect(magproc, limbmask)
        thisar, pslmap = ar_detect_core(magproc, thissm.data)
        print('SMART detections found: ', np.max(np.unique(thisar.data)))

    # Get properties
        print('Calculating properties')
        posprop = ar_posprop(magproc, thisar.data, cosmap)
        magprop = ar_magprop(magproc, thisar.data, cosmap)
        pslprop = ar_pslprop(magproc, thisar.data, doproj=False, projmaxscale=1024)

    # Output to json
        smartdate = thismap.date.strftime('%Y%m%d_%H%M')
        out = {'meta': {'dateobs': smartdate,
                        'dimension': thismap.dimensions[0].value,
                        'instrument': thismap.instrument}}
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
#                  sort_keys=True,
                  indent=4, separators=(',', ': '))

    # Visualise
    ## Just something simple for my testing - to be replaced by propert SolarMonitor stuff eventually...
        figure = plt.figure()
        bottom_left = SkyCoord(-1000 * u.arcsec, -1000 * u.arcsec, frame=magproc.coordinate_frame)
        top_right = SkyCoord(1000 * u.arcsec, 1000 * u.arcsec, frame=magproc.coordinate_frame)
        submag = magproc.submap(bottom_left, top_right)
        subar = thisar.submap(bottom_left, top_right)
        subpsl = pslmap.submap(bottom_left, top_right)
        axes = wcsaxes_compat.gca_wcs(submag.wcs)
        image = submag.plot(vmin=-500, vmax=500, axes=axes)
        axes.coords.grid(False)
        #limb nans
        limbmask[np.where(limbmask == 0.)] = np.nan
        limbmap = sunpy.map.Map(limbmask, magproc.meta)
        sublimb = limbmap.submap(bottom_left, top_right)
        submag.data = submag.data*sublimb.data
        subar.data = subar.data * sublimb.data
        subpsl.data = subpsl.data * sublimb.data
        # Draw solar lat/lon grid
        overlay = grid_overlay(axes, grid_spacing=10 * u.deg)
#    plt.colorbar(label='B [G]')
    # Overlay PILs and SMART detections
 #       plt.contour(subpsl.data, origin='lower',
 #               colors='yellow', linewidths=0.5,
 #               vmin=0., vmax=np.max(np.unique(thisar.data))+1)
        plt.contour(subar.data>0., origin='lower',
                colors='blue', linewidths=1.0,
                vmin=0., vmax=np.max(np.unique(subar.data))+1)
        plt.savefig(data_dir+smartdate+'.eps')
        plt.close()

    # Output SMART detection and map
        thisar.save(data_dir+smartdate+'_detections.fits')
        magproc.save(data_dir+smartdate+'_map.fits')

    # How long did that take?
        print(str(smartdate))
        print('Runtime:', round(time.time() - start_time),'seconds')