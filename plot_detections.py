'''
    First created 2018-02-01
    Sophie A. Murray

    Python Version:    Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)

    Description:
    Create a basic plot of SMART detections over processed magnetogram.
    Detections are numbered plus PILs also displayed.
'''

import matplotlib.pylab as plt
from sunpy.visualization import wcsaxes_compat
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

def main(processedmap, coredetectionmap, pslmap, data_dir, smartdate):
    """
    Plot and save a simple SMART detection image

    Parameters
    ----------
    processedmap : `~sunpy.map.sources.sdo.HMIMap`
        SMART-processed magnetogram onto which detections will be overlayed
    coredetectionmap : `~sunpy.map.sources.sdo.HMIMap`
        map of SMART detections to plot
    pslmap : `~sunpy.map.sources.sdo.HMIMap`
        SMART-detected PIL map
    data_dir : `~str`
        output directory location of where to save image
    smartdate : `~str`
        date of observation for image name
    """
    figure = plt.figure()
    # Get same axes
    bottom_left = SkyCoord(-1000*u.arcsec, -1000*u.arcsec,
                           frame=processedmap.coordinate_frame)
    top_right = SkyCoord(1000*u.arcsec, 1000*u.arcsec,
                         frame=processedmap.coordinate_frame)
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


def grid_overlay(axes, grid_spacing):
    """
    Create a heliographic overlay using wcsaxes.
    Also draw a grid and label the top axes.

    Parameters
    ----------
    axes : `~matplotlib.axes._subplots.WCSAxesSubplot` object.
        The `~astropy.visualization.wcsaxes.WCSAxes` object to create the HGS overlay on.

    grid_spacing: `~astropy.units.Quantity`
        Spacing for longitude and latitude grid in degrees.

    Returns
    -------
    overlay : `~astropy.visualization.wcsaxes.WCSAxes` overlay
        The overlay object.

    """
    lon_space = lat_space = grid_spacing
    overlay = axes.get_coords_overlay('heliographic_stonyhurst')
    lon = overlay[0]
    lat = overlay[1]
    lon.coord_wrap = 180
    lon.set_major_formatter('dd')
    lon.set_ticks_position('tr')
    lat.set_ticks_position('tr')
    grid_kw = {'color': 'white', 'zorder': 100, 'alpha': 0.5}#, linestyle: 'dashed', linewidth: 0.1}
    lon.set_ticks(spacing=lon_space, color=grid_kw['color'])
    lat.set_ticks(spacing=lat_space, color=grid_kw['color'])
    overlay.grid(**grid_kw, linestyle='dashed', linewidth=0.1)
    return overlay