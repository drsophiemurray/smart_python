# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:09:53 2018

@author: laura
"""

import sunpy.map
import os
from matplotlib import pyplot as plt
import json
import datetime
import itertools
import matplotlib as mpl
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunpy.visualization import wcsaxes_compat
import imageio
import detect_core
import numpy as np


#plotting labels on SMART images

'''
This code takes in outputted json files, maps and detections and outputs images labelled
by Area IDs. 

'''

def main(input_folder='C:/Users/laura/Documents/SMART/labelled_data/',
         output_folder='C:/Users/laura/Documents/SMART/labelled_data/'):
    
    
    filenames = sorted(os.listdir(input_folder))
    filenames_json = [x for x in filenames if ".json" in x]
    
    
    
    filename_dates = [datetime_from_file_string(x) for x in filenames_json]
    
    
    start_date, end_date = filename_dates[0],  filename_dates[len(filename_dates)-1]

    date_strings = []
    for index, value in enumerate(filename_dates):
        if start_date <= value <= end_date:
            date_strings.append([filenames_json[index][:13], value])
    
    print(date_strings)
    
    property='totarea'
    group='magprop'
    property_values_totarea = {}
    x_position, y_position = {}, {}
    for date_string in date_strings:
        json_filename = input_folder + date_string[0] + "_properties.json"
        json_data = json.load(open(json_filename))
        for key, value in json_data['posprop']['arid'].items():
            if str(value) in property_values_totarea:
                property_values_totarea[str(value)][0].append(date_string[1])
                property_values_totarea[str(value)][1].append(json_data[group][property][str(key)])
            else:
                property_values_totarea[str(value)] = [[date_string[1]], [json_data[group][property][str(key)]]]
            if str(value) in x_position:
                x_position[str(value)].append(json_data['posprop']['xcenarea'][str(key)])
                y_position[str(value)].append(json_data['posprop']['ycenarea'][str(key)])
            else:
                x_position[str(value)] = [json_data['posprop']['xcenarea'][str(key)]]
                y_position[str(value)] = [json_data['posprop']['ycenarea'][str(key)]]
    
    
    
    
    for date_string in date_strings:
        detection_filename = input_folder + date_string[0] + "_detections.fits"
        detection_map = sunpy.map.Map(detection_filename)

        #----------------------------------------
        # get actual image of sun
        magnetogram_filename = input_folder + date_string[0] + "_map.fits"
        magnetogram_map = sunpy.map.Map(magnetogram_filename)
        
        psl_filename= input_folder + date_string[0] + "_psl.fits"
        psl_map = sunpy.map.Map(psl_filename)

        #----------------------------------------
        # read in numbers and centroids from json
        json_data = json.load(open(input_folder + date_string[0] + "_properties.json"))
        # smart id
        number_json = list(json_data['posprop']['arid'].keys())
        # the tracked ids in the data
        number_json_values = [json_data['posprop']['arid'][i] for i in number_json]

        json_centx, json_centy = [], []
        for i in number_json:
            json_centx.append(json_data['posprop']['xcenarea'][i])
            json_centy.append(json_data['posprop']['ycenarea'][i])

        #----------------------------------------
        
        fig = plt.figure()
        
        bottom_left = SkyCoord(-1000*u.arcsec, -1000*u.arcsec,
                           frame=magnetogram_map.coordinate_frame)
        top_right = SkyCoord(1000*u.arcsec, 1000*u.arcsec,
                         frame=magnetogram_map.coordinate_frame)
        submap = magnetogram_map.submap(bottom_left, top_right)
        axes = wcsaxes_compat.gca_wcs(magnetogram_map.wcs)
        image = magnetogram_map.plot(vmin=-500, vmax=500, axes=axes)
        axes.coords.grid(False)
        
    # Draw solar lat/lon grid
        overlay = grid_overlay(axes, grid_spacing=10 * u.deg)
        plt.colorbar(label='B [G]')
    
        plt.contour(psl_map.data, origin='lower',
                colors='blue', linewidths=0.5,
                vmin=0., vmax=np.max(np.unique(detection_map.data))+1)
   
        plt.contour(detection_map.data, origin='lower',
                    colors='red', linewidths=0.5)
#                
#     #Add number labels
        for x, y, numb in zip(json_centx, json_centy, number_json_values):
            plt.text(x+10, y+10, str(numb),
                     color='yellow')
        


        plt.savefig(output_folder + date_string[0] + "_tracking.png",
                    dpi=1000)
        plt.close()
        
    
    
    
    
    
def datetime_from_file_string(a):
    """
    Extracting a datetime object from the SMART json file
    - convert timestring to datetime object
    """
    year = int(a[:4])
    month = int(a[4:6])
    day = int(a[6:8])
    hour = int(a[9:11])
    time1 = datetime.datetime(year, month, day, hour)
    return time1

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

if __name__ == '__main__':
    main()
