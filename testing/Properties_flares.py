# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 13:47:57 2018

@author: laura
"""


'''
This code overplots flare data from a flare_list file to tracked properties plots
'''

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


mpl.rc('font', size = 10, family = 'serif', weight='normal')
mpl.rc('legend', fontsize = 8)
mpl.rc('lines', linewidth = 1.5)


def main(input_folder='C:/Users/laura/Documents/SMART/tracking/output/',
         
         smart_folder='C:/Users/laura/Documents/SMART/tracking/input/', 

         output_folder='C:/Users/laura/Documents/SMART/tracking/properties_w_flares/'):
    """
    Parameters
    ----------
    input_folder : Folder string location of .json files with true_id's
    group : Indicate what group in .json file the property you wish to plot is, e.g. 'magprop'
    property : Indicate what property you wish to plot, e.g. 'totarea'
    smart_folder : Optional folder location to load maps and detections from if not in same folder as json
    output_folder: Optional output folder location of images created with algorithm
    Returns
    -------
    """
#    if not input_folder:
#        input_folder = os.getcwd() + '/'
    if not smart_folder:
        smart_folder = input_folder
    if not output_folder:
        output_folder = input_folder

    # load json files
    filenames = sorted(os.listdir(input_folder))
    filenames_json = [x for x in filenames if ".json" in x]

    filename_dates = [datetime_from_file_string(x) for x in filenames_json]
    start_date, end_date = filename_dates[0],  filename_dates[len(filename_dates)-1]

    date_strings = []
    for index, value in enumerate(filename_dates):
        if start_date <= value <= end_date:
            date_strings.append([filenames_json[index][:13], value]) #todo get rid of hardcoded 13
            
    #getting flare times and dates
    file_location='C:/Users/laura/Documents/SMART/flare_list/mx_list.txt'

    f=open(file_location,"r")
    lines=f.readlines()
    flaredatesM=[]
    flaredatesX=[]
    for x in lines:
        flareclass=(x.split(' ')[6])
        if flareclass[0]=='M':    #Pulling out the M Class flares 
            yearM=(x.split(' ')[0])
            monthM=(x.split(' ')[1])
            dayM=(x.split(' ')[2])
            hourM=(x.split(' ')[3])
    
            timeM = yearM+monthM+dayM+hourM
             
#        
            yearM2=(timeM[:4])
            monthM2=(timeM[4:6])
            dayM2=(timeM[6:8])
            hourM2=(timeM[8:10])
            minuteM2=(timeM[10:12])
#        
#        
            flaretimeM=str(yearM2)+'-'+str(monthM2)+'-'+str(dayM2)+' '+str(hourM2)+':'+str(minuteM2)
            flaredatesM.append(flaretimeM)
            
   
        if flareclass[0]=='X':    #Pulling out the X Class flares
            yearX=(x.split(' ')[0])
            monthX=(x.split(' ')[1])
            dayX=(x.split(' ')[2])
            hourX=(x.split(' ')[3])
    
            timeX = yearX+monthX+dayX+hourX
#        
#        
            yearX2=(timeX[:4])
            monthX2=(timeX[4:6])
            dayX2=(timeX[6:8])
            hourX2=(timeX[8:10])
            minuteX2=(timeX[10:12])
#        
#        
            flaretimeX=str(yearX2)+'-'+str(monthX2)+'-'+str(dayX2)+' '+str(hourX2)+':'+str(minuteX2)
            flaredatesX.append(flaretimeX)
   
    print(flaredatesX)
        

        
    #get properties (time and value for each ID)
    #property=totarea
    property='totarea'
    group='magprop'
    property_values_totarea = {}
    x_position, y_position = {}, {}
    for date_string in date_strings:
        json_filename = input_folder + date_string[0] + "_properties.json"
        json_data = json.load(open(json_filename))
        for key, value in json_data['posprop']['trueid'].items():
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

    #property=bmax
    property='bmax'
    group='magprop'
    property_values_bmax = {}
    x_position, y_position = {}, {}
    for date_string in date_strings:
        json_filename = input_folder + date_string[0] + "_properties.json"
        json_data = json.load(open(json_filename))
        for key, value in json_data['posprop']['trueid'].items():
            if str(value) in property_values_bmax:
                property_values_bmax[str(value)][0].append(date_string[1])
                property_values_bmax[str(value)][1].append(json_data[group][property][str(key)])
            else:
                property_values_bmax[str(value)] = [[date_string[1]], [json_data[group][property][str(key)]]]
            if str(value) in x_position:
                x_position[str(value)].append(json_data['posprop']['xcenarea'][str(key)])
                y_position[str(value)].append(json_data['posprop']['ycenarea'][str(key)])
            else:
                x_position[str(value)] = [json_data['posprop']['xcenarea'][str(key)]]
                y_position[str(value)] = [json_data['posprop']['ycenarea'][str(key)]]
                
#property=psllength
    property='rvalue'
    group='pslprop'
    property_values_psllength = {}
    x_position, y_position = {}, {}
    for date_string in date_strings:
        json_filename = input_folder + date_string[0] + "_properties.json"
        json_data = json.load(open(json_filename))
        for key, value in json_data['posprop']['trueid'].items():
            if str(value) in property_values_psllength:
                property_values_psllength[str(value)][0].append(date_string[1])
                property_values_psllength[str(value)][1].append(json_data[group][property][str(key)])
            else:
                property_values_psllength[str(value)] = [[date_string[1]], [json_data[group][property][str(key)]]]
            if str(value) in x_position:
                x_position[str(value)].append(json_data['posprop']['xcenarea'][str(key)])
                y_position[str(value)].append(json_data['posprop']['ycenarea'][str(key)])
            else:
                x_position[str(value)] = [json_data['posprop']['xcenarea'][str(key)]]
                y_position[str(value)] = [json_data['posprop']['ycenarea'][str(key)]]
                
                

    #----------------------------------------
     #get detection outlines as outline_edges

    count = 1
    for date_string in date_strings:
        detection_filename = smart_folder + date_string[0] + "_detections.fits"
        detection_map = sunpy.map.Map(detection_filename)

        #----------------------------------------
        # get actual image of sun
        magnetogram_filename = smart_folder + date_string[0] + "_map.fits"
        magnetogram_map = sunpy.map.Map(magnetogram_filename)

        #----------------------------------------
        # read in numbers and centroids from json
        json_data = json.load(open(input_folder + date_string[0] + "_properties.json"))
        # smart id
        number_json = list(json_data['posprop']['trueid'].keys())
        # the tracked ids in the data
        number_json_values = [json_data['posprop']['trueid'][i] for i in number_json]

        json_centx, json_centy = [], []
        for i in number_json:
            json_centx.append(json_data['posprop']['xcenarea'][i])
            json_centy.append(json_data['posprop']['ycenarea'][i])

        #----------------------------------------
        
        fig = plt.figure(figsize=(10, 12))
        
        ax1 = fig.add_subplot(3, 1, 3)
        
        #plotting evolution of psllength
        
        for key, value in property_values_psllength.items():
            if key =='8':
                ax1.plot(value[0], value[1],
                     label=key,
                     marker= 'o', markersize=3.0)
                plt.legend(loc='upper left')

        import matplotlib.dates as mdates
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H'))
        plt.gcf().autofmt_xdate()
        plt.xlabel('Date and time [UT]')

        plt.ylabel('R-value [Maxwell]')
        
        plt.title('R-value v Time')
        
        for date in flaredatesX:
                plt.axvline(date,
                            linestyle = "dashed", color = "red", label= "X Flare" if date=='2017-09-10 15:35' else "")
    
                
        for date in flaredatesM:
                plt.axvline(date,
                            linestyle = "dashed", color = "purple", label= "M Flare" if date=='2017-09-09 22:04' else "")
        
                plt.legend()
        
        #plot evolution of totarea
        ax1 = fig.add_subplot(3, 1, 2)


        for key, value in property_values_totarea.items():
            if key =='8':
                ax1.plot(value[0], value[1],
                     label=key,
                     marker= 'o', markersize=3.0)
                plt.legend(loc='upper left')

        import matplotlib.dates as mdates
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H'))
        plt.gcf().autofmt_xdate()
        plt.xlabel('Date and time [UT]')

        plt.ylabel('Total Area [M.S.H]') #todo should be meta data
        
        plt.title('Total Area v Time')
        
        for date in flaredatesX:
                plt.axvline(date,
                            linestyle = "dashed", color = "red", label= "X Flare" if date=='2017-09-10 15:35' else "")
    
                
        for date in flaredatesM:
                plt.axvline(date,
                            linestyle = "dashed", color = "purple", label= "M Flare" if date=='2017-09-09 22:04' else "")
        
                plt.legend()
        
        
        #plot evolution of bmax
        ax1 = fig.add_subplot(3, 1, 1)

        for key, value in property_values_bmax.items():
            if key =='8':
                ax1.plot(value[0], value[1],
                     label=key,
                     marker= 'o', markersize=3.0)
                plt.legend(loc='upper left')

        import matplotlib.dates as mdates
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H'))
        plt.gcf().autofmt_xdate()

        plt.ylabel('Maximum Magnetic Field [G]') #todo should be meta data
        

#        # plot detections on sunpy map of magnetogram
#        ax1 = fig.add_subplot(2, 1, 1, projection=magnetogram_map)
#        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95,
#                            wspace=None, hspace=None)
#
#        # Get same axes
#        bottom_left = SkyCoord(-1000 * u.arcsec, -1000 * u.arcsec,
#                               frame=magnetogram_map.coordinate_frame)
#        top_right = SkyCoord(1000 * u.arcsec, 1000 * u.arcsec,
#                             frame=magnetogram_map.coordinate_frame)
#        submap = magnetogram_map.submap(bottom_left, top_right)
#        axes = wcsaxes_compat.gca_wcs(magnetogram_map.wcs)
#
#        image = magnetogram_map.plot(vmin=-500, vmax=500, axes=axes)
#
#        # Draw solar lat/lon grid
#        axes.coords.grid(False)
#        overlay = grid_overlay(axes, grid_spacing=10 * u.deg)
##        plt.colorbar(label='B [G]')
#   
#        plt.contour(detection_map.data, origin='lower',
#                    colors='lightblue', linewidths=0.5)
#        
#        
#                
#        # add numbers
#        plt.plot(json_centx, json_centy, 'or',
#                 color='yellow', markersize=2.0)
#        for x, y, numb in zip(json_centx, json_centy, number_json_values):
#            plt.text(x+10, y+10, str(numb),
#                     color='yellow')
        plt.title('Maximum Magnetic Field v Time')
        
        for date in flaredatesX:
                plt.axvline(date,
                            linestyle = "dashed", color = "red", label= "X Flare" if date=='2017-09-10 15:35' else "")
    
                
        for date in flaredatesM:
                plt.axvline(date,
                            linestyle = "dashed", color = "purple", label= "M Flare" if date=='2017-09-09 22:04' else "")
        
                plt.legend()
                     
        
        plt.savefig(output_folder + date_string[0] + "_tracking.png",
                    dpi=150)
        plt.close()
#
#    # convert to gif
#    images = []
#    filenames = sorted(os.listdir(input_folder))
#    filenames_images = [x for x in filenames if "_tracking.png" in x]
#    filenames_images = [input_folder + x for x in filenames_images]
#    for filename in filenames_images:
#        images.append(imageio.imread(filename))
#    imageio.mimwrite(input_folder+'SMART_evolution.gif',
#                    images,
#                    fps=1.)

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
