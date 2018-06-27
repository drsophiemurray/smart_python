'''
    Solar Monitor Active Region Tracker Visualation
    ===============================================
    Written by Sean Blake, with help from Tadhg Garton.
    Code later integrated into smart_python Github repository by Sophie Murray

    Provide a folder location that contains all the files you want to run, and a property you want to plot.
    The code will then make some evolution plots of the tracked regions (images and 2D).

    Last tested under:
    - Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)
    - Sunpy 0.9.0

    Steps:
    - Load .json file with true ids, and detections and maps
    - Plot images
    - Plot properties over time

    TODO:
    - Sort out the reading in of files - repetitive and uses hardcoded numbers.
    - Naming scheme, one giant function.
    - Plot colours for the numbers.
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


def main(input_folder='/Users/sophie/data/smart/track_test/',
         group='magprop', property='totarea',
         *smart_folder, **output_folder):
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


    # get properties (time and value for each id)
    property_values = {}
    x_position, y_position = {}, {}
    for date_string in date_strings:
        json_filename = input_folder + date_string[0] + "_properties.json"
        json_data = json.load(open(json_filename))
        for key, value in json_data['posprop']['trueid'].items():
            if str(value) in property_values:
                property_values[str(value)][0].append(date_string[1])
                property_values[str(value)][1].append(json_data[group][property][str(key)])
            else:
                property_values[str(value)] = [[date_string[1]], [json_data[group][property][str(key)]]]
            if str(value) in x_position:
                x_position[str(value)].append(json_data['posprop']['xcenarea'][str(key)])
                y_position[str(value)].append(json_data['posprop']['ycenarea'][str(key)])
            else:
                x_position[str(value)] = [json_data['posprop']['xcenarea'][str(key)]]
                y_position[str(value)] = [json_data['posprop']['ycenarea'][str(key)]]

    #----------------------------------------
    # get detection outlines as outline_edges

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
        # plot evolution of property
        fig = plt.figure(figsize=(10, 12))
        ax1 = fig.add_subplot(2, 1, 2)

        colors = itertools.cycle(["black", "grey",
                                  "brown", "orange",
                                  "red", "pink",
                                  "purple", "blue",
                                  "turquoise", "green"])

        for key, value in property_values.items():
            ax1.plot(value[0], value[1],
                     label=key,
                     marker= 'o', markersize=3.0)
            plt.legend(loc='upper left')

        plt.gcf().autofmt_xdate()
        plt.axvline(date_string[1],
                    linestyle = "dashed", color = "black")
        plt.ylabel('Total Area [m.s.h]') #todo should be meta data
#        plt.xlabel('Date time') #todo should be meta data

        # plot detections on sunpy map of magnetogram
        ax1 = fig.add_subplot(2, 1, 1, projection=magnetogram_map)
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95,
                            wspace=None, hspace=None)

        # Get same axes
        bottom_left = SkyCoord(-1000 * u.arcsec, -1000 * u.arcsec,
                               frame=magnetogram_map.coordinate_frame)
        top_right = SkyCoord(1000 * u.arcsec, 1000 * u.arcsec,
                             frame=magnetogram_map.coordinate_frame)
        submap = magnetogram_map.submap(bottom_left, top_right)
        axes = wcsaxes_compat.gca_wcs(magnetogram_map.wcs)

        image = magnetogram_map.plot(vmin=-500, vmax=500, axes=axes)

        # Draw solar lat/lon grid
        axes.coords.grid(False)
        overlay = grid_overlay(axes, grid_spacing=10 * u.deg)
#        plt.colorbar(label='B [G]')

        plt.contour(detection_map.data, origin='lower',
                    colors='lightblue', linewidths=0.5)

        # add numbers
        plt.plot(json_centx, json_centy, 'or',
                 color='yellow', markersize=2.0)
        for x, y, numb in zip(json_centx, json_centy, number_json_values):
            plt.text(x+10, y+10, str(numb),
                     color='yellow')
        plt.title(date_string[0])

        plt.savefig(output_folder + date_string[0] + "_tracking.png",
                    dpi=100)
        plt.close()

    # convert to gif
    images = []
    filenames = sorted(os.listdir(input_folder))
    filenames_images = [x for x in filenames if "_tracking.png" in x]
    filenames_images = [input_folder + x for x in filenames_images]
    for filename in filenames_images:
        images.append(imageio.imread(filename))
    imageio.mimsave(input_folder+'SMART_evolution.gif',
                    images,
                    fps=1.)

def datetime_from_file_string(a):
    """convert timestring to datetime object
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











































