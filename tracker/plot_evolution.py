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
'''

import sunpy.map
import os
from matplotlib import pyplot as plt
#import numpy as np
#from scipy import ndimage
import json
import datetime
import matplotlib.gridspec as gridspec


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

    ims = []
    fig = plt.figure(1)

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
        # plotting shite

        gs1 = gridspec.GridSpec(2, 1,
                                height_ratios=[1, 2],
                                )

        #evolution
        ax1 = plt.subplot(gs1[0])

        for key, value in property_values.items():
            ax1.plot(value[0], value[1])

        plt.axvline(date_string[1], lw = 2, linestyle = "dashed", color = "black")
        plt.title(date_string[0], fontsize = 12)

        #image
        ax2 = plt.subplot(gs1[1])

        plt.imshow(magnetogram_map.data,  origin='lower',
                   vmin=-1000., vmax=1000.,
                   cmap='Greys')
        plt.contour(detection_map.data, origin='lower',
                    colors='blue', linewidths=1.0)

        #these need to be different colours
        plt.plot(json_centx, json_centy, 'or',
                 color='blue', markersize=2)
        for x, y, numb in zip(json_centx, json_centy, number_json_values):
            plt.text(x, y, str(numb), fontsize = 14)



        count += 1
        plt.savefig(output_folder + str(count) + ".png", dpi = 100, figsize = (80, 40))
        plt.clf()
    
    # aborted attempts to get the above to animate
    #ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=1000)
    #ani.save('sdo_aia.mp4', writer='ffmpeg')



def datetime_from_file_string(a):
    """convert timestring to datetime object
    """
    year = int(a[:4])
    month = int(a[4:6])
    day = int(a[6:8])
    hour = int(a[9:11])
    time1 = datetime.datetime(year, month, day, hour)
    return time1

if __name__ == '__main__':
    main()











































