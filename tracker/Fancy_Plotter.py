import sunpy.map
import os
from matplotlib import pyplot as plt
import numpy as np
from scipy import ndimage
import json
import datetime

################################################################################
################################################################################

input_folder = '/Users/sophie/data/smart/track_test/'
json_folder = '/Users/sophie/data/smart/track_test_result/'
movie_folder = '/Users/sophie/data/smart/track_test_result/'

###############################################################################
def main():

    filenames = sorted(os.listdir(json_folder))
    filename_dates = [datetime_from_file_string(x) for x in filenames]

    #cuts off last day? dunno whats going on...
    start_date, end_date = filename_dates[0],  filename_dates[len(filename_dates)-1] #datetime.datetime(2011, 9, 1, 0, 0, 0, 0), datetime.datetime(2012, 9, 4, 1, 1)

    date_strings = []
    for index, value in enumerate(filename_dates):
        if start_date < value < end_date:
            date_strings.append([filenames[index][:13], value])

    # get properties
    sunspot_data = {}
    for date_string in date_strings:
        json_filename = json_folder + date_string[0] + "_properties.json"
        data = json.load(open(json_filename))

        for key, value in data['posprop']['trueid'].items():
            if str(value) in sunspot_data:
                sunspot_data[str(value)][0].append(date_string[1])
                sunspot_data[str(value)][1].append(data['magprop']['totarea'][str(key)])
            else:
                sunspot_data[str(value)] = [[date_string[1]], [data['magprop']['totarea'][str(key)]]]


    #----------------------------------------
    # get detection outlines as outline_edges

    ims = []
    fig = plt.figure(1)

    count = 1
    for date_string in date_strings:
        filename = input_folder + date_string[0] + "_detections.fits"
        yy1 = sunpy.map.Map(filename).data

        # get outlines of sunspot detections
        outline_edges = np.zeros((1024, 1024))
        num_of_ss = np.max(yy1.flatten())

        for i in np.arange(1, num_of_ss + 1):
            yy_copy = np.copy(yy1)

            mask = yy_copy==i
            yy_copy[~mask] = 0

            # rank 2 structure with full connectivity
            struct = ndimage.generate_binary_structure(2, 2)
            erode = ndimage.binary_erosion(mask, struct)
            edges = mask ^ erode

            outline_edges += edges
        outline_edges = np.ma.masked_where(outline_edges == 0, outline_edges)

        #----------------------------------------
        # get actual image of sun
        filename = input_folder + date_string[0] + "_map.fits"
        yy2 = sunpy.map.Map(filename).data

        #----------------------------------------
        # read in numbers and centroids from json
        # read in json
        json_data = json.load(open(json_folder + date_string[0] + "_properties.json"))
        time2 = datetime_from_json(json_data)
        number_json = list(json_data['posprop']['trueid'].keys())
        number_json_values = [json_data['posprop']['trueid'][i] for i in number_json]


        json_centx, json_centy = [], []
        for i in number_json:
            json_centx.append(json_data['posprop']['xcenarea'][i])
            json_centy.append(json_data['posprop']['ycenarea'][i])

        #----------------------------------------
        # plotting shite

        ax1 = plt.subplot2grid((5, 4), (0, 0), colspan = 4, rowspan = 3)
        im1 = plt.imshow(yy2)
        im2 = plt.imshow(outline_edges, cmap = "Greys", interpolation = "nearest")

        plt.plot(json_centx, json_centy, 'or')
        for x, y, numb in zip(json_centx, json_centy, number_json_values):
            plt.text(x, y, str(numb), fontsize = 24, color = "red")

        plt.title(date_string[0], fontsize = 24)

        ax2 = plt.subplot2grid((5, 4), (3, 0), colspan = 4, rowspan = 2)

        for key, value in sunspot_data.items():
            ax2.plot(value[0], value[1])

        plt.axvline(date_string[1], lw = 3, linestyle = "dashed", color = "black")

        count += 1
        plt.savefig(movie_folder + str(count) + ".png", dpi = 100, figsize = (80, 40))
        plt.clf()
    
    # aborted attempts to get the above to animate
    #ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=1000)
    #ani.save('sdo_aia.mp4', writer='ffmpeg')



def datetime_from_json(data):
    # convert timestring to datetime object

    a = data['meta']['dateobs']
    year = int(a[:4])
    month = int(a[4:6])
    day = int(a[6:8])
    hour = int(a[9:11])
    time1 = datetime.datetime(year, month, day, hour)
    return time1

def datetime_from_file_string(a):
    # convert timestring to datetime object

    year = int(a[:4])
    month = int(a[4:6])
    day = int(a[6:8])
    hour = int(a[9:11])
    time1 = datetime.datetime(year, month, day, hour)
    return time1

if __name__ == '__main__':
    main()











































