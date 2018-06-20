import json
import datetime
import os
import sunpy.map
from copy import deepcopy
import SS_tracker_module as SS

###############################################################################
################################################################################

input_folder = '/Users/sophie/data/smart/track_test/'
output_folder = '/Users/sophie/data/smart/track_test_result/'

###############################################################################

def main():
    # first image
    # read in json
    filenames = sorted(os.listdir(input_folder))
    filename_stems = [x[:14] for x in filenames if "_detections.fits" in x]
    extensions = ["properties.json", "detections.fits", "map.fits"]

    json_data = json.load(open(input_folder + filename_stems[0] + extensions[0]))
    time1 = datetime_from_json(json_data)

    # read in map
    yy = sunpy.map.Map(input_folder + filename_stems[0] + extensions[1]).data
    master_num_of_ss, old_SS = SS.get_sunspot_data(yy, time1)

    # write out json with updated 'trueid' numbers
    true_id = {}
    for index in range(len(old_SS)):
        true_id[str(index)] = int(old_SS[index].number)

    json_data['posprop']['trueid'] = true_id
    with open(output_folder + filename_stems[0] + extensions[0], 'w') as outfile:
        json.dump(json_data, outfile)

    ###############################################################################
    # now for the rest of the images
    count = 0
    for filename_stem in filename_stems[1:]:
        # read in json
        json_data = json.load(open(input_folder + filename_stem + extensions[0]))
        time2 = datetime_from_json(json_data)

        # read in map
        yy2 = sunpy.map.Map(input_folder + filename_stem + extensions[1]).data
        num_of_ss2, new_SS = SS.get_sunspot_data(yy2, time2)

        overlap_matrix = SS.make_overlap_matrix_V2(old_SS, new_SS)
        old_SS, new_SS, master_num_of_ss = SS.assign_numbers(old_SS, new_SS, overlap_matrix, master_num_of_ss)

        # write out json with updated 'trueid' numbers
        true_id = {}
        for index in range(len(new_SS)):
            true_id[str(index)] = int(new_SS[index].number)

        json_data['posprop']['trueid'] = true_id
        with open(output_folder + filename_stem + extensions[0], 'w') as outfile:
            json.dump(json_data, outfile)

        old_SS = deepcopy(new_SS)

        count += 1
        print(count)


def datetime_from_json(data):
    """convert timestring to datetime object
    """
    a = data['meta']['dateobs']
    year = int(a[:4])
    month = int(a[4:6])
    day = int(a[6:8])
    hour = int(a[9:11])
    time1 = datetime.datetime(year, month, day, hour)
    return time1



if __name__ == '__main__':
    main()







