'''
    First created 2017-09-15
    Sophie A. Murray

    Python Version:    Python 3.6.1 |Anaconda custom (x86_64)| (default, May 11 2017, 13:04:09)

    Description:
    Read NRT JSOC magnetograms for use with SMART. Code adapated from MOSWOC.
    Pulling from: http://jsoc.stanford.edu/data/hmi/fits/latest_fits_time

    Notes:
    -urllib2 was removed and replaced with urlib.request for Python3 compatibility

    CHECK CODE OK FOR PYTHON 3
    ADD FITS READING
    CONFIG PARSER: config.ini
    https://wiki.python.org/moin/ConfigParserExamples
    https://docs.python.org/3/library/configparser.html
    https://community.canvaslms.com/groups/canvas-developers/blog/2016/10/18/create-a-python-config-file-for-api-scripts
    https://stackoverflow.com/questions/5055042/whats-the-best-practice-using-a-settings-file-in-python

    all pointless: from astropy.utils.data import download_file

'''

import os
import urllib.request
from configparser import ConfigParser
import sunpy.map
import astropy.units as u
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def ar_readmag():
    """Grab latest fits file time from JSOC.
    First find the latest url, download that using the fits file name.
    Then read it with SunPy and make a plot.
    """
    ## Load configuration file
    config = ConfigParser()
    config.read("config.ini")
    ## Get latest fits file
    link = get_link(config.get('input','fits_source'), int(config.get('input','strip_no')))
    data_dir = "".join([os.path.expanduser('~'), config.get('paths','data_dir')])
    fits = download(link,folder=data_dir)
    ## Read file and plot it
    map = plot_map(fits, data_dir)
    map.save(data_dir + 'latest.fits', filetype='auto')
    return map

def get_link(url,strip_no):
    """Standard URL grab, making sure the latest version is downloaded.
    Strip number included as currently there is descriptive JSOC text before the http is written.
    Also have to decode from bytes.
    """
    web_file = web_request(url)
    link = web_file.read()
    link = link.strip()[strip_no::]
    return link.decode("utf-8")

def download(url,folder):
    """Copy the contents of a file from a given URL
    to a defined folder.
    """
    web_file = web_request(url)
    ## Create with fits file name
    file_loc = "".join([folder, url.split('/')[-1]])
    if not os.path.isdir(folder):
        os.mkdir(folder)
    save_file = open(file_loc, 'wb')
    save_file.write(web_file.read())
    web_file.close()
    save_file.close()
    return file_loc

def web_request(url):
    """Set up request for weblink,
    inlcuding getting latest version (not cached)
    """
    web_file = urllib.request.Request(url)
    web_file.add_header('Cache-Control', 'max-age=0')
    web_file = urllib.request.build_opener().open(web_file)
    return web_file

def plot_map(fits, folder):
    """Quickplot of HMI data just downloaded
    Data is in map.data"""
    map = sunpy.map.Map(fits)
    ## Need to downsample if 4096x4096
    if map.meta['naxis1']==4096:
        map = map.resample(u.Quantity([1024, 1024], u.pixel))
    ## Normalise
    map.plot_settings['norm'] = colors.Normalize(vmin=-500., vmax=500.)
    #map.peek()
    # fsz = 5
    # f, (ax) = plt.subplots(1, figsize=[fsz, fsz])
    # ax.set_xlabel('X (pixels)')
    # ax.set_ylabel('Y (pixels)')
    # ax.imshow(map.data)
    # axpos = ax.get_position()
    # dpi = map.data.shape[0] / (axpos.x1 - axpos.x0) / fsz
    # plt.savefig(folder+'aia-figsave.png', dpi=dpi)
    return map


if __name__ == '__main__':
    ar_readmag()
