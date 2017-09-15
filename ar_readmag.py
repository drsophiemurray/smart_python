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
'''

import os
import urllib.request
from configparser import ConfigParser

def main():
    """Grab latest fits file time from JSOC.
    """
    config = ConfigParser()
    config.read("config.ini")
    link = get_link(config.get('input','fits_source'), int(config.get('input','strip_no')))
    download(link,folder=config.get('paths','data_dir'))

def get_link(url,strip_no):
    """Standard URL grab, making sure the latest version is downloaded.
    Strip number included as currently there is descriptive text before the http is written.
    """
    web_file = urllib.request.Request(url)
    web_file.add_header('Cache-Control', 'max-age=0')
    web_file = urllib.request.build_opener().open(web_file)
    link = web_file.read()
    return link.strip()[strip_no::]

def download(url,folder):
    """Copy the contents of a file from a given URL.
    """
    web_file = urllib.Request(url)
    web_file.add_header('Cache-Control', 'max-age=0')
    web_file = urllib.build_opener().open(web_file)
    folder = "".join([os.path.expanduser('~'), OUT_FOLDER])
    file_loc = "".join([os.path.expanduser('~'),
                        OUT_FOLDER, url.split('/')[-1]])
    if not os.path.isdir(folder):
        os.mkdir(folder)
    save_file = open(folder, 'w')
    save_file.write(web_file.read())
    web_file.close()
    save_file.close()
    del folder
    return file_loc


if __name__ == '__main__':
    main()
