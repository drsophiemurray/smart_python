#from ar_readmag import ar_readmag
from ar_processmag import ar_processmag
from ar_detect import ar_detect
from ar_detect_core import ar_detect_core
from ar_posprop import ar_posprop
from ar_magprop import ar_magprop
from ar_pslprop import ar_pslprop
import astropy.units as u
import sunpy
import matplotlib.pylab as plt

# First load the data
#map = ar_readmag()
thismap = sunpy.map.Map('/Users/sophie/data/smart/latest.fits')
thismap = sunpy.map.Map('/Users/sophie/sunpy/data/hmi_m_45s_2011_10_17_00_01_30_tai_magnetogram.fits')

## Need to downsample if 4096x4096
if thismap.meta['naxis1'] == 4096:
    thismap = thismap.resample(u.Quantity([1024, 1024], u.pixel))

# Process it
imgsz = len(thismap.data)
original_time = thismap.date
thisdatafile = original_time.strftime('%Y%m%d_%H%M')

# Now process magnetogram
magproc, cosmap, limbmask = ar_processmag(thismap)
#TO DO- deprecated coordinates in sunpy

#Create AR masks
thissm = ar_detect(magproc, limbmask)
thisar, pslmap = ar_detect_core(magproc, thissm.data)
thismask = thisar.data
#TO DO -ridgmask! Use skeleton for now maybe?
#TO DO ar_core2mask thismask,smartmask,coresmartmask=ar_core2mask(thisar.data)

#Get properties
posprop = ar_posprop(magproc, thismask, cosmap)
magprop = ar_magprop(magproc, thismask, cosmap)
pslprop = ar_pslprop(magproc, thismask, doproj=False, projmaxscale=1024)#dproj=True, projmaxscale=1024)
#TO DO-projection!


import sunpy.map

thismap_idl = sunpy.map.Map('/Users/sophie/data/smart/latest.fits')
magproc_idl = sunpy.map.Map('magproc.fits')
thissm_idl = sunpy.map.Map('thissm.fits')
thisar_idl = sunpy.map.Map('thisar.fits')
thismaskmap_idl = sunpy.map.Map('thismaskmap.fits')
thismask_idl = thismaskmap_idl.data