#from ar_readmag import ar_readmag
from ar_processmag import ar_processmag
from ar_detect import ar_detect
import sunpy

# First load the data
#map = ar_readmag()
map = sunpy.map.Map('/Users/sophie/data/smart/latest.fits')

# Process it
imgsz = len(map.data)
original_time = map.date
thisdatafile = original_time.strftime('%Y%m%d_%H%M')

# Now process magnetogram
magproc, cosmap, limbmask = ar_processmag(map)

#Create AR masks

thissm = ar_detect(magproc, limbmask)
thisar, pslmap = ar_detect_core(magproc, thissm.data)

import matplotlib.pylab as plt
plt.ion()
plt.imshow(thissm.data)