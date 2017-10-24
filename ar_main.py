from ar_readmag import ar_readmag
from ar_processmag import ar_processmag
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

thissm = ar_detect(magproc, params = params, status = smartstatus, $
                        cosmap = cosmap, limbmask = limbmask) ;prob dont need this?


