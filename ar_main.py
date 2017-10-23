from ar_readmag import ar_readmag
from ar_processmag import ar_processmag


# First load the data
map = ar_readmag()

# Process it
imgsz = len(map.data)
original_time = map.date
thisdatafile = original_time.strftime('%Y%m%d_%H%M')

# Now process magnetogram
magproc cosmap, limbmask = ar_processmag(map)


