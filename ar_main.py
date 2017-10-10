
# First load the data
map = ar_readmag()

# Process it
    imgsz = len(map.data)
    original_time = map.date
    thisdatafile = original_time.strftime('%Y%m%d_%H%M')

# Now processmag
;Turn on median filtering due to the noise...
params.DOMEDIANFILT = 0
params.DOCOSMICRAY = 0

# Process magnetogram
magproc = ar_processmag(thismap, cosmap = cosmap, limbmask = limbmask, $
                        params = params, /nofilter, /nocosmicray)


