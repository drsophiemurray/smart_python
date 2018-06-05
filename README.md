smart_python
============

*A line-by-line translation of SMART IDL. Not yet prettified...*

Take an SDO magnetogram and run SMART to detect magnetically interesting regions. In default mode, a .json property file and .png image of detections will be output.
Currently to run use the wrapper, e.g.,

    import wrapper
    wrapper.main(fits_file='hmi.M_720s.20170906_120000_TAI.fits')

Optional keyword inputs:
- `fits_file` run SMART on specifically defined fits file (default will take latest available fits from JSOC)


JSON output
-----------
yyyymmdd_HHMM_properties.json contains SMART detection properties, including position, magnetic, and PIL properties* 

| Section | Subsection | Description |
| :------ | :--------- | :---------- |
| meta | dateobs | time of observation | 
|  | dimension | pixel resolution |
|  | instrument | observation source |
| posprop | arid | SMART detection no. |
|  | x(y)bnd | x(y) pixel coordinates of box boundary of detection |
|  | x(y)cenflx(area) | centre of detection in pixel coordinates based on flux(area) |
|  | hcx(y)bnd | heliocentric coordinate boundaries in arcseconds |
|  | hcx(y)flx(area) | heliocentric coordinates of centre of detection based on flux(area) in arcseconds |
|  | hglat(lon)bnd | solar heliographic latitude(longitude) coordinate boundaries |
|  | hglat(lon)flx(area) | solar heliographic latitude(longitude) coordinates of centre of detection based on flux(area) |
|  | carlonbnd | carrington longitude of boundary |
|  | carlonflx(area) | carrington longitude of centre of detection based on flux(area) |
| magprop | arid | SMART detection no. |
|  | areabnd | total area of boundary box region |
|  | posareabnd | area of positive magnetic field part of boundary box region |
|  | negareabnd | area of negative magnetic field part of boundary box region |
|  | totarea | total area of detection (units are millionths of a solar hemisphere) |
|  | posarea | area of positive magnetic field part of detection |
|  | negarea | area of negative magnetic field part of detection |
|  | totflx | total magnetic flux of detection (units are maxwell) |
|  | posflx | magnetic flux of positive magnetic field part of detection |
|  | negflx | magnetic flux of negative magnetic field part of detection |
|  | imbflx | |
|  | frflx | flux fraction of detection, i.e. (posflx - negflx) / totflx |
|  | bmin | total negative magnetic field strength of detection |
|  | bmax | total positive magnetic field strength of detection |
|  | bmean | mean magnetic field strength of detection (units are gauss) |
| pslprop | arid | SMART detection no. |
|  | psllength | polarity inversion line length of detection (units are megameters) |
|  | pslsglength | strong gradient PIL length (megameters) |
|  | pslcurvature | curvature of PIL |
|  | rvalue | r value (maxwell) |
|  | wlsg | gradient-weighted integral length of PIL (gauss) |
|  | bipolesep_mm | bipole separation of detection (units are megameters) |
|  | bipolesep_px | bipole separation of detection (units are pixels) |

*note PIL properties might have zeros since deprojection not implemented (see issue [#3](https://github.com/sophiemurray/smart_python/issues/3)).

External dependencies
---------------------

    astropy
    matplotlib
    numpy
    pandas
    scikit_image
    scipy
    sunpy

See `requirements.txt` for latest development versions created with pipreqs.

License
-------
The content of this project is licensed under the [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/), and the underlying source code used to format and display that content is licensed under the [MIT license](https://opensource.org/licenses/mit-license.php). Please reference the original Higgins et al publication in all instances and relevant GitHub repositries if using the code.