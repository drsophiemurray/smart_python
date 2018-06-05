smart_python
============

*A line-by-line translation of SMART IDL. A work in progress...*

Take an SDO magnetogram and run SMART to detect magnetically interesting regions. In default mode, a .json property file and .png image of detections will be output.
Currently to run use the wrapper, e.g.,

    import wrapper
    wrapper.main(fits_file='hmi.M_720s.20170906_120000_TAI.fits')

Optional keyword inputs:
- `fits_file`: run SMART on specifically defined fits file (default will take latest available fits from JSOC)

Optional outputs:
- Save fits files of magnetogram and detections (to be added)

JSON output
-----------
The code will output a JSON file in the format `yyyymmdd_HHMM_properties.json`, which contains property values for all detected SMART regions (blank file if none), including position, magnetic, and polarity separation line information*

| Section | Subsection          | Description                                                                                    |
| :------ | :---------          | :----------                                                                                    |
| meta    | dateobs             | time of observation (yyyymmdd_HHMM)                                                            |
|         | dimension           | pixel resolution                                                                               |
|         | instrument          | observation source                                                                             |
| posprop | arid                | SMART detection no.                                                                            |
|         | x(y)bnd             | x(y) pixel coordinates of box boundary of detection                                            |
|         | x(y)cenflx(area)    | centre of detection in pixel coordinates based on flux(area)                                   |
|         | hcx(y)bnd           | heliocentric coordinate boundaries in arcseconds                                               |
|         | hcx(y)flx(area)     | heliocentric coordinates of centre of detection based on flux(area) in arcseconds              |
|         | hglat(lon)bnd       | solar heliographic latitude(longitude) coordinate boundaries                                   |
|         | hglat(lon)flx(area) | solar heliographic latitude(longitude) coordinates of centre of detection based on flux(area)  |
|         | carlonbnd           | carrington longitude of boundary                                                               |
|         | carlonflx(area)     | carrington longitude of centre of detection based on flux(area)                                |
| magprop | arid                | SMART detection no.                                                                            |
|         | areabnd             | total area of boundary box region (millionths of a solar hemisphere)                           |
|         | posareabnd          | area of positive magnetic field part of boundary box region (millionths of a solar hemisphere) |
|         | negareabnd          | area of negative magnetic field part of boundary box region (millionths of a solar hemisphere) |
|         | totarea             | total area of detection (millionths of a solar hemisphere)                                     |
|         | posarea             | area of positive magnetic field part of detection (millionths of a solar hemisphere)           |
|         | negarea             | area of negative magnetic field part of detection (millionths of a solar hemisphere)           |
|         | totflx              | total magnetic flux of detection (maxwell)                                                     |
|         | posflx              | magnetic flux of positive magnetic field part of detection (maxwell)                           |
|         | negflx              | magnetic flux of negative magnetic field part of detection (maxwell)                           |
|         | imbflx              | magnetic flux imbalance (maxwell)                                                              |
|         | frflx               | flux fraction of detection, i.e. (posflx - negflx) / totflx                                    |
|         | bmin                | total negative magnetic field strength of detection (gauss)                                    |
|         | bmax                | total positive magnetic field strength of detection (gauss)                                    |
|         | bmean               | mean magnetic field strength of detection (gauss)                                              |
| pslprop | arid                | SMART detection no.                                                                            |
|         | psllength           | polarity inversion line length of detection (megameters)                                       |
|         | pslsglength         | strong gradient PIL length (megameters)                                                        |
|         | pslcurvature        | curvature of PIL                                                                               |
|         | rvalue              | r value (maxwell)                                                                              |
|         | wlsg                | gradient-weighted integral length of PIL (gauss)                                               |
|         | bipolesep_mm        | bipole separation of detection (megameters)                                                    |
|         | bipolesep_px        | bipole separation of detection (pixels)                                                        |

*note PIL properties might have zeros since deprojection not implemented (see issue [#3](https://github.com/sophiemurray/smart_python/issues/3)).

Tracking
--------

SMART can be run in 'tracking' mode, where the evolution of SMART properties will be plotted over time since the last detection (default six hours previous). The code, developed by Sean Blake, can also be run stand alone, as long as default json and optional processed fits files are available (see optional outputs above - to be added).
More description to be added once integration is complete!

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