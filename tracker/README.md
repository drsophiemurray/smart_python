# SMART_Python
My unfinished contribution to the sunspot tracking dealy...

As it stands:

(1) tracking_modules.py contains the functions used by the other two scripts.

(2) tracker.py actually does the tracking. You point the script to the folder with all of the input .fits and .json files, and an output folder for .json files (see below). 

It will loop through these files, and figure out which sunspot numbers should be assigned to which sunspots. It then writes out a copy of each input .json file to the output .json folder, which includes a ['posprop']['trueid']  field. These are the SMART assigned sunspot numbers.

(3) plot_evolution.py makes the pretty pictures, and saves them one-by-one to an images folder (couldnt figure out how to animate them). I specify ['magprop']['trueid'] to plot in the bottom, but this should in theory be changeable for whatever property you want to plot over time.
