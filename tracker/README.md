# SMART Tracking Code

Python Script to track SMART processed sunspots, and assign likely identifying numbers. 
Main code tracker.py uses functions from SS_tracker_module.py, and plot results with plot_evolution.py

Method
---------------
 +    A list of sunspot objects are created from SMART processed .fits files for some time t. Each sunspot is given an id number. At some later time t+1, another list of sunspot numbers is created. 

 +    Every possible sunspot pair from the two lists are compared after, the older sunspots have been rotated (to estimate new position due to solar rotation).
 
 +    If the rotated sunspots overlap with any of the new sunspots, they are probably the same sunspot, and the new sunspot number is updated with the old sunspot number.

 +    If the old, rotated sunspot does not overlap with any new sunspots, it is likely retired (or not detected by the SMART program. This could be updated in the future)

 +    If there is a new sunspot which does not overlap with any of the old sunspots, it is given a new sunspot number.

General Workflow
----------------
Point this code to three folders:
 +    input_folder = contains the SMART processed .fits and properties .json files.
 +    output_folder = where updated .json files (with new id numbers for sunspots) will be written
 +    image_folder = where images will be saved

As it stands:

 +    tracking_modules.py contains the functions used by the other two scripts.
 +    tracker.py actually does the tracking. You point the script to the folder with all of the input .fits and .json files, and an output folder for .json files (see below). It will loop through these files, and figure out which sunspot numbers should be assigned to which sunspots. It then writes out a copy of each input .json file to the output .json folder, which includes a ['posprop']['trueid']  field. These are the SMART assigned sunspot numbers.
 +    plot_evolution.py makes the pretty pictures, and saves them one-by-one to an images folder (couldnt figure out how to animate them). I specify ['magprop']['trueid'] to plot in the bottom, but this should in theory be changeable for whatever property you want to plot over time.

License
-------
The content of this project is licensed under the [Creative Commons Attribution 4.0 license](https://creativecommons.org/licenses/by/4.0/), and the underlying source code used to format and display that content is licensed under the [MIT license](https://opensource.org/licenses/mit-license.php). Please reference the original Higgins et al publication in all instances and relevant GitHub repositries if using the code.