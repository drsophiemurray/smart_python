# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:48:50 2018

@author: laura
"""

'''
Created on Jul 23, 2015

@author: smurray

This code grabs M and X flare data from a list of flare data. 

Python Version:    2.7.2 (default, Oct  1 2012, 15:56:20)
Working directory:     /home/smurray/workspace/verification

Extract XRA rows from GOES event lists into a nice observations text file.

***THIS NEEDS TO BE A CRON JOB - TO BE DONE!***

http://legacy-www.swpc.noaa.gov/ftpdir/warehouse/2015/2015_events/
http://legacy-www.swpc.noaa.gov/ftpdir/warehouse/2016/2016_events/
#YYYY MM DD Event    Begin    Max       End  Obs  Q  Type  Loc/Frq   Particulars       Reg#
   0  1   2   3        4       5         6    7   8    9     10        11         12    13
'''


DIR_TXT = "C:/Users/laura/Documents/SMART/swpc_event_data/"
GRAB = "XRA"

import os
import datetime 

def main():
    """Grab list of M and X class flare events from the SWPC event lists
       First extract all X ray events, then be selective..."""
    flare_list = []
    for file in os.listdir(DIR_TXT): 
        # Ignore hidden files      
        if file in ['.access^','.txt~']:  # also need to ignore .txt~ -> maybe glob was better!
            print('Skipping hidden file')
        else:
            # Extract date info from file name
            filename = file.strip()
            year = str(filename)[0:4]
            month = str(filename)[4:6]
            day = str(filename)[6:8]
            print(year+month+day)
            try:
                with open(DIR_TXT + file, "r") as inp:
                    for line in inp:
                        #if any(s in line for s in GRAB):
                        if GRAB in line:
                            todays_flares = year + ' ' + month + ' ' + day + ' ' + line[0:80]   #avoiding \n at end of these lines
                            flare_list.append(todays_flares)
            except IOError:
                print("Cannot find XRA", file)
    # Order in terms of date
    sorted_flare_list = sorted(flare_list, key=lambda x: datetime.datetime.strptime(x.strip()[0:10], '%Y %m %d'))
    # Save list
    outfile = open('C:/Users/laura/Documents/SMART/flare_list/flare_list.txt', 'w')
    for item in sorted_flare_list:
        outfile.write("%s\n" % item)
 
    # Now lets just get the M and X class flares       
    #   f = open("/home/h05/smurray/flare_list.txt", 'r')   #loadtxt(comments = "#")
    #    x = f.readlines() # then split each element  
    mx_list = []
    for line in sorted_flare_list:
        split_line = line.split()
        # Extract the flares to create a list that looks like:
        # YYYY MM DD Start Peak End Magnitude Region
        # If no region specified, list as NaN to edit later
        if '+' in split_line:
            if 'M' in split_line[12] or 'X' in split_line[12]:
                if len(split_line) >=15 :
                    todays_mx = (split_line[0] + ' ' + split_line[1] + ' ' + split_line[2] + ' ' + 
                                 split_line[5] + ' ' + split_line[6] + ' ' + split_line[7] + ' ' + 
                                 split_line[12] + ' ' + split_line[14])
                    mx_list.append(todays_mx)           
                else:
                    todays_mx = (split_line[0] + ' ' + split_line[1] + ' ' + split_line[2] + ' ' + 
                                 split_line[5] + ' ' + split_line[6] + ' ' + split_line[7] + ' ' + 
                                 split_line[12] + ' ' + 'NaN')
                    mx_list.append(todays_mx)           
        else:
            if 'M' in split_line[11] or 'X' in split_line[11]:
                if len(split_line) >=14 :
                    todays_mx = (split_line[0] + ' ' + split_line[1] + ' ' + split_line[2] + ' ' + 
                                 split_line[4] + ' ' + split_line[5] + ' ' + split_line[6] + ' ' + 
                                 split_line[11] + ' ' + split_line[13])
                    mx_list.append(todays_mx)           
                else:
                    todays_mx = (split_line[0] + ' ' + split_line[1] + ' ' + split_line[2] + ' ' + 
                                 split_line[4] + ' ' + split_line[5] + ' ' + split_line[6] + ' ' + 
                                 split_line[11] + ' ' + 'NaN')
                    mx_list.append(todays_mx)
    # Again, order in terms of date
    sorted_mx_list = sorted(mx_list, key=lambda x: datetime.datetime.strptime(x.strip()[0:10], '%Y %m %d'))
    # Again, save list!
    outfile = open('C:/Users/laura/Documents/SMART/flare_list/mx_list.txt', 'w')
    for item in sorted_mx_list:
        outfile.write("%s\n" % item)

if __name__ == '__main__':
    main()