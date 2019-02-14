#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 16:54:40 2018

@author: LACORBET
"""

'''
This code plots numerous SMART detected longterm properties in one figure
'''

import sunpy.map
import os
from matplotlib import pyplot as plt
import json
import datetime
import itertools
import matplotlib as mpl
from astropy.coordinates import SkyCoord
import astropy.units as u
from sunpy.visualization import wcsaxes_compat
import pandas
import matplotlib.dates as mdates
import numpy as np


#Define input folders
def main(input_folder='/home/LACORBET_PROJECT2018/smart_harddrive/smart/data/daily_json_files/', 
         smart_folder='/home/LACORBET_PROJECT2018/smart_harddrive/smart/data/maps_detections/',
         output_folder='/home/LACORBET_PROJECT2018/smart_harddrive/smart/data/longterm_graphs/'):

    #Take in the json files, sort and name them filenames_json
    filenames = sorted(os.listdir(input_folder))
    filenames_json = [x for x in filenames if ".json" in x]
     
    #Take in the filename dates
    filename_dates = [datetime_from_file_string(x) for x in filenames_json]
    start_date, end_date = filename_dates[0],  filename_dates[len(filename_dates)-1]

    date_strings = [] #Adding my dates to a list
    for index, value in enumerate(filename_dates):
        if start_date <= value <= end_date:
            date_strings.append([filenames_json[index][:13], value])

        
#Defining my empty lists
    psl_length=[] 
    tot_area=[]
    new_tot_area=[]
    bmax=[]
    new_bmax=[]
    new_psl_length=[]
    dates=[]
    for date_string in date_strings: #getting active regions for each day
        
         json_filename =input_folder + date_string[0] + "_properties.json" #taking in my json files
         json_data = json.load(open(json_filename)) #loading/opening the data
         number_json1 = list(json_data['pslprop']['psllength'].values())#taking the key values and making them into lists
         number_json2=list(json_data['magprop']['totarea'].values())
         number_json3=list(json_data['magprop']['bmax'].values())
         
         tot_area.append(number_json2)
         psl_length.append(number_json1)
         bmax.append(number_json3)
         dates.append(date_string[1])
         
    for i in tot_area:
        y=np.sum(i)
        new_tot_area.append(y)
    #dictionary=dict(zip(dates, number_ar))
    for i in psl_length:
       x=np.sum(i)
       new_psl_length.append(x)
       
    for i in bmax:
        z=np.sum(i)
        new_bmax.append(z)
       
    kernel=np.ones(50)/50
#    
#    
    fig=plt.figure(figsize=(10,12))
    
        
        
    ax1=fig.add_subplot(3, 1, 1)
    convolve2=np.convolve(new_tot_area,kernel,'same')
    ax1.plot(dates,convolve2)
    plt.ylabel('Total AR Area Per Day [M.S.H]')
    
    
    ax1=fig.add_subplot(3, 1, 2)
    convolve3=np.convolve(new_bmax,kernel,'same')

    ax1.plot(dates,convolve3)
    plt.ylabel('Total AR Max. Magnetic Field Per Day [G]')
    
    ax1= fig.add_subplot(3, 1, 3)
    convolve1=np.convolve(new_psl_length,kernel,'same')  
    
    plt.ylabel('Total AR PSL Length Per Day [Mm]')
    ax1.plot(dates,convolve1)

    plt.xlabel('Year')
    #plt.savefig(output_folder + 'multiple_properties')
    plt.show()
    plt.close()
        
#    plt.plot(dates,number_ar)
#    plt.xlabel('Year')
#    plt.ylabel('# Active Regions')
#    plt.show()
#    plt.savefig(output_folder + "longterm")
#    plt.close()
    
def datetime_from_file_string(a):
    """
    Extracting a datetime object from the SMART json file
    - convert timestring to datetime object
    """
    year = int(a[:4])
    month = int(a[4:6])
    day = int(a[6:8])
    hour = int(a[9:11])
    time1 = datetime.datetime(year, month, day, hour)
    return time1