#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 15:01:52 2018

@author: LACORBET
"""



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

'''
This code imports NOAA daily sunspot numbers and daily SMART detections and adds them to lists. 
The data is then smoothed using the convolve function and plotted against each other.
'''


def main(file_location = '/home/LACORBET_PROJECT2018/smart_harddrive/smart/data/NOAA/Daily_Sunspot_Number.txt',
         input_folder = '/home/LACORBET_PROJECT2018/smart_harddrive/smart/data/daily_json_files/'):
    
    #NOAA Data ---------------------------------------
    #-------------------------------------------------
      
    NOAA_dates=[] #Empty list for my NOAA_dates
    daily_sunspots=[] #Empty list for daily sunspots
    
    #Read in my file
    f=open(file_location, "r")
    #Split my file into lines
    lines=f.readlines()
    
    lines=[x.strip() for x in lines]
    
    #Strip lines for date information
    for el in lines:
        year=int(el[0:4])
        month=int(el[5:7])
        day=int(el[8:10])
        
        time=datetime.datetime(year, month, day)
        NOAA_dates.append(time) #Append list with dates
        
        
        sunspots=int(el[22:24])
        daily_sunspots.append(sunspots)


    
    kernel=np.ones(50)/50

    convolve=np.convolve(daily_sunspots,kernel,'same')
    
    
    
    #SMART Data--------------------------------------------
    #------------------------------------------------------
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
    
    SMART_sunspots=[]
    SMART_dates=[]
    
    for date_string in date_strings: #getting active regions for each day
        
         json_filename =input_folder + date_string[0] + "_properties.json" #taking in my json files
         json_data = json.load(open(json_filename)) #loading/opening the data
         number_json = list(json_data['posprop']['arid'].keys()) #taking the key values and making them into lists
         SMART_sunspots.append(len(number_json))
         SMART_dates.append(date_string[1])
         
    kernal2=np.ones(50)/50
    convolve2=np.convolve(SMART_sunspots,kernal2,'same')
    

#   #plot figures with different axes
    fig, ax1 = plt.subplots()
#     
    color='tab:blue'
    ax1.set_xlabel('Years')
    ax1.set_ylabel('Daily Sunspot Number', color=color)
    ax1.plot(SMART_dates, convolve, color=color)
    ax1.tick_params(axis='y', labelcolor=color)
     
    ax2=ax1.twinx()
     
    color='tab:red'
    ax2.set_ylabel('SMART Daily Detections', color=color)
    ax2.plot(SMART_dates, convolve2, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    
    plt.show()

    
def datetime_from_file_string(a):
    
    year = int(a[:4])
    month = int(a[4:6])
    day = int(a[6:8])
    hour = int(a[9:11])
    time1 = datetime.datetime(year, month, day, hour)
    return time1
        
        
        