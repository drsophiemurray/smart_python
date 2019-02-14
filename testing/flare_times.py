# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 15:24:32 2018

@author: laura
"""
'''
This file imports M and X flare data and exports time and dates of those flares to a new list. 
'''


file_location='C:/Users/laura/Documents/SMART/flare_list/mx_list.txt'

f=open(file_location,"r")
lines=f.readlines()
result=[]
for x in lines:
    result.append(x.split(' ')[0:4])
    
print(result)

outfile = open('C:/Users/laura/Documents/SMART/flare_list/mx_time.txt', 'w')
for item in result:
    outfile.write("%s\n" % item)