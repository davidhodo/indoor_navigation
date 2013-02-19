# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 00:00:10 2012

@author: hododav
"""

datafile='/Users/hododav/Documents/Data/2012-05-02/atrv_log_2012-05-02-19-52-09_filt_front_encoder.csv'

import csv
with open( datafile, "rb" ) as theFile:
    reader = csv.DictReader( theFile )
    for line in reader:
        print(line)
        # line is { 'workers': 'w0', 'constant': 7.334, 'age': -1.406, ... }
        # e.g. print( line[ 'workers' ] ) yields 'w0'