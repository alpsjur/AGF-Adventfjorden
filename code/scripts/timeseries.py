#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 18:01:11 2022

@author: miriamsterl
"""

import numpy as np
import pandas as pd
import gsw
from datetime import datetime

df = pd.read_csv('../../data/ISA_CTD/ctd_isa_all_processed.csv',index_col=0)





df2 = pd.read_excel('../../data/ISA_CTD/ctd_isa_all.xlsx',skiprows=1)
df2.rename(columns={'Pressure (db)': 'Pressure', 'Salinity (psu)': 'Salinity', 'Temperature (degC)': 'Temperature'},
         inplace=True)

# Add full date as datetime
df2['Date'] = df2.apply(lambda row: datetime(int(row.Year),int(row.Month),int(row.Day),
                                           int(row.Hour),int(row.Minute),int(row.Second)),axis=1)
# Add absolute salinity
df2['SA'] = df2.apply(lambda row: gsw.SA_from_SP(row.Salinity,row.Pressure,row.Longitude,row.Latitude),axis=1)

# Add conservative temperature
df2['CT'] = df2.apply(lambda row: gsw.CT_from_t(row.SA,row.Temperature,row.Pressure),axis=1)

df2 = df2.drop(['Cruise', 'Station', 'Type', 'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second',
       'Bot. Depth (m)', 'Conductivity (S/m)',
       'Ship', 'Data Owner', 'File name', 'Citation', 'Instrument',],axis=1)