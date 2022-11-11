"""
Computing the climatological atmospheric forcing
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



#%% Functions to make climatologies

def getMonths(df):
    """
    Extracts months for each date in dataframe
    """
    months = np.empty(12, dtype=object)
    for j in range(months.shape[0]):
        months[j] = [] # initialize as empty lists
    for i in range(len(df)):
        m = int(df['t'][i][3:5])
        months[m-1].append(i)
    return months


def makeClim(df,variable):
    months = getMonths(df)
    clim = np.zeros(12)
    for i in range(12):
        clim[i] = np.nanmean([float(df[variable][j]) for j in months[i]],axis=0)
    return clim


#%% Copy data from model for snow albedo and ice/wind ratio

snow_albedo = np.array([0.85,0.84,0.83,0.81,0.82,0.78,0.64,0.69,0.84,0.85,0.85,0.85])
ice_wind_ratio = np.array([.01,.01,.01,.01,.01,.02,.02,.02,.01,.01,.01,.01])


#%% Radiation data from Janssonhaugen Vest

df_Janssonhaugen = pd.read_csv('../../data/atmospheric_forcing/Janssonhaugen.csv',
                               skipfooter=1,delimiter=';')
df_Janssonhaugen.rename(columns={'Time(norwegian mean time)': 't', 
                             'Mean global radiation (1 t)': 'SW',
                             'Mean longwave radiation (1 h)': 'LW'},
                        inplace=True)

SW = makeClim(df_Janssonhaugen,'SW')
LW = makeClim(df_Janssonhaugen,'LW')


#%% Air temp, humidity, precipitation and wind data from Svalbard Lufthavn

df_Lufthavn = pd.read_csv('../../data/atmospheric_forcing/Lufthavn.csv',
                               skipfooter=1,delimiter=';')
df_Lufthavn.rename(columns={'Time(norwegian mean time)': 't', 
                             'Air temperature': 'Ta',
                             'Relative air humidity': 'RH',
                             'Precipitation (1 h)': 'precip',
                             'Mean wind speed': 'wind'},
                   inplace=True)

# TODO:
    # - check what's up with this dataset. Seems to be missing data and weird dates??
    # - if not working, use other station, or shorter time range
    # - for RH: divide numbers by 100
    # - for precip: convert to snowfall in m/s (convert density!)
    # - for wind speed: check units + compute stdev



