"""
Computing the climatological atmospheric forcing
"""

import numpy as np
import pandas as pd



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


def makeClim_stdev(df,variable):
    months = getMonths(df)
    clim = np.zeros(12)
    for i in range(12):
        clim[i] = np.nanstd([float(df[variable][j]) for j in months[i]],axis=0)
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
df_Lufthavn = df_Lufthavn.replace(to_replace='-',value=np.nan)


Ta = makeClim(df_Lufthavn,'Ta')
RH = makeClim(df_Lufthavn,'RH')/100 # convert to number between 0 and 1
wind_speed = makeClim(df_Lufthavn,'wind')
wind_stdev = makeClim_stdev(df_Lufthavn,'wind')

# precipitation: convert from mm/hour to m/s and from water to snow density
snowfall = makeClim(df_Lufthavn,'precip')*(1e-3/3600)*(1000/330)
# for months where Ta > 0, set snowfall to 0
snowfall[np.where(Ta>0)] = 0


#%% Combine the climatologies into 1 DataFrame

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
data = {'Shortwave radiation': SW, 'Longwave radiation': LW, 
        'Air temperature': Ta, 'Relative humidity': RH, 
        'Snow albedo': snow_albedo,'Snowfall': snowfall, 
        'Wind speed': wind_speed, 'Wind standard deviation': wind_stdev,
        'Ice/wind ratio': ice_wind_ratio}
df_forcing = pd.DataFrame(data=data,index=months)
df_forcing.to_csv('../../data/atmospheric_forcing/forcing.csv')




