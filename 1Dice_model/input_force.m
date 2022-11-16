% Forcing file with monthly mean forcing

% Based on the original forcing input file 1DICE_forc.d
% FTAB(:,1:12) monthly mean forcing, JANUARY => DECEMBER
% FTAB(:,13) scaling factor to multiply forcing by
% CL 02/2013, LHS 27/4-16
% Sign convention is POSITIVE FLUX WARMS OCEAN,  NEGATIVE FLUX COOLS OCEAN

forcing=[...
     0.0    5.1   32.90 142.41 256.81 302.03 232.64 132.89  47.63   9.61   0.0    0.0    1.000e0 ;... % 1 ShortWave rad [W/m2]- AARI
   162.60 154.99 158.49 179.53 236.88 278.20 279.11 277.79 259.57 213.67 179.82 164.78   1.000e0 ;... % 2 LongWave rad [W/m2]
   -29.4  -32.0  -31.2  -22.7  -10.5   -1.2   -0.6   -1.4   -9.9  -19.8  -27.8  -29.2    1.000e0 ;... % 3 Ta [deg C]
   0.85   0.85   0.88   0.88   0.89   0.93   0.93   0.94   0.88   0.85   0.85   0.85   1.000e0 ;... % 4 Humidity [ratio 0-1] RH
   0.85   0.84   0.83   0.81   0.82   0.78   0.64   0.69   0.84   0.85   0.85   0.85   1.000e0 ;... % 5 albedo_Snow [ratio 0-1]
   0.83   0.83   0.83   0.83   5.00   0.00   0.00   0.00  12.87  12.87   0.83   0.83   4.244e-9;... % 6 snowfall [mm/d => mm/s] (Maykut, 1982)
   5.6    5.7    5.3    5.1    5.0    5.2    5.2    5.4    6.2    6.2    5.8    5.5    1.000e0 ;... % 7 Wind_vel  [m/s] Maykut (1982)
   3.3    3.3    3.0    3.0    2.7    2.9    3.1    3.2    3.8    3.5    3.5    3.2    1.000e0 ;... % 8 Wind_std_dev
    .01    .01    .01    .01    .01    .02    .02    .02    .01    .01    .01    .01   1.000e0 ;... % 9. Ice/wind ratio (Not used?)
];


