% Set initial Ice and snow thickness and temperature
% STV is 'STate Variable' with four columns
% One row for each ice class 'N'
% Column#1 = Ice Concentration
% Column#2 = Ice thickness
% Column#3 = Snow thickness
% Column#4 = Ice temperature

% Initial Open Water and no snow;
STV=[...
    1.0    0.00    0.00    -1.69 ;...
];

% Initial 50% Ice concentration and 50% open water
% 25% 0.5 m  ice thickness 0.1 m snow
% 25% 1.0 m ice thickness 0.2 m snow
% N = 3, and the sum of STV(1,:) must be 1.0

% STV=[...
%     0.5    0.00    0.00    -1.69 ;...
%     0.25   0.50    0.10    -1.69 ;...
%     0.25   1.00    0.20    -1.69 ;...
% ];

N=size(STV,1); % ice classes

