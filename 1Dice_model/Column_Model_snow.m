%********************************************************************
%
%     1DICE MODEL
%
%     NON ADVECTIVE AIR-SEA-ICE COLUMN MODEL
%
%     1DICE_SNOW  - SPECIFIED ATMOSPHERE (This version)
%     1DICE_AIR   - SIMPLE ATMOSPHERE
%     1DICE_CLOUD - NCAR ATMOSPHERE
%
%     ICE THICKNESS DISTRIBUTION DUE TO EXPORT AND RIDGING
%
%     Bjork G. (1989) J.Phys.Oceanography 19(1):52-67
%     Bjork G. (1997) J. Geoph. Res. 102 (C8) 18681-19698
%
%     1DICE_AIR : Bjork and Soderkvist (2002)
%                 J. Geoph. Res. 107 (C10) 10.1029/2000JC000723
%     1DICE_CLOUD : Soderkvist and Bjork (2004)
%                   Clim. Dyn. 22: 57-68 10.1007/s00382-003-0363-z
%
%      Advection of Heat added: by Lars H. Smedsrud, BCCR, 2008
%      Smedsrud et al, GRL (2008)  35, L20503
%      Smedsrud et al, Ocean Science (2010), 6, 10.5194/os-6-219-2010
%
%      Extension of depth and active Atlantic Layer:
%      Aleksi Nummelin, GFI, UiB, spring  2013
%      Nummelin et al, JGR (2015)119, doi:10.1002/2014JC010571
%
%      Conversion to Matlab started by Camille Li, GFI, UiB, spring  2013
%
%      Search for "F2Mconv" for notes and comments regarding
%      Fortran-to-Matlab conversion
%
%    ********************************************************************
%
%     *** STV *** (STATE VARIABLE) *****************
%     STV(N,1) ice concentration     STV(N,2) ice thickness
%     STV(N,3) snow thickness        STV(N,4) ice temperature
%
%     Number of ice classes, N ( N = 1 to NMAX )
%
%     *** FTAB (FORCING TABLE - MOTHLY MEANS - FROM FORCINGFILE) ***************
%       1: Short Wave Radiation 2: Long Wave Rad.  3: Ta - Atmosph. Temp   4: Humidity r
%       5: albedo of snow       6: snowfall (mm/d) 7: velocity_wind        8: variance_wind
%
%
%    F2Mconv: is 14 really ice/wind ratio? Because it is later set to something
%    called AICE in what used to be the forcing subroutine (now part of main
%    loop).
%
%     *** FHET ****(12 FLUXES OF HEAT FOR ALL ICE CLASSES )*********************
%
%     FHET(1,N) Short Wave IN                FHET(2,N) Short Wave NET IN
%     FHET(3,N) Short Wave Penetration(F0)   FHET(4,N) Long Wave IN
%     FHET(5,N) Long Wave OUT                FHET(6,N) Sensible heat
%     FHET(7,N) Latent heat                  FHET(8,N) Net Upward Flux 'Air' from surface (Fa)
%     FHET(9,N) Ice Upward Flux 'Snow' (Fs)  FHET(10,N) Flux to ice 'Bottom'(Fb)
%     FHET(11,N) Flux from Water(Fw)         FHET(12,N) Total flux to ice and snow (F0+Fa)
%
% ***************************************************************************
%      "Illustration of the flows" - Johanna Nilsson, Gothenburg,
%
%                            Ts      /\
%    _____                ____.______||Fa___________   _
%         |              |     \                    |  |          |
%     snow|              |      \    /\      snow   |  h          |
%    _____|  ||Fa        |_||____\___||Fs___________| _|_ _  _    |
%         |~~\/~~~~~~~~~~| \/(F0) \_                |     |  |    V z-direction
%         |  ||          |          \_ T            |     |  |
%      ice|  \/          |            \.       ice  |     H  |
%         |  (1-aice)Fsw |    /\Fb      \   Tf      |     |  |
%    _____|              |____||_________\._________|    _|_ |
%                             /\          |         |        |
%                             ||Fw        |         |        Hm
%                                         |         |        |
%                                         .Tm       |        |
%    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _|       _|_
%            ||                                     |
%            \/(1-aice)Fsw*Iw*exp(-ae*Hm)           |
%                                                   |
% ****************************************************************************

clear all % Clears all previous variables
%close(figure(1), figure(2), figure(3))       % Clears existing figures

global S T SS TS RA
global DZ HMAX IMAX
global DT
global AS BS AREA
global RAI RAS RAW CPA CPW CPI RLN VFS RLI RKI RKW RKS

% *** Read parameters that used to be in parameter file ***
% Each block is one row from the original parameter input file
% 1DICE_air_param.d.
% Original parameter file variable names are indicated if different.
% CL 02/2013

RMY=1;      % NRY number of years
NDT=360;      % timesteps per year (NDT = 8640 for 1h resolution)
%NDT=720;      % timesteps per year (NDT = 8640 for 1h resolution)

HMAX=300.;    % vertical dimension of model domain (m)
DZ=2.0;       % vertical resolution (m)
Depth=(0:DZ:HMAX); % Used for plotting only

% Specify Sea Ice Area EXport
% Arctic Basin mean annual export is 883.000 km�, i.e. 28.0e3 m�/s
% AEX=28.0e3;      % Ice Export [m�/s] , removed EXCON (Flag to switch off Export)
AEX=0.0;         % Rember that ice export removes freshwater ...

OINI(1,1)=32.0;  % SI1 Arctic salinity
OINI(1,2)=34.8;  % SI2 Atlantic salinity in Fram Strait
OINI(2,1)=-1.7;  % TI1 Arctic temperature
OINI(2,2)=0.5;   % TI2 Atlantic temperature in Fram Strait

NMAX=8;          % Maximum number of ice classes
HTH=0.25;        % some ice thickness limit (see checkmelt.m, for example) \CL
SICE=5.;         % ice salinity
GFAC=0.15;       % Value from parameter file, Factor in Shelf/Deep distribution NB: Not found in 1DICE_param.d

% Moved parameters for calculation of solar forcing to main routine - Lars H. 27. April 2016

% Also removed ICI(1) => ICI(8). These do not make sense in Matlab, because
% it is no longer needed to compile the code, so one can switch off fluxes
% and export directly in the code instead - Lars H. 27. April 2016

  SECY=3.1104e7;            % seconds in a year (360 days per year)
  DT=SECY/NDT;              % seconds per timestep
  TMAX=RMY*SECY;            % total simulation time (s): RMY is numbers of years
  IMAX=floor(HMAX/DZ)+1;    % total number of vertical level interfaces in model domain (depth/dz+1)
  NTMAX=RMY*NDT;            % total number of timesteps (NDT is timesteps per year)

  ENR=0.001;                % convergens limit newton-r
  AICE=0.01;                % wind/ice speed ratio
  AREA=8.0E+12;             % Area of model domain 8 mill km**2 is Arctic Basin
  DKF=7.0E-6;               % Diffusion coefficient, LHS, sets mixing below mixed layer
                            % 'Standard' value was 5.0E-6

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initial time step set below, should NOT be dependant on DT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  day1=30*11+1;           % initialization Julian day (August 15) - 30-day months /CL
  TIME=(day1-1)*SECY/NDT+SECY/(NDT*2);  % halfway through first timestep on day1 (s) /CL
                          % TIME starts at 0s on January 1 /CL

  SECM=SECY/12;         % seconds per month
  SECD=24*3600;         % seconds per day
  NDLY=NTMAX-(NDT-1);   % first time step (day) of the last year

  disp(['Running 1DICE_snow ......',num2str(RMY),' years, ',num2str(NTMAX),' timesteps']);
  disp(['Vertical ocean domain: ',num2str(HMAX),'m, ',num2str(DZ),'m levels, ',num2str(IMAX-1),'  levels']);

 % Fluxes depends on some constants and can be switched off HERE (set to zero)
 % Short Wave radiation is monthly varying in the forcing file
 SIG=5.67e-8;        % Outgoing Long Wave Radiation - Stefan-Bolzmann Coeff. 5.67e-8 is normal value
 CCS=1.75e-3;        % Sensible Heat Flux Coeff -  1.75e-3 is normal value
 CCE=1.75e-3;        % Latent Heat Coeff. -   1.75e-3 is normal value
 CCH=3.0e-4;         % Ocean Heat Flux Coeff. - 3.0e-4 is normal value

 % Solar forcing and penetration into ocean
  I0=0.17;         % Intensity of light at surface
  IW=0.5;          % Short Wave reflection
 % ACOFF=0.08;      % Extinction coeff - Arctic Ocean is very clear
  ACOFF=0.5;        % Extinction coeff - Muddy Svalbard fjords
  FRCUT=0.0001;    % S-W CUT
  EXPO=exp(-[1:IMAX]'*ACOFF*DZ);      % Beer's Law: Iz = I0*exp(-k*z)
  IPEN=-round(log(FRCUT)/(ACOFF*DZ)); % Shortwave penetration depth
  % IPEN=min(IPEN,HMAX-4); % Limits to bottom depth  % LHS 1. Nov 2016
  EXIP=exp(-ACOFF*(IPEN+.5)*DZ);

% *** Model Constants, should not change these
  SATL=OINI(1,2);     % Atlantic salinity in Fram Strait (used in outflow calculation) < param
  TATL=OINI(2,2);     % Atlantic temperature in Fram Strait (used in outflow calculation) < param
  RAATL=density(SATL,TATL);
  TREF=0;
  EPS=0.9;            %ICEDENSITY/WATERDENSITY
  RF=0.05;            %RICHARDSON FLUX NUMBER
  CDICE=5.5e-3;       %ICE/WATER DRAG COEFF
  CDINT=1.e-4;        %INTERNAL DRAG COEFF
  RM0=RF/sqrt(CDICE);
  G=9.81;             %GRAVITATIONAL ACCELERATION
  OMEG=7.29e-5;       %ANGULAR VELOCITY
  FIG=80;             % mean latitude
  FIR=pi*FIG/180;
  F=2.*OMEG*sin(FIR); %CORIOLIS PARAMETER
  RHC=0.2;            %CONSTANT FOR EKMAN-DEPTH
  EEI=0.99;           %LONG-W EMISSIVITY ICE
  EES=0.97;           %LONG-W EMISSIVITY SNOW

  RAA=1.3;            %ATM. DENSITY (kg/m^3)
  RAS=330;            %SNOW DENSITY (kg/m^3)
  RAI=900;            %ICE DENSITY (kg/m^3)
  RAW=1000;           %WATER DENSITY (kg/m^3)
  CPA=1.0044e+3;      %SPECIFIC HEAT OF AIR
  CPI=2.11e+3;        %SPECIFIC HEAT OF PURE ICE
  CPW=4.18e+3;        %SPECIFIC HEAT OF WATER

  RP0=1013;           %SURFACE PRESSURE (hPa)
  RLN=2.49E+6;        %HEAT OF VAPORIZATION
  RLI=334.8E+3;       %HEAT OF FUSION OF ICE
  VFS=RAS*RLI;        %VOLUMETRIC HEAT OF FUSION SNOW
  RKI=2.034;          %THERMAL CONDUCTIVITY ICE
  RKS=0.3097;         %  "          "       SNOW
  RKW=0.60;           %  "          "       WATER
  TFFW=273;           %FRESH WATER FREEZING TEMP (K)
  BCON=0.0543;        %TF=-BCON*S
  BETA=0.13;          %CONST. IN CONNEC. WITH HEAT CONDUCTION

  % Some useful FORTAN constants - these should be "removed" later
  ALF(1)=0.44;        %
  ALF(2)=0.28;        %ALFA=MIN(ALF(1)*HICE**ALF(2)+ALF(3),ALF(4))
  ALF(3)=0.08;        %
  ALF(4)=0.64;        %

  AS(1)=2.7798202e-6; % related to saturation pressure
  AS(2)=-2.6913393e-3;
  AS(3)=0.97920849 ;
  AS(4)=-158.63779;
  AS(5)=9653.1925;
  AS(6)=4.0*AS(1);
  AS(7)=3.0*AS(2);
  AS(8)=2.0*AS(3);

  BS(1)=EES*SIG; % Emmisivity of snow * Stefan Boltz
  BS(2)=EEI*SIG; % Emmisivity of ice * Stefan Boltz
  BS(3)=RAA*CPA*CCS;
  BS(4)=0.622*RAA*RLN*CCE/RP0;
  BS(5)=CPW*RAW*CCH;
  BS(6)=4*BS(1);
  BS(7)=BS(2)+BS(3)*AS(4);
  BS(8)=BS(3)*AS(7);
  BS(9)=BS(3)*AS(8);
  BS(10)=BETA*SICE;
  BS(11)=CPI*RAI;
  BS(12)=CPW*RAW;
  BS(13)=RAI*RLI;
  BS(14)=-BCON*SICE*RLI/CPI;
  BS(15)=BCON*SICE;
  BS(16)=RAS/RAW;
  BS(17)=TFFW^4;
  BS(18)=TFFW*(TFFW*(TFFW*(AS(1)*TFFW+AS(2))+AS(3))+AS(4));
  BS(19)=RAI/RAW;
  BS(20)=G*DZ/(RAATL*F); % F = Coriolis
  BS(21)=sqrt(CDICE);
  BS(23)=exp(ACOFF*DZ/2.)*(1.-exp(-ACOFF*DZ))/(DZ*RAW*CPW);

% Read Forcing File
  input_force;
  % Multiply by scaling factors in column forcing(:,13)
  FTAB=forcing(1:8,1:12).*repmat(forcing(1:8,13),[1 12]);
  % Calculate Wind variance from wind + STD (W^3+3W*std^2)^(1/3)
  FTAB(10,:)=(FTAB(7,:).^3 + 3*FTAB(7,:).*FTAB(8,:).^2).^(1/3);

% Initial Conditions
  disp('Reading initialization files ... ');

  % *** read ocean profile ....
  input_ocean='initial_ocean';
  eval(input_ocean);
  T=TMP1(:,2);
  S=TMP1(:,3);
  %T=ones(HMAX/DZ+1,1); T(1:5)=T(1:5)*+4.0; T(6:end)=T(6:end)*2.0;
  %S=ones(HMAX/DZ+1,1); S(1:5)=S(1:5)*32.5; S(6:end)=S(6:end)*34.5;
  %HM=11;	       % mixed layer depth
  %TS=3.0;  % split layer temperature
  %SS=33.0; % split layer salinity

  disp(['10 dz ocean temperature: ',num2str(T(10))]);
  disp(['10 dz ocean salinity: ',num2str(S(10))]);

  % *** read ice profile ....
  input_ice=  'initial_ice';
  eval(input_ice);
  % sum over ice classes
  HICE=STV(:,1)'*STV(:,2); % concentration * ice thickness
  HSNO=STV(:,1)'*STV(:,3); % concentration * snow thickness

  disp(['Initial T and S at the surface: ',num2str(T(1)),'deg C, and  ',num2str(S(1)),' psu']);
  disp(['and at the bottom: ',num2str(T(IMAX)),' deg C, and  ',num2str(S(IMAX)),' psu']);
  disp(['    ']);
  disp([' Initial mean ice thickness:  ', num2str(HICE),'m, and snow thickness:  ',num2str(HSNO),' m']);
  disp(['    ']);
  disp(['  Integrating over time .....  ']);
  disp(['    ']);

% Matrix for storing variables
HICE_timeseries = zeros(NTMAX,1);
HSNO_timeseries = zeros(NTMAX,1);
T_timeseries = zeros(NTMAX,length(T));
S_timeseries = zeros(NTMAX,length(S));

HICE_timeseries(1) = HICE;
HSNO_timeseries(1) = HSNO;
T_timeseries(1,:) = T;
S_timeseries(1,:) = S;

%--------------------------------------------
% *** MAIN LOOP ***
%--------------------------------------------
  for itim=1:NTMAX % time-stepping (originally itim = IC)

  % calculate month, day of simulation (itim = timestep, NDT= timesteps/year)
    iyr =floor((itim-1)/NDT); % run starts at iyr=0
    iday=itim-iyr*NDT;
    imon=floor(((itim-1)/NDT-iyr)*12)+1;

    TIME_D=TIME-floor(TIME/SECY)*SECY;      % Julian day of time step, regardless of year (s)
    imonth=floor((TIME_D+SECM/2)/SECM);     % index of "previous" month
    TDEL=(TIME_D-(imonth-0.5)*SECM)/SECM;   % fraction of month ahead of previous mid-month forcing value

    % Reads the correct forcing for this day
    vnames={'FR','FL','TA','RH','ASN','SRATE','WIND','STD','WIR','MWIN'};
    for iv=1:length(vnames)
        if imonth>=1 && imonth<12       % For January => November
          vtemp=FTAB(iv,imonth) + TDEL*diff(FTAB(iv,imonth:(imonth+1)));
        else % last half of December and first half of January
          vtemp=FTAB(iv,12)     + TDEL*(FTAB(iv,1)-FTAB(iv,12));
        end
        assignin('base',char(vnames(iv)),vtemp);
        clear vtemp;
    end

    clear STD WIR % Clear Wind standard dev (STD) and Wind/Ice speed ratio (WIR) - not used

  % Add snow on top of all ice classes
    for ins=1:N
       STV(ins,3)=STV(ins,3)+(SRATE*DT);      %  Let it snow ..
    end

  % Calculate friction velocity
    [UST,UICE,UFRA]=ustar(RAW,RAA,STV,AICE,WIND,MWIN);

  % Calculate heat fluxes FHET(I,CLASS)
    %     1: SW IN         2: SW NET IN   3: SW PEN.     4: LW IN
    %     5: LW OUT        6: SENSIBEL    7: LATENTN     8: NET. UPP. SURF
    %     9: ICE UP SURF. 10: LO SUR ICE 11: WAT TO SUR 12: TOT ATM
    %     HEAT FLOW AT THE SNOW/ICE SURF. EQUALS FHET(9,I)
    [FHET,TSF,OPGR,FWAT] = heatflux(N,ENR,HTH,TFFW,BCON,ALF,STV,I0,IW,...
                               FR,FL,TA,WIND,UICE,RH,ASN);

    % Calculate lower surface growth/melt and upper surface snow or ice melt;
    % return new ice temperature and mixed layer salinity/temperature
    [STV,HICE,HSNO,HM,WIME,WSME,DHS,DHU,DHL]=icegrowth(N,HM,SICE,OPGR,ACOFF,EXPO,IPEN,EXIP,...
                                                       GFAC,STV,TSF,FHET,HICE,HSNO);

    % disp([' Time step Delta T: (',num2str(itim),'): Mean ice thickness = ',num2str(HICE)])
    % disp([' Time step Delta T: (',num2str(itim),'): sum(STV(:,1)) ice conc = ',num2str([sum(STV(:,1))])])
    if sum(STV(:,1))>1.1; error('total ice concentration > 1'); end

    % Sort ice classes by thickness
    [STV]=sortice(STV);

    % Merge open water classes if necessary
    [N,STV]=compress(N,STV);

    % Calculate surface buoyancy flux
    [BF]=buoyflux(G,WIME,SICE,FWAT);

    %   CALL MIXING(DT,IC,G,F,RHC,RM0,UST,BF,HSNO,HICE,HM,HMON,HEK,PAR)
    [HM]=mixing(DT,G,F,RHC,RM0,UST,BF,HSNO,HICE,HM);

    % CALL DIFFUSION, mixing below mixed layer - added Oct 14, 2022 LHS
    [T,S]=diffusion(DT,DKF,HM);

    %   Advection could be added at a later stage like done for the Barents Sea:
    %         CALL ADVECTHEAT(DT,AREA,IMAX,T,ADVT,VOLFLUX)  !    LHS
    %         CALL ADVECTSALT(DT,AREA,IMAX,S,ADVS,VOLFLUX)  !    LHS

    % Export part of each ice class (can increase N)
    if AEX>0          %  AEX is constant export specified as forcing
      [N,STV,HSNO,HICE]=iceexp(AEX,N,STV);
    end

    if isempty(find(STV(:,1)<0, 1))~=1; error('C Negative ice concentrations found'); end

    % Sort and, if necessary, merge two most similar ice classes of other ice
    while N>NMAX
      [N,STV]=mergeice(N,STV);
    end

    % Store
    HICE_timeseries(itim) = HICE;
    HSNO_timeseries(itim) = HSNO;
    T_timeseries(itim,:) = T;
    S_timeseries(itim,:) = S;

    %disp(num2str(N));


  % Plotting every 10 days
      if rem(itim,10)==0

      % Figure: Ocean column with T and S
      figure(1)
      subplot(1,2,2), hold on,
      plot(S(:,1),Depth,'g--')
      subplot(1,2,1), grid on, hold on,h=gca;
      plot(T(:,1),Depth,'r--')

      % Figure: Ice and Snow Thickness
      figure(2), hold on,
      h1=plot(itim/NDT,HICE,'ko'); % Mean ice thickness
      set(h1,'LineWidth',2)
      h2=plot(itim/NDT,HSNO+HICE,'b*'); % Mean snow thickness on top
      plot(itim/NDT,DHU+HICE,'c.')  % Change in ice thickness due to Upper boundary
      %plot(itim,DHL+HICE,'m.') % Change in ice thickness due to Lower boundary

      % Figure: Surface fluxes over time
    %     1: SW IN         2: SW NET IN   3: SW PEN.     4: LW IN
    %     5: LW OUT        6: SENSIBEL    7: LATENTN     8: NET. UPP. SURF
    %     9: ICE UP SURF. 10: LO SUR ICE 11: WAT TO SUR 12: TOT ATM
    %     HEAT FLOW AT THE SNOW/ICE SURF. EQUALS FHET(9,I)

      figure(3), hold on,
      plot(itim/NDT,mean(FHET(1,:)),'r*');                 % Solar in - positive
      plot(itim/NDT,mean(FHET(4,:))+mean(FHET(5,:)),'ko'); % Net longwave - negative
      plot(itim/NDT,mean(FHET(6,:)),'b<');                 % Sensible - negative
      plot(itim/NDT,mean(FHET(7,:)),'cs');                 % Latent - negative

      % disp('    ');          % Creates blank line on screen
      % disp(['*** month ',num2str(itim/30),' completed']);
      % disp(['*** time step ',num2str(itim)]);
      % disp(['*** ice classes N = ',num2str(N)]);
      end

   TIME=TIME+DT; % Next time step (1 day in seconds)
   end  % for itim=1:NTMAX, % time-stepping (originally itim = IC)
%--------------------------------------------
% *** END OF MAIN LOOP ***
%--------------------------------------------

% Fresh up figures at the end
      figure(1)
      subplot(1,2,2), grid on, h=gca;
      set(h,'YDir','reverse');
      axis([ 28.0 35.0 0 HMAX])
      set(h,'Fontsize',16,'yTick',[0,50,100,150,200,250])
      set(h,'Fontsize',16,'xTick',[28.0:1.0:35.0])
      h=xlabel('Salinity');
      set(h,'Fontsize',16,'FontWeight','bold')
      % Replot last profile
      h=plot(S(:,1),Depth,'g-');
      set(h,'LineWidth',3)

      subplot(1,2,1), grid on, h=gca;
      set(h,'YDir','reverse')
      axis([ -4.5 4.5 0 HMAX])
      set(h,'Fontsize',16,'yTick',[0,50,100,150,200,250])
      set(h,'Fontsize',16,'xTick',[-4.0:2.0:4.0])
      h=xlabel('Temperature');
      set(h,'Fontsize',16,'FontWeight','bold')
      % Replot last profile
      h=plot(T(:,1),Depth,'r-');
      set(h,'LineWidth',3)

      figure(2), grid on,
      h=ylabel('Ice Thickness (m)');
      set(h,'Fontsize',16,'FontWeight','bold')
      h=xlabel('Years');
      set(h,'Fontsize',16,'FontWeight','bold')
      h=gca;
      set(h,'Fontsize',16,'yTick',[0,0.5,1,1.5,2.0,2.5,3.0,3.5])
      set(h,'Fontsize',16,'xTick',[0.2:0.2:max(itim)/NDT])
      axis([0 max(itim)/NDT 0 2.5])  % Normal ...

      figure(3), grid on,
      h=ylabel('Surface Fluxes [W/m�]');
      set(h,'Fontsize',16,'FontWeight','bold')
      h=xlabel('Years');
      set(h,'Fontsize',16,'FontWeight','bold')
      h=gca;
      set(h,'Fontsize',16,'yTick',[-200,-100,0,100,200,300])
      set(h,'Fontsize',16,'xTick',[0.2:0.2:max(itim)/NDT])
      axis([0 max(itim)/NDT -100 300])  % Normal ...

disp('    ');          % Creates blank line on screen
disp(['Finished Running 1DICE_air ......',num2str(RMY),' years, ',num2str(NTMAX),' timesteps']);
disp(['    ']);
disp([' Final T and S at the surface: ',num2str(T(1)),'deg C, and  ',num2str(S(1)),' psu']);
disp([' and at the bottom: ',num2str(T(IMAX)),' deg C, and  ',num2str(S(IMAX)),' psu']);
disp(['    ']);
disp([' Final mean ice thickness:  ', num2str(HICE),'m, and snow thickness:  ',num2str(HSNO),' m']);

% Save to file
%writematrix(HICE_timeseries,'../data/modeloutput/HICE.txt')


dlmwrite('../data/modeloutput/HICE.txt', HICE_timeseries)
dlmwrite('../data/modeloutput/HSNO.txt', HSNO_timeseries)
dlmwrite('../data/modeloutput/T.txt', T_timeseries)
dlmwrite('../data/modeloutput/S.txt', S_timeseries)




