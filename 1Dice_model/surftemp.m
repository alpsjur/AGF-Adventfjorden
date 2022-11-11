function [FLO,FSE,FLA,FCO,T0GN] = surftemp(ENR,TFFW,RJA,EEE,TIG,I0,TAG,CTA,WIND,FRN,FL,T0G)

% Calculates surface temperature and surface fluxes through iteration 
% FLO = FHET(5,I) i.e. Longwave OUT
% FSE = FHET(6,I) i.e. Sensible heat flux
% FLA = FHET(7,I) i.e. Latent heat flux
% FCO = FHET(9,I) i.e. Conductive Heat from top of snow/ice 
% T0GN = TSURF  i.e. new calculated surface temperature 
% T0G = existing value before iteration
% 
% Fluxes abd values calculated earlier: 
% RJA = RJAM conductive heat through ice
% FL = FHET(4,I) i.e. Long Wave IN
% FRN = FHET(2,I) i.e. Short Wave IN
% TIG = STV(I,4) i.e. ice temperature for this ice class
% I0  is light intensity at surface (constant parameter) 
% TAG = TA   i.e. Atmospheric temperature (forcing)
% CTA is black body radiation values
% WIND is wind value 
 
%  SUBROUTINE SURFTEMP(ENR,TFFW,RJA,EEE,TIG,I0,TAG,AS,CTA,
%      >           WIND,FRN,FL,FLO,FSE,FLA,FCO,T0G)
%       REAL*8 RJA,TI,I0,FRN,FRT,FL,TA,CTA,T0,FSE,FLA,T1,FT0,FPT0
%       REAL*8 T02,T03,T04,FLO,AS(8),BS(0:25),FCO,ENR,EEE,TFFW
%       REAL*8 TIG,TAG,T0G,WIND,C1,C2,C3,C4,D2,D3
%       COMMON/C/ BS

    global AS BS

      FRT=(1.-I0)*FRN;                   % (1-I0)*Sw_net_in
      D2=BS(3)*WIND;                     % rho_air*Cpa*Ccs*W  Sens
      D3=BS(4)*WIND;                     % 0.622*rho_a*Rln*Cce*W/RPO Lat
      C1=4*(EEE+D3*AS(1));               % 4(sig*emiss_s +D3*2.72*10^-6)
      C2=BS(8)*WIND;
      C3=BS(9)*WIND;
      C4=BS(7)*WIND+RJA;
      TI=TIG+TFFW;                      % Internal Ice temp or water T(1) [K] 
      TA=TAG+TFFW;                      % Atm temp [K]
      T0=T0G+TFFW;                      % Surface temp "guess" [K]
      T02=T0*T0;
      T03=T02*T0;                       
      T04=T03*T0;   
      
      % Initial calculation of T1 
      FLO=-EEE*T04;                     % Long Wave (sig*T^4)
      FSE=D2*(TA-T0);                    % Sensible
      FLA=D3*(CTA-(AS(1)*T04+AS(2)*T03+AS(3)*T02+AS(4)*T0)); % Latent
      FCO=RJA*(T0-TI);                   % RJA: combined heat conduction
      FT0=FRT+FL+FLO+FSE+FLA-FCO;        % Total Flux
      FPT0=-C1*T03-C2*T02-C3*T0-C4;
      T1=T0-FT0/FPT0; 
      
      while(abs(T1-T0) > ENR)        % convergens error newton-r  
          T0=T1;                       % Set to calcualted value 
          T02=T0*T0;
          T03=T02*T0;                       
          T04=T03*T0;   
          FLO=-EEE*T04;                     % Long Wave (sig*T^4)
          FSE=D2*(TA-T0);                    % Sensible
          FLA=D3*(CTA-(AS(1)*T04+AS(2)*T03+AS(3)*T02+AS(4)*T0)); % Latent
          FCO=RJA*(T0-TI);                   % RJA: combined heat conduction
          FT0=FRT+FL+FLO+FSE+FLA-FCO;        % Total Flux
          FPT0=-C1*T0^3-C2*T0^2-C3*T0-C4;
          T1=T0-FT0/FPT0; 
      end
      
      if(T0 > TFFW)    % T_surf > 0 deg C
        T0=TFFW;
      end
      
      T0GN=T0-TFFW;      % Calcualte T0G in Kelvin
end

