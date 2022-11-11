function [TIN] = icetemp(TI,DHET,HI)
% function [TIN] = icetemp(TI,DHET,HI);
%
% Calculate internal ice temperature with brine correction
%      TI = ice temperature (K)
%      DHET = heat flux (W/m2) 
%      HI = ice thickness (m)

global BS
global CPI RLI

%      REAL*8 RAI,RAS,RAW,CPA,CPW,CPI,RLN,VFS,RLI,RKI,RKW,RKS
%      REAL*8 BS(0:25),TMEL,DF,A,B,TI,DHE,HI,TIN,EPS
%      COMMON/C/ BS
%      COMMON/ICE/ RAI,RAS,RAW,CPA,CPW,CPI,RLN,VFS,RLI,RKI,RKW,RKS
%      PARAMETER (EPS=1.E-6)

%      BS(11)=CPI*RAI; specific heat of ice * ice density [should be (J/K/m^3)]
%      BS(15)=BCON*SICE; 
%      RLI; heat of fusion of ice (J/kg)

eps=1e-6;

if abs(DHET)<eps,
  TIN=TI;
else
  TMEL=-BS(15);
  DF=DHET/(BS(11)*HI); % K?
  A=TI+RLI*TMEL/(CPI*TI)+DF; % RLI/(CPI*TI) is ratio (heat of fusion):(heat in ice of temperature TI)
  B=RLI*TMEL/CPI;
  TIN=A/2.0-sqrt(A*A/4-B);
  if TIN> TMEL,
    disp(['     icetemp warning TIN>TMEL: TIN=',num2str(TIN),', TMEL=',num2str(TMEL),', HI=',num2str(HI)]);
    TIN=TMEL-eps;
  end
end
