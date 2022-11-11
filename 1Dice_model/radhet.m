function [T]=radhet(DT,STV,HM,ACOFF,FHET,EXPO,EXHM,IPEN)
% function [T]=radhet(DT,STV,HM,ACOFF,FHET,EXPO,EXHM,IPEN);
%
% Calculate heating below mixed layer due to shortwave radiation


%function [temp,temp2]=radhet(DT,STV,HM,ACOFF,FHET,EXPO,EXHM,IPEN);

global DZ HMAX IMAX
global BS
global RAI RAS RAW CPA CPW CPI RLN VFS RLI RKI RKW RKS
global S T SS TS

%      CALL RADHET(DT,STV,HM,ACOFF,FHET,EXPO,EXHM,IPEN)
%      SUBROUTINE RADHET(DT,STV,HM,ACOFF,FHET,EXPO,EXHM,IPER)
%      REAL*8 S(0:500),T(0:500),RA(0:500),SS,TS,FHET(16,0:80)
%      REAL*8 BS(0:25),EXPO(0:500),HM,EXHM,D1,D2,HMAX,DZ
%      REAL*8 RAI,RAS,RAW,CPA,CPW,CPI,RLN,VFS,RLI,RKI,RKW,RKS
%      REAL*8 ACOFF,CON1,STV(4,0:241),DT
%      INTEGER I,IPER,IM,IMAX
%      COMMON/A/ S,T,RA,SS,TS
%      COMMON/B/ HMAX,DZ,IMAX
%      COMMON/C/ BS
%      COMMON/ICE/ RAI,RAS,RAW,CPA,CPW,CPI,RLN,VFS,RLI,RKI,RKW,RKS

IM=floor(HM/DZ);      % bottom level of mixed layer (original: IM=INT((HM-0.5*DZ)/DZ))
D1=HM-(IM+0.5)*DZ;   % vertical distance to mid-level 
D2=DZ-D1;            % vertical distance to next mid-level
TS=TS+DT*STV(1,1)*FHET(3,1)*(EXHM-exp(-ACOFF*(HM+D2)))/(D2*RAW*CPW);
T(IM+1)=(T(1)*D1+TS*(DZ-D1))/DZ;
CON1=DT*STV(1,1)*FHET(3,1)*BS(23);

Ttemp=T((IM+2):IPEN)+CON1*EXPO((IM+2):IPEN);
T((IM+2):IPEN)=Ttemp;

