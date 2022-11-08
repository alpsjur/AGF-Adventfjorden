function [HM]=npconv(HM,HSNO,HICE)
% function [HM]=npconv(HM,HSNO,HICE);
%
% Original comment: Calculate non-penetrative mixing and heat flux 
% Adjust mixed layer depth, salinity, temperature if not stable \CL
% 
% This function is only called when convectiopn into the next layer happens
% Gradual increase into D1 and D2 does not trigger this function

global S T SS TS
global DZ IMAX
global BS

%      SUBROUTINE NPCONV(HM,HSNO,HICE)
%      REAL*8 S(0:500),T(0:500),RA(0:500),SS,TS,HM
%      REAL*8 HMAX,DZ,HSNO,HICE,D1,D2,HW,BS(0:25),RAM,RS,SM,TM

IM=ceil((HM-0.5*DZ)/DZ);    % number of whole layers within HM  (original: IM=INT((HM-0.5*DZ)/DZ))
D1=HM-(IM-0.5)*DZ;          % upper delta Z within HM (original: D1=HM-(IM+0.5)*DZ)
D2=DZ-D1;                   % lower delta Z below HM to next layer
HW=HM-BS(19)*HICE-BS(16)*HSNO;
RAM=2;
RS=1;
SM=S(1);  % Present mixed layer Salt
TM=T(1);  % Present mixed layer temp

disp(['     npconv: new mixed layer convecting,  TM, SM, = ',num2str([TM SM])]);

while RAM>RS, % mixed layer deepens to where density of mixed layer = density of level
    
  if IM >= IMAX-1,     % Bottom fix LHS jan 2008
     SM=SM;
     TM=TM;
     HM=HM;
     IM=IM;
     SS=S(IMAX);
     TS=T(IMAX);
     HW=HW;
     RS=RAM+0.01;     % Stop mixing 
     disp ([ 'Convecting to bottom.'])
  else    
  SM=(HW*SM+D2*SS+D1*S(IM+2))/(HW+DZ);  % LHS Changed from S(IM+2)
  TM=(HW*TM+D2*TS+D1*T(IM+2))/(HW+DZ);  % LHS Changed from T(IM+2)
  HM=HM+DZ;  % Mixed layer grows one more DZ
  IM=IM+1;   % add one more layer in 
  SS=S(IM+1); % adjust salinity below ML 
  TS=T(IM+1); % adjust temperature below ML
  HW=HW+DZ;
  RAM=density(SM,TM);
  RS=density(SS,TS);
  end
end

S(1:IM)=SM;
T(1:IM)=TM;

S(IM+1)=(SM*D1+SS*D2)/DZ;
T(IM+1)=(TM*D1+TS*D2)/DZ; 
      
