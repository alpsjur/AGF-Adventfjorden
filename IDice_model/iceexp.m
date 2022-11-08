function [N,STV,HSNO,HICE]=iceexp(AEX,N,STV)
%
% Export EXS part of each ice class
%      SUBROUTINE ICEEXP(AREA,AEX,DT,TM,N,STV,HSNO,HICE)
%       REAL*8 STV(4,0:241),AEX,EFA,TM,DT,EXT,AREA,HICE,HSNO
%       REAL*8 RAI,RAS,RAW,CPA,CPW,CPI,RLN,VFS,RLI,RKI,RKW,RKS
%       INTEGER N,I,J
%       COMMON/ICE/ RAI,RAS,RAW,CPA,CPW,CPI,RLN,VFS,RLI,RKI,RKW,RKS

global T DT    
global AREA  

EXT=AEX*DT/AREA;    % New concentration of open water due to export

if STV(1,2)>0             % No open water at this time
  N=N+1;
  tempSTV(2:N,:)=STV(1:(N-1),:);  % Create an extra ice class
  STV=tempSTV; 

  % GENERATE OPEN WATER
  STV(1,1)=EXT;  % New open water Concentration
  STV(1,2)=0;    % Ice thickness
  STV(1,3)=0;    % snow thickness
  STV(1,4)=T(1);   % ice temperature

  STV(2:N,1)=STV(2:N,1)-EXT*STV(2:N,1); % Reduce the other areas equally

  % disp([' iceexp: new open water created, total open water now [percent] :  ',num2str(STV(1,1)*100)]);

else                % Open water is present

  EFA=EXT/(1-STV(1,1));                   % Part of the ice area exported
  STV(1,1)=STV(1,1)+EXT;                  % Add to open water fraction 
  STV(2:N,1)=STV(2:N,1)-EFA*STV(2:N,1);   % Reduce the other areas equally

end 

HICE = STV(1:N,1)'*STV(1:N,2);  % Concentration * Ice Thickness
HSNO = STV(1:N,1)'*STV(1:N,3);  % Concentration * Snow Thickness

