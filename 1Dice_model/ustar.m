function [UST,UICE,UFRA]=ustar(RAW,RAA,STV,AICE,WIND,MWIN);
% function [UST,UICE,UFRA]=ustar(RAW,RAA,STV,AICE,WIND,FWIN,MWIN);
% Calculate friction velocity
% Lars H removed variable in call FWIN which was constant from table =5.0

global BS

%      CALL USTAR(RAW,RAA,STV,AICE,WIND,FWIN,MWIN,
%     >                 UST,UICE,UFRA)
%      SUBROUTINE USTAR(RAW,RAA,STV,AICE,WIND,FWIN,MWIN,
%     >                 UST,UICE,UFRA)
%      REAL*8 BS(0:25),WIND,AICE,UICE,UST,MWIN,FWIN,UFRA
%      REAL*8 RAW,RAA,CDW,USTI,USTW,STV(4,0:241)
%      COMMON/C/ BS

if STV(1,2)<=0, % if ice thickness < = 0
  if MWIN<6,    % MWIN = wind forcing FTAB(17,:) calculcated in main program
    CDW = 1.1e-3;
  else,
    CDW = .61e-3+MWIN*.063e-3;
  end
  USTW = sqrt(RAA/RAW*CDW)*MWIN*STV(1,1);
  USTI = BS(21)*AICE*MWIN*(1.-STV(1,1)); % BS(21) is ice/water drag coefficient 
  UST  = USTW+USTI; 
else,
  UST=BS(21)*AICE*MWIN; 
end

UICE=AICE*WIND;

%UFRA=-AICE*FWIN  % UFRA also seems not to be in use ...
UFRA=-0.02*5.0;   % Lars H. replaced FWIN with constant value in table 5.0 
