function [BF]=buoyflux(G,WIME,SICE,FWAT)
% Calculate surface buoyancy flux

global BS
global S T SS TS RA

RA(1)=density(S(1),T(1));
[CT,CS]=densder(S(1),T(1));

ALF=CT/RA(1);
BET=CS/RA(1);

BF=G*(-ALF*FWAT/BS(12)-BET*(BS(19)*WIME*(S(1)-SICE)));

end



