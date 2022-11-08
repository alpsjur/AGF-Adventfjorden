function [CT,CS]=densder(S,T)
% function [CT,CS]=densder(S,T);
      
C1=6.793952e-2; C2=-9.09529e-3; C3=1.001685e-4; 

D1=-4.0899e-3;  D2=7.6438e-5;   D3=-8.2467e-7; 

E1=1.0227e-4;   E2=.824493;     E3=-4.0899e-3; E4=7.6438e-5; 
E5=-8.2467e-7;  E6=-5.72466e-3; E7=1.0227e-4;  E8=4.8314e-4; 
E9=3/2; 

C22=2*C2; C33=3*C3; D22=2*D2; D33=3*D3; E82=2*E8;

SQS=sqrt(S);
A=C1+T*(C22+T*C33);
B=S*(D1+T*(D22+T*D33));
CT=A+B+S*SQS*E1;

CS=E2+T*(E3+T*(E4+T*E5))+E9*SQS*(E6+E7*T)+E82*S;

