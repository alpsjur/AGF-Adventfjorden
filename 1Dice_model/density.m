function rho = density(S,T)
% function rho = density(S,T)
%          returns rho (kg/m^3) as single value or row vector
%          S (psu) and T (C) can be vectors

if size(S,2)~=1; S=S'; end
if size(T,2)~=1; T=T'; end

C1=999.842594;
C2=  6.793952e-2;
C3= -9.095290e-3;
C4=  1.001685e-4;
D1=  0.824493;
D2= -4.0899e-3;
D3=  7.6438e-5;
D4= -8.2467e-7;
E1= -5.72466e-3;
E2=  1.0227e-4;
E3=  4.8314e-4;
                                  
ROW=C1+T.*(C2+T.*(C3+C4*T));
ROS=S.*(D1+T.*(D2+T.*(D3+D4*T)));

rho=ROW+ROS+S.*(sqrt(S).*(E1+E2*T)+E3*S);

