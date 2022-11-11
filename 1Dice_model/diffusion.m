function [T,S]=diffusion(DT,DKF,HM)
% Simulating diffusion, or mixing below the mixed layer 
% Converted to Matlab from FORTAN subroutine INTMIX.f by Lars H. Smedsrud, Oct 2022 

% DT is time step [s]
% DKF is 'Diffusion Constant' set as parameter
% HM is mixed layer thickness

global S T
global DZ IMAX

IM = ceil((HM-0.5*DZ)/DZ);        % number of whole layers inside HM (original: IM=INT((HM-0.5*DZ)/DZ))

SL(IM+2)=S(IM+2)+DKF*(S(IM+3)-S(IM+2))*DT/(DZ*DZ); % S for first layer below HM
TL(IM+2)=T(IM+2)+DKF*(T(IM+3)-T(IM+2))*DT/(DZ*DZ); % T for first layer below HM
          
   for i=IM+3:IMAX-1  % Calculate new values below mixed layer
       SL(i)=S(i)+DKF*(S(i+1)-2.*S(i)+S(i-1))*DT/(DZ*DZ);
       TL(i)=T(i)+DKF*(T(i+1)-2.*T(i)+T(i-1))*DT/(DZ*DZ);
   end
   for i=IM+2:IMAX-1  % Update to new S and T values 
       S(i)=SL(i);
       T(i)=TL(i);
   end

