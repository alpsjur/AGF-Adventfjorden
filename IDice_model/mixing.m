function [HM]=mixing(DT,G,F,RHC,RM0,UST,BF,HSNO,HICE,HM)
% Simulating mxing at the bottom of the mixed layer (No diffusion further down)
% Mixing entrained water and calculating new mixed layer thickness and salinity.
% If advection is added, and should also go into the mixed layer,
% this needs to be added, see MIXING2 subroutine in
% ngfls/modeller/1Dice_snow/Barnilsen/1DICE_snow.f 
% 
% HM is mixed layer depth, given initially by ocean column values
% Bottom fix added by Lars H. Smedsrud, June 2016 
% Thus convection can reach Bottom, HMAX, and salinity will increase for whole column.

% **************************************************************************
%                                               * -----Surface -----------
%                                               * T(1) is from Z=0 to Z=0.5 DZ
%                                               * T(2) is 0.5 - 1.5 DZ
%                                               * T(3) is from 1.5 - 2.5 DZ
%                                               * _______________________  *
%      MIXING THE ENTRAINED WATER AND CALC.     *        |                 *
%      NEW MIXED LAYER THICKNESS AND SALINITY.  * _______|dz_____________  *
%                                               * -------|--------|D1----  *
%                                               * _______|dz__________|D2  *
% **************************************************************************

global S T SS TS
global DZ HMAX IMAX
global BS

%  CALL MIXING2(DT,IC,G,F,RHC,RM0,UST,BF,HSNO,HICE,HM,HMON,HEK,PAR) 
% SUBROUTINE MIXING2(DT,IT,G,F,RHC,RM0,UST,BF,HSNO,HICE,HM,HMON,HEK,PAR)

  HW = HM-BS(19)*HICE-BS(16)*HSNO;  % height of water layer
  IM = ceil((HM-0.5*DZ)/DZ);        % number of whole layers inside HM (original: IM=INT((HM-0.5*DZ)/DZ))
  D1 = HM-(IM-0.5)*DZ;              % upper delta Z, inside HM 
  D2 = DZ-D1;                       % lower delta Z, below HM 

  if BF>0  % Positive bouyancy, freshening/warming
      HMON = 2.*RM0*UST*UST*UST/BF+HICE*BS(19);      % Monin-Obukhov
      HEK  = RHC*UST/F+HICE*BS(19);                  % Ekman
      
      if HM>HMON || HM>HEK 
         HMN  = min(HMON,HEK);  % New mixed layer set by smallest value
         
         % Old subroutine "RETREAT" added here 
         IMN=ceil(HMN/DZ);
         D1=HM-(IM+.5)*DZ;
         D1N=HMN-(IM+.5)*DZ;
         D2=DZ-D1;

         if IMN >= IM                  % SAME GRID LEVEL
           SS=((D1-D1N)*S(1)+D2*SS)/((D1-D1N)+D2);
           TS=((D1-D1N)*T(1)+D2*TS)/((D1-D1N)+D2);
         else                          % NEW GRID LEVEL
           SS=S(1);
           TS=T(1);
         end
         % End of previous subroutine RETREAT   
         
         HM=HMN; return,  % Return out of mxing.m function
      else
         HMN=min(HMON,HEK);    % New mixed layer is deeper, and set by smallest value
         DH=HMN-HM;            % Mixed layer increase
         IMN=ceil(HMN/DZ);    % Mixed layer increased below next DZ?
         
         if IMN==IM           % New mixed layer not across next DZ
            SM=(HW*S(1)+DH*SS)/(HW+DH);
            TM=(HW*T(1)+DH*TS)/(HW+DH);
         else                  % New mixed cross into next DZ
            SM=(HW*S(1)+D1*SS)/(HW+D1);
            TM=(HW*T(1)+D1*TS)/(HW+D1);
            HWI=HW+D1;
            I=IM+2;
            while I < IMN
              SM=(HWI*SM+DZ*S(I))/(HWI+DZ);
              TM=(HWI*TM+DZ*T(I))/(HWI+DZ);
              HWI=HWI+DZ;
              I=I+1;
            end
            DH=HMN-(HWI+(HM-HW));  % DH set to new mixed layer - (mixed to now+ice) 
            SM=(HWI*SM+DH*S(I))/(HWI+DH);
            TM=(HWI*TM+DH*T(I))/(HWI+DH);
            DH=HMN-HM;             % Added by LHS as bottom fix, new mixed layer - old
         end
          HM=HM+DH;            % New mixed layer depth
          D1=HM-(IMN+.5)*DZ;   % New upper fraction 
          D2=DZ-D1;            % New lower fraction
            for I=1,IMN;                       
              S(I)=SM;          % Update HM Salt
              T(I)=TM;          % Update HM Temp
            end
          S(IMN+1)=(D1*S(1)+D2*SS)/DZ;  % Update one level down Salt
          T(IMN+1)=(D1*T(1)+D2*TS)/DZ;  % Update one level down Temp
          return  
      end
     
  else  % i.e. Boyancy flux is negative at this DT   LHS
      EKIN=DT*1000.*(RM0*UST*UST*UST-0.05*BF*HW/2.);    % AVAIL. KINETIC ENERGY
  end

  % Mixed Layer is loosing bouyancy and deepening below here ...
    RAM=density(S(1),T(1));    % density of upper cell
    RAS=density(SS,TS);        % mean density of mixed layer 
    EWHO=D2*G*(RAS-RAM)*HW/2.0;  % Energy needed for mxinig entire D2

     if EKIN < EWHO   %  ! i.e. mixing part of next DZ LHS
        PAR=1;
        
        if IM >= IMAX-1     % Bottom solution  LHS jan 2008
          DH=HMAX-HM;
%          WRITE(*,*) 'Mixing to bottom, one cell down, day',IT
          SM=(HW*S(1)+DH*SS)/(HW+DH);    %  Mixing to bottom
          TM=(HW*T(1)+DH*TS)/(HW+DH);    %  mixed T
          HM=HM+DH;                      %  New mixed layer H
          HW=HW+DH;
        else
          DH=EKIN*2./(HW*G*(RAS-RAM));
          SM=(HW*S(1)+DH*SS)/(HW+DH);
          TM=(HW*T(1)+DH*TS)/(HW+DH);
          HM=HM+DH;
          HW=HW+DH;
          D1=HM-(IM+.5)*DZ;
          D2=DZ-D1;
        end
        % Calculating new mixed layer values, advection must be added here if implemeted     
        for I=1,IM;                  
            S(I)=SM;
            T(I)=TM;
        end
     if IM >= IMAX-1      % Bottom solution  LHS jan 2008
        S(IM+1)=S(1);     % All mixed, no layer below
        T(IM+1)=T(1);     % D2 is below bottom 
     else   % Normal  
        S(IM+1)=(D1*S(1)+D2*SS)/DZ;
        T(IM+1)=(D1*T(1)+D2*TS)/DZ;
     end  % Bottom  test
        
     else    %  i.e. mixing all of next DZ  LHS
        PAR=2;
        
     if IM >= IMAX-1 %  Bottom solution  LHS jan 2008
        DH=HMAX-HM;   % bottom 
     else 
        DH=D2;
     end   
        SM=(HW*S(1)+DH*SS)/(HW+DH);
        TM=(HW*T(1)+DH*TS)/(HW+DH);
        HM=HM+DH;
        HW=HW+DH;
        IC=0;
        EUSE=EWHO;    		            % USED ENERGY
     while EUSE < EKIN,
        PAR=PAR+1;
     
     if IM+IC >= IMAX-1   % Bottom solution  LHS jan 2008
        DH=HMAX-HM;        % Mixing to bottom     
        SM=(HW*SM+DH*SS)/(HW+DH);    %  Mixed S to bottom
        TM=(HW*TM+DH*TS)/(HW+DH);    %  mixed T
        HM=HM+DH;              %  New mixed layer H
        HW=HW+DH;
        
        EUSE=EKIN+0.1;        % Rest of energy used, NOT used to mix further
        IC=IC+1;
     else                             % Normal solution 
             
        RAM=density(SM,TM);                   % density of mixed layer
        RAS=density(S(IM+2+IC),T(IM+2+IC));   % density below 
        EWHO=DZ*G*(RAS-RAM)*HW/2.0;
        if EKIN < EUSE+EWHO
            DH=(EKIN-EUSE)*2./(HW*G*(RAS-RAM));
        else
            DH=DZ;  % Mixing whole next delta Z as well
        end
        SM=(HW*SM+DH*S(IM+2+IC))/(HW+DH);
        TM=(HW*TM+DH*T(IM+2+IC))/(HW+DH);
        HM=HM+DH;
        HW=HW+DH;
        EUSE=EUSE+EWHO;
        IC=IC+1;   % Counting DZ's
      end
    end
     
      IM=ceil((HM-0.5*DZ)/DZ);   % Subtract surface layer 0.5*DZ
      D1=HM-(IM-0.5)*DZ;         % upper delta Z, inside HM  
      D2=DZ-D1;
          for I=1,IM;                      
            S(I)=SM;
            T(I)=TM;
          end
        SS=S(IM+1);        % New salt "split" value, Old values below HM                           
        TS=T(IM+1);        % Temp split value
        
        if IM+IC >= IMAX-1  % Bottom solution  LHS jan 2008
          S(IM+1)=S(1);     % All mixed, no layer below
          T(IM+1)=T(1);     % D2 is below bottom 
        else         % Normal  
          S(IM+1)=(D1*S(1)+D2*SS)/DZ;  
          T(IM+1)=(D1*T(1)+D2*TS)/DZ;
        end   % Bottom test
     end  %  End of mixing all below DZ  LHS
end  % return
    
    
 