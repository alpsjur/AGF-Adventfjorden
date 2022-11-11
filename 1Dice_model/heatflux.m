function [FHET,TSF,OPGR,FWAT] = heatflux(N,ENR,HTH,TFFW,BCON,ALF,STV,I0,IW,...
                               FR,FL,TA,WIND,UICE,RH,ASN)
                           
%      SUBROUTINE HEATFLUX(N,HM,DT,CE,ENR,HTH,TFFW,BCON,ALF,AS,STV,I0,
%     >         IW,ALCU,FR,FL,TA,WIND,UICE,RH,ASN,FHET,TSF,OPGR,FWAT)
%                                     
%     *** FHET ****(12 FLUXES OF HEAT FOR ALL ICE CLASSES )*********************
%                                                                               
%     FHET(1,N) Short Wave IN                FHET(2,N) Short Wave NET IN  
%     FHET(3,N) Short Wave Penetration(F0)   FHET(4,N) Long Wave IN            
%     FHET(5,N) Long Wave OUT                FHET(6,N) Sensible heat
%     FHET(7,N) Latent heat                  FHET(8,N) Net Upward Flux 'Air' from surface (Fa) 
%     FHET(9,N) Ice Upward Flux 'Snow' (Fs)  FHET(10,N) Flux to ice 'Bottom'(Fb) 
%     FHET(11,N) Flux from Water(Fw)         FHET(12,N) Total flux to ice (F0+Fa)
%             
% ***************************************************************************
%      "Illustration of the flows" - Johanna Nilsson, Gothenburg,
%
%                            Ts      /\
%    _____                ____.______||Fa___________   _
%         |              |     \                    |  |          |  
%     snow|              |      \    /\      snow   |  h          |
%    _____|  ||Fa        |_||____\___||Fs___________| _|_ _  _    |
%         |~~\/~~~~~~~~~~| \/(F0) \_                |     |  |    V z-direction
%         |  ||          |          \_ T            |     |  |
%      ice|  \/          |            \.       ice  |     H  |
%         |  (1-aice)Fsw |    /\Fb      \   Tf      |     |  |
%    _____|              |____||_________\._________|    _|_ |
%                             /\          |         |        |   
%                             ||Fw        |         |        Hm
%                                         |         |        |
%                                         .Tm       |        |
%    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _|       _|_
%            ||                                     |
%            \/(1-aice)Fsw*Iw*exp(-ae*Hm)           |    
%                                                   |        
% ****************************************************************************                                              

global S T
global AS BS
global RKI RKS

TB=TA+TFFW;   % Changes to Kelvin
CTA=RH*TB*(TB*(TB*(AS(1)*TB+AS(2))+AS(3))+AS(4))+AS(5)*(RH-1);
IICE=1;

% ****  CALC OCEAN HEAT FLUX
TF=-BCON*S(1);                % Freezing temp of sea water
FWAT=BS(5)*UICE*(TF-T(1));    % Ocean Heat Flux
OPGR=0; % false

if (STV(1,2)<=0), % i.e. open water is present
  IICE=2;
  TX1=T(1)+TFFW;
  TX2=TX1*TX1;
  TX3=TX2*TX1;
  TX4=TX3*TX1;
  ALFA=ALF(3);
  FHET(1,1) = FR;                   % Short wave radiation in
  FHET(2,1) = (1-ALFA)*FHET(1,1);   % Net short wave in
  FHET(3,1) = IW*FHET(2,1);         % Short wave penetration
  FHET(4,1) = FL;                   % Long wave in
  FHET(5,1) = -BS(2)*TX4;           % Long Wave out
  FHET(6,1) = BS(3)*(TA-T(1))*WIND; % Sensible
  FHET(7,1) = BS(4)*(CTA-(AS(1)*TX4+AS(2)*TX3+AS(3)*TX2+AS(4)*TX1))*WIND;  % latent 
  FHET(8,1) = FHET(2,1)*(1-IW)+FHET(4,1)+FHET(5,1)+FHET(6,1)+FHET(7,1);  % Net upper surf
  FHET(9,1) = FHET(8,1);   % Ice upper surf
  FHET(10,1) = FHET(8,1);  % Ice lower surf

  if (FHET(8,1) < FWAT),  % open water ice growth
    FHET(11,1)=FWAT;
    OPGR=1; % true
  else 
    % no open water production
    FHET(11,1)=FHET(8,1);
  end

  TSF(1)=T(1);
  FHET(12,1)=FHET(8,1)+FHET(3,1);

end %if (STV(1,2)<=0), % ice thickness

TSURF=0;

% Calculate Heat Flux through present Ice Cover 
for I=IICE:N,

  RKT=RKI+BS(10)/STV(I,4);

  if (STV(I,3)<=0),  % No snow
    ALFA=min(ALF(1)*STV(I,2)^ALF(2)+ALF(3),ALF(4));
    
    EI0=I0;
    EEE=BS(2);
  else            % Snow present
    ALFA=ASN;
    EI0=0;
    EEE=BS(1);
  end

  FHET(1,I)=FR;  % Short wave
  FHET(2,I)=(1-ALFA)*FHET(1,I);
  FHET(4,I)=FL; 
  
  if (STV(I,2)>HTH),  % Ice thicker than HTH 

    FHET(3,I)=EI0*FHET(2,I);
    RJAM=RKT*RKS/(RKS*STV(I,2)/2.+RKT*STV(I,3));
    % Calculate T suf based on incoming fluxes 
    % disp(['   T surf: ',num2str([TSURF])]);
    
    [FHET(5,I),FHET(6,I),FHET(7,I),FHET(9,I),TSURF] = ...  
        surftemp(ENR,TFFW,RJAM,EEE,STV(I,4),EI0,TA,CTA,WIND,FHET(2,I),FHET(4,I),TSURF);
     FHET(10,I)=2.*RKT*(STV(I,4)-TF)/STV(I,2);
  else % if (STV(I,2)<HTH
    EI0=0;
    FHET(3,I)=EI0*FHET(2,I);
    RJAM=RKT*RKS/(RKS*STV(I,2)+RKT*STV(I,3));
    % Calculate T suf based on incoming fluxes 
    [FHET(5,I),FHET(6,I),FHET(7,I),FHET(9,I),TSURF] = ...   
        surftemp(ENR,TFFW,RJAM,EEE,T(1),EI0,TA,CTA,WIND,FHET(2,I),FHET(4,I),TSURF);
     FHET(10,I)=FHET(9,I);
  end % HTH check 

  TSF(I)=TSURF;
  if (TSF(I)<0),
      FHET(8,I)=FHET(9,I);
  else 
      FHET(5,I)=-EEE*BS(16);
      FHET(6,I)=BS(3)*WIND*(TA-TSF(I));
      FHET(7,I)=BS(4)*WIND*(CTA-BS(18));
      FHET(8,I)=(1-EI0)*FHET(2,I)+FHET(4,I)+FHET(5,I)+FHET(6,I)+FHET(7,I);
  end
  FHET(11,I)=FWAT;
  FHET(12,I)=FHET(8,I)+FHET(3,I);

end 
 % Sum up the ocean heat flux over all ice classes 
 FWAT=0.;
 for I=IICE:N,
     FWAT=FWAT+STV(I,1)*FHET(11,I);
 end

%disp(['     heatflux: ',num2str([IICE FR FL])]);
%disp(FHET(:,1));
