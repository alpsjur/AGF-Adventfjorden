function [STV,HICE,HSNO,HM,WIME,WSME,DHS,DHU,DHL]=icegrowth(N,HM,SICE,OPGR,ACOFF,EXPO,IPEN,EXIP,...
                                                            GFAC,STV,TSF,FHET,HICE,HSNO)
% DHU = change in ice thickness at upper surface
% Uses: density, icetemp, npconv, radhet
%     SUBROUTINE ICEGROWTH(N,HM,SICE,OPGR,ACOFF,EXPO,IPEN,EXIP,
%     >           GFAC,STV,TSF,FHET,HICE,HSNO,WIME,WSME,DHS,DHU,DHL)

global S T SS TS
global DT DZ BS
global RAI CPW CPI VFS RLI

IICE=1;
TMEL=-BS(15);
TFRE=-BS(15)*S(1)/SICE;
EXHM=exp(-ACOFF*HM);
C2=RLI*TMEL/CPI;
HWME=HM-BS(19)*HICE-BS(16)*HSNO;

% *** OPEN WATER SPECIAL TREATEMENT
if (STV(1,2) <= 0), % ice thickness < 0
  % disp(['     icegrowth: open water found']);
  IICE=2;     % marks open water
  DHU(1)=0;
  %DHS(1)=0;           % AL: old code
  DHS(1)=-STV(1,3);    % AL: melt snow in open water, see line 97->, do something like this?
  if (OPGR),
    % disp(['     icegrowth: OPGR==1 (open water production)']);
    STV(1,4)=0;
    DTI(1)=TFRE;
    DHL(1)=(FHET(8,1)-FHET(11,1))*DT/(RAI*((CPI*TFRE-RLI)*...
            (1.-TMEL/TFRE)+CPW*(TMEL-T(1))));
    DTW(1)=DT*(FHET(11,1)+FHET(3,1)*(1.-EXHM+EXIP))/(BS(12)*...
    (HM-BS(19)*DHL(1)));
    % disp(['     icegrowth: OPGR==0 FHET([11 3],1) DHL(1) EXHM DTW(1) ',num2str([FHET(11,1) FHET(3,1) DHL(1) EXHM DTW(1)])]);
  else
    % disp(['     icegrowth: OPGR==0 (no open water production)']);
    DHL(1)=0;
    DTI(1)=0;
    DTW(1)=DT*(FHET(8,1)+FHET(3,1)*(1.-EXHM+EXIP))/(BS(12)*HM);
  end
%  disp(['     icegrowth: debug a, T(1) = ',num2str(T(1))]);
%  disp(['     icegrowth: debug x, DTW(1), FHET([3 8],1)= ',num2str([DTW(1) FHET([3 8],1)'])]);
  [T]=radhet(DT,STV,HM,ACOFF,FHET,EXPO,EXHM,IPEN);
%  disp(['     icegrowth: debug b, T(1) = ',num2str(T(1))]);
end   % Open water

% loop through rest of ice classes
for I=IICE:N,
  if STV(I,3) <= 0,   % snow thickness <= 0
    DHET=FHET(8,I)-FHET(11,I)+FHET(3,I);  % Net heat goes to ice melt
  else
     DHET=-1e20;              % Goes to snow melt instead ...
  end

  if (DHET > -RAI*((CPI*STV(I,4)-RLI)*(1-TMEL/STV(I,4)))*STV(I,2)), % More heat than total ice melt?
    DHU(I)=-STV(I,2);      % melting all ice thicness
    DHL(I)=0;
    DHET=FHET(12,I)*DT+STV(I,2)*RAI*((CPI*STV(I,4)-RLI)*(1.-TMEL/STV(I,4))+CPW*(TMEL-T(1)));
    DTW(I)=DHET/(HM*BS(12));  % remaining heat goes to warm ocean

  else  			   %!!!NEW

    % *** LOWER SURFACE GROWTH/MELT

    DHET=(FHET(11,I)-FHET(10,I))*DT; % flux from water (Fw) - flux to bottom of ice (Fb)
      %if I==1,
        %disp(['     icegrowth lower surface melt DHET= ',num2str(DHET)]);
        %disp(['     icegrowth lower surface melt FHET(10:11,1)= ',num2str(FHET(10:11,1)')]);
      %end
    if (DHET> 0), %GROWING LOWER SURF
      DHL(I)=-DHET/(RAI*((CPI*TFRE-RLI)*(1-TMEL/TFRE)+CPW*(TMEL-T(1))));
        %disp(['     icegrowth lower surface growth: ',num2str(DHL(1))]);
      LGR=1; % true  - lower surface growth?
    else
      DHL(I)=-DHET/(RAI*((CPI*STV(I,4)-RLI)*(1-TMEL/STV(I,4))));
%      if I==1,
%        disp(['     icegrowth lower surface melt: ',num2str([-DHET  RAI*(CPI*STV(1,4)-RLI)*(1-TMEL/STV(1,4))])]);
%        disp(['     icegrowth lower surface melt: ',num2str([-DHET/(RAI*(CPI*STV(1,4)-RLI)*(1-TMEL/STV(1,4)))])]);
%        disp(['     icegrowth lower surface melt: ',num2str(DHL(1))]);
%      end
      LGR=0; % false  - lower surface growth?
    end

    %*** UPPER SURFACE MELT

    UME=0; % false
    SME=0; % false
    DHS(I)=0;
    DHU(I)=0;
    if (TSF(I) >= 0),   % Temp at surface > zero?
      DHET=(FHET(8,I)-FHET(9,I))*DT; % net upward surface flux (Fa)- ice upward flux 'snow' (Fs)
      if (STV(I,3) > 0.),    % Snow thickness > zero
        DHS(I)=-DHET/VFS;
        SME=1; % true  - some snow melts
        if (STV(I,3)+DHS(I) < 0),  % All snow melts
          DHS(I)=-STV(I,3);
          DHET=DHET-STV(I,3)*VFS;
          DHU(I)=DHET/(RAI*((CPI*STV(I,4)-RLI)*(1.-TMEL/STV(I,4))));
          UME=1; % true - all 'ultimate' snow melts
        end
      else               % snow thickness has been zero last time step
        DHU(I)=DHET/(RAI*(CPI*STV(I,4)-RLI)*(1.-TMEL/STV(I,4)));  % C  rhoi(cpi*Ti-Li)(1-Tmelt/Tfreez)
        UME=1; % true
      end
    end
    DHWU=-BS(19)*DHU(I);
    DHWL=-BS(19)*DHL(I);
    DHWS=-BS(16)*DHS(I);
    HW=HM-BS(19)*STV(I,2)-BS(16)*STV(I,3);

    %*** ICE AND WATER TEMPERATURE CHANGE
    if (STV(I,2)+DHU(I)+DHL(I) < 0),  %MELTING ALL ICE
      DHU(I)=-STV(I,2);
      DHL(I)=0;
      DHET=FHET(12,I)*DT-STV(I,3)*VFS+STV(I,2)*RAI*((CPI*STV(I,4)-RLI)*(1.-TMEL/STV(I,4))+CPW*(TMEL-T(1)));
      DTW(I)=DHET/(HM*BS(12));
      DTI(I)=0.0;  % LHS added - to prevent crashing when all ice melts

    else

      DHET=(FHET(9,I)-FHET(10,I)+FHET(3,I))*DT;   % heat leaving ice W/m² * time [s] = [J/m²]
       if (UME==1 && SME==1),  %SNOW MELT AND UPPER SURF. MELT
          if (LGR),
              TIN=icetemp(STV(I,4),DHET,STV(I,2)+DHU(I));
              C1=((CPI*TFRE+RLI*TMEL/TFRE)*DHL(I)+(CPI*TIN+...
                 RLI*TMEL/TIN)*(STV(I,2)+DHU(I)))/...
                 (CPI*(STV(I,2)+DHL(I)+DHU(I)));
              DTI(I)=C1/2-sqrt(C1*C1/4.-C2)-STV(I,4);
              DTW(I)=(FHET(11,I)*DT/BS(12)-T(1)*DHWS+(TMEL-T(1))*DHWU)/(HW+DHWL+DHWU+DHWS);
          else % LGR=0
              TIN=icetemp(STV(I,4),DHET,STV(I,2)+DHU(I)+DHU(I));
              DTI(I)=TIN-STV(I,4);
              DTW(I)=(FHET(11,I)*DT/BS(12)+(TMEL-T(1))*(DHWL+DHWU)-T(1)*DHWS)/(HW+DHWL+DHWU+DHWS);
          end
       elseif SME==1 && UME==0,     % snow melts partly
          if (LGR),   % still growth at ice bottom
              TIN=icetemp(STV(I,4),DHET,STV(I,2));
              C1=((CPI*TFRE+RLI*TMEL/TFRE)*DHL(I)+(CPI*TIN+RLI*TMEL/TIN)*STV(I,2))/(CPI*(STV(I,2)+DHL(I)));
              DTI(I)=C1/2-sqrt(C1*C1/4.-C2)-STV(I,4);
              DTW(I)=(FHET(11,I)*DT/BS(12)-T(1)*DHWS)/(HW+DHWL+DHWS);
          else        % no growth at bottom
              TIN=icetemp(STV(I,4),DHET,STV(I,2)+DHL(I));
              DTI(I)=TIN-STV(I,4);
              DTW(I)=(FHET(11,I)*DT/BS(12)+(TMEL-T(1))*DHWL-T(1)*DHWS)/(HW+DHWL+DHWS);
          end
       elseif SME==0 && UME==1,         % snow melts totally 'ultimate melt'
          if (LGR),
            TIN=icetemp(STV(I,4),DHET,STV(I,2)+DHU(I));
            C1=((CPI*TFRE+RLI*TMEL/TFRE)*DHL(I)+(CPI*TIN+RLI*TMEL/...
               TIN)*(STV(I,2)+DHU(I)))/(CPI*(STV(I,2)+DHL(I)+DHU(I)));
            DTI(I)=C1/2-sqrt(C1*C1/4.-C2)-STV(I,4);
            DTW(I)=(FHET(11,I)*DT/BS(12)+(TMEL-T(1))*DHWU)/...
                   (HW+DHWL+DHWU);
          else  % now growth below ice
            TIN=icetemp(STV(I,4),DHET,STV(I,2)+DHL(I)+DHU(I));
            DTI(I)=TIN-STV(I,4);
            DTW(I)=(FHET(11,I)*DT/BS(12)+(TMEL-T(1))*(DHWL+DHWU))/...
                   (HW+DHWL+DHWU);
          end
       elseif SME==0 && UME==0,         % no snow melt
          if (LGR),   % Growth of ice
            TIN=icetemp(STV(I,4),DHET,STV(I,2));
            C1=((CPI*TFRE+RLI*TMEL/TFRE)*DHL(I)+...
               (CPI*TIN+RLI*TMEL/TIN)*STV(I,2))/(CPI*(STV(I,2)+DHL(I)));
            DTI(I)=C1/2.-sqrt(C1*C1/4.-C2)-STV(I,4);
            DTW(I)=FHET(11,I)*DT/(BS(12)*(HW+DHWL));
          else        % No growth
            TIN=icetemp(STV(I,4),DHET,STV(I,2)+DHL(I));
            DTI(I)=TIN-STV(I,4);
            DTW(I)=(FHET(11,I)*DT/BS(12)+(TMEL-T(1))*DHWL)/(HW+DHWL);
          end
      end
     end
  end
end % for ice all ice classes loop


%     *** CALC. TOTAL CHANGES INCLUDING M.L. SALINITY

%DHICE=0;
%DHSNO=0;
%HICE=0;
%HSNO=0;
%TM=0;

STV(1:N,2)=STV(1:N,2)+DHL(1:N)'+DHU(1:N)'; % Ice Thickness
STV(1:N,3)=STV(1:N,3)+DHS(1:N)';           % snow thickness
STV(1:N,4)=STV(1:N,4)+DTI(1:N)';           % Ice temp
HICE =sum(STV(1:N,1).*STV(1:N,2));         % New mean ice thickness
HSNO =sum(STV(1:N,1).*STV(1:N,3));         % New mean snow thickness
DHICE=sum(STV(1:N,1).*(DHL(1:N)+DHU(1:N))'); % Total change in ice thickness
DHSNO=sum(STV(1:N,1).*DHS(1:N)');            % Change in snow thickness

%for I=1:N,
%  TM=TM+STV(I,1)*(HM-BS(19)*STV(I,2)-BS(16)*STV(I,3))*(T(1)+DTW(I));
%end
TM=sum(STV(1:N,1).*(HM-BS(19)*STV(1:N,2)-BS(16)*STV(1:N,3)).*(T(1)+DTW(1:N)'))/(HM-BS(19)*HICE-BS(16)*HSNO);

%disp(['     icegrowth: debug n, STV(1,1:3), HM, DTW(1) = ',num2str([STV(1,1:3) HM DTW(1)])]);
%disp(['     icegrowth: debug e, TM, T(1) = ',num2str([TM T(1)])]);
%disp(['     icegrowth: debug e, HM, HICE, HSNO = ',num2str([HM HICE HSNO])]);
%disp(['     icegrowth: debug b, term1 = ',num2str([ HM BS(19)*STV(1,2)-BS(16)*STV(1,3)+BS(16)*STV(1,3)])]);
%disp(['     icegrowth: debug b, term2 = ',num2str([T(1) DTW(1) T(1)+DTW(1)])]);
%disp(['     icegrowth: debug b, term3 = ',num2str([HM-BS(19)*HICE-BS(16)*HSNO])]);

%TM1=STV(1:N,1).*(HM-BS(19)*STV(1:N,2)-BS(16)*STV(1:N,3)).*(T(1)+DTW(1:N)');
%TM2=(HM-BS(19)*HICE-BS(16)*HSNO)/N;

% plots some diagnostics
%figure(2);
%subplot(1,2,1);
%  plot(TM1,1:N,100*STV(1:N,1),1:N,HM-BS(19)*STV(1:N,2)-BS(16)*STV(1:N,3),1:N,T(1)+DTW(1:N),1:N);
%  title(sum(TM1)/(HM-BS(19)*HICE-BS(16)*HSNO)); drawnow;
%  legend('TM1','STV(:,1) %','term1','term2');
%subplot(2,2,2);
%  plot(T(1:100),-[1:100],'.-');
%  title('T'); axis([-7 3 -100 0]);
%  drawnow;
%subplot(2,2,4);
%  plot(T(100:451),-[100:451],'.-'); axis([-6 3 -451 -100]);
%  drawnow;

if (DHICE > 0),
  SM=(S(1)*HWME-SICE*BS(19)*DHICE*(1.-GFAC))/(HWME-BS(19)*DHICE*(1.-GFAC)-BS(16)*DHSNO);
  %disp(['     icegrowth a ',num2str([S(1) SM])]);
else
  SM=(S(1)*HWME-SICE*BS(19)*DHICE)/(HWME-BS(19)*DHICE-BS(16)*DHSNO);
  %disp(['     icegrowth b ',num2str([DHL(1) DHU(1) STV(1,1)])]);
  %disp(['     icegrowth b ',num2str([DHICE DHSNO])]);
  %disp(['     icegrowth b ',num2str([S(1) SM])]);
end

WIME=DHICE/DT;
WSME=DHSNO/DT;
IM=ceil((HM-0.5*DZ)/DZ);   % bottom level of mixed layer (original: IM=INT((HM-0.5*DZ)/DZ))
D1=HM-(IM-0.5)*DZ;         % vertical distance to mid-level (original: D1=HM-D1=HM-(IM+0.5)*DZ)

% set new mixed layer salinity and temperature
S(1:IM)=SM;
T(1:IM)=TM;

% Update layer below mixed layer
S(IM+1)=(SM*D1+SS*(DZ-D1))/DZ;
T(IM+1)=(TM*D1+TS*(DZ-D1))/DZ;

RAM=density(SM,TM); % new mixed layer density
RS=density(SS,TS);  % density of split layer directly below

% disp(['     icegrowth: debug h, T(1), T(IM+1), IM = ',num2str([T(1) T(IM+1) IM])]);
% disp(['     icegrowth: debug i, old T, new T = ',num2str([TS TM])]);
% disp(['     icegrowth: debug i, old S, new S = ',num2str([SS SM])]);
% disp(['     icegrowth: debug j, old density, new density = ',num2str([RS RAM])]);
if (RAM > RS),
  disp(['     icegrowth: Unstable mixed layer: old density ',num2str(RS),', new density ',num2str(RAM),' - adjust']);
  HM=npconv(HM,HSNO,HICE);
end

%if T(1)<-2; error('     icegrowth: Error = T(1) < -2.2 C'); end

