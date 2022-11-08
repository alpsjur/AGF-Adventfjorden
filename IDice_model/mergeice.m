function [Nmerge,STVmerge]=mergeice(Nmerge,STVmerge0);
% function mergeice(Nmerge,STVmerge);
%
% Merge two most similar ice classes

%          CALL MERGEICE(NRIDGE,STVR)
%      SUBROUTINE MERGEICE(N,STV)
%      REAL*8 STV(4,0:241),DHM,DH,AT,HE1,HE2,TMEL,A1,BS(0:25)
%      REAL*8 RAI,RAS,RAW,CPA,CPW,CPI,RLN,VFS,RLI,RKI,RKW,RKS
%      COMMON/C/ BS
%      COMMON/ICE/ RAI,RAS,RAW,CPA,CPW,CPI,RLN,VFS,RLI,RKI,RKW,RKS
%      INTEGER I,J,IM,N

global BS
global CPI RLI

TMEL=-BS(15);
STVmerge=STVmerge0;

%DHM=STVmerge(2,2)-STVmerge(1,2);

% Find smallest DH (ice thickness difference) between ice classes 
% Careful to check if there is open water in the first ice class, which will be
% the case for the Nlevel case, but not for the Nridge case

if STVmerge(1,2)==0, % open water present
  DHs=STVmerge(3:Nmerge,2)-STVmerge(2:(Nmerge-1),2);
  [DH,iDH]=min(DHs);
  if length(DH)>1; DH=DH(1); end
  IM=iDH+1; % e.g., if iDH=1, then we want DH between ice class 2 (IM) and 3 (IM+1)
else
  DHs=STVmerge(2:Nmerge,2)-STVmerge(1:(Nmerge-1),2);
  [DH,iDH]=min(DHs);
  if length(DH)>1; DH=DH(1); end
  IM=iDH; % e.g., if iDH=1, then we want DH between ice class 1 (IM) and 2 (IM+1)
end

% check which ice classes to merge (i.e., which ice classes have smallest
% thickness difference)
%if DH<DHM | (DH>DHM & STVmerge(1,2)<=0), % don't want to merge open water
%  IM=iDH+1; % e.g., if iDH=1, then we want DH between ice class 2 (IM) and 3 (IM+1)
%elseif DH>DHM & STVmerge(1,2)>0,
%  IM=1; % if DH>DHM and there's no open water, we want DH between first two ice classes (index 2 and 3)
%end

% disp(['Nmerge=',num2str(Nmerge), ', DH=',num2str(DH),', IM=',num2str(IM)]);

ATm=STVmerge(IM,1)+STVmerge(IM+1,1);
HE1=STVmerge(IM,1)*STVmerge(IM,2)*(CPI*STVmerge(IM,4)+RLI*TMEL/STVmerge(IM,4));
HE2=STVmerge(IM+1,1)*STVmerge(IM+1,2)*(CPI*STVmerge(IM+1,4)+RLI*TMEL/STVmerge(IM+1,4));
STVmerge(IM,2)=(STVmerge(IM,1)*STVmerge(IM,2)+STVmerge(IM+1,1)*STVmerge(IM+1,2))/ATm;
STVmerge(IM,3)=(STVmerge(IM,1)*STVmerge(IM,3)+STVmerge(IM+1,1)*STVmerge(IM+1,3))/ATm;
A1=(HE1+HE2)/(CPI*STVmerge(IM,2)*ATm);
STVmerge(IM,4)=A1/2-sqrt(A1*A1/4-BS(14));
STVmerge(IM,1)=ATm;

STVmerge((IM+1):(Nmerge-1),:)=STVmerge0((IM+2):Nmerge,:);
STVmerge=STVmerge(1:(Nmerge-1),:);
Nmerge=Nmerge-1;

%plot(1:IM,STVmerge(1:IM,2),'bx-',(IM+1):Nmerge,STVmerge((IM+1):Nmerge,2),'ro-',(IM+2):Nmerge,STVmerge0((IM+2):Nmerge,2),'g.-');
%title(['IM=',num2str(IM),', thickness(IM:IM+1)=',num2str(STVmerge(IM:(IM+1),2)')]);
%drawnow

