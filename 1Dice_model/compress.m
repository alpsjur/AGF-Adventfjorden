function [N,STV]=compress(N,STV)

% Merge open water classes if necessary

while (N>1 && (STV(1,2)==0 && STV(2,2)==0)), % two open water classes

  clear tempSTV tempIV

  disp(['     compress: more than one open-water class found - merging']);
  disp(['               # open-water classes: ',num2str(length(find(STV(:,2)==0)))]);

  tempSTV(1,:)      =[STV(1,1)+STV(2,1) STV(1,2:4)];  % Add open water fraction, rest fractions remains

  tempSTV(2:(N-1),:)=STV(3:N,:);  % Loose values for N=2
  STV=tempSTV;                    % Overwrite with new values
  N=N-1;                          % Cut one ice class

end

if size(STV,1)~=N; error('     compress: N does not equal length of STV'); end


