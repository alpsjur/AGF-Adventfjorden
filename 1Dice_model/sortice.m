function [STV]=sortice(STV);

% Function that sorts into increasing thickness classes. 

if issorted(STV(:,2))==0,

  % disp('     sort: Out-of-order ice classes found - sort');
  [vsort,isort]=sort(STV(:,2));
  tempSTV=STV(isort,:);
  STV    =tempSTV;
    
end

