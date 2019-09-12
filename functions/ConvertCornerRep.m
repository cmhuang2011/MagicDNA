function CornerRep = ConvertCornerRep( BaseArray )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
RemainInd = zeros(size(BaseArray ,1) ,1) ;
RemainInd(1)=1 ; RemainInd(end)=1 ;

for k=2: length(RemainInd)-1
    
    if BaseArray(k,1)~= BaseArray(k+1,1)
        RemainInd(k)=1 ;
    elseif BaseArray(k-1,1)~= BaseArray(k,1)
        RemainInd(k)=1 ;
    elseif abs(diff(BaseArray(k:k+1,2) )) >2
        RemainInd(k)=1 ;
    elseif abs(diff(BaseArray(k-1:k,2) )) >2
        RemainInd(k)=1 ;
    end
    
end
% [size(RemainInd)  size(BaseArray)  ]

% BaseArray
CornerRep= BaseArray(RemainInd==1 ,:) ;

% BaseRout = interpolateBase( CornerRep ) ;
% if  sum(sum(BaseArray==BaseArray)) == numel(BaseArray)
% sdfs=3;
% else
%     sdfsf=333
% end

% opposite :  BaseRout = interpolateBase( CornerRoute )
end

