function IndS =QuerryXoverInd_BCB( Xover,BaseByBaseRouting )
% debugging use, multiscaffold split
%   Detailed explanation goes here

IndS=zeros(size(Xover,1) ,size(Xover,2)/3 ) ;

for i= 1 : size(IndS,1)
    for j= 1 : size(IndS,2)
     [~,bi] =ismember(Xover(i,3*j-2:3*j),  BaseByBaseRouting,'rows' );

        IndS(i,j) =bi ;
        
    end
end
end

