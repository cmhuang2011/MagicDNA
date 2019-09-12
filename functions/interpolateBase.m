function BaseRout = interpolateBase( CornerRoute )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% tic
BaseRout = zeros(50000, 2 )    ; % allocate
n=1 ;
for k= 1 :2: size(CornerRoute,1)-1  % fixed Apr 25
    PiPj=  CornerRoute(k:k+1,:) ;
    if PiPj(1,1)==PiPj(2,1)
    Add =  [PiPj(1,1)*ones(   abs(diff(PiPj(:,2)))+1 ,1) , linspace(PiPj(1,2) ,PiPj(2,2), abs(diff(PiPj(:,2)))+1)' ] ;
    na =size(Add,1) ;
    BaseRout(n:n+na-1, :) = Add;
    n=n+na ;
    end
end

BaseRout= BaseRout(sum(BaseRout,2)~=0 ,:) ;
% toc

% opposite : CornerRep = ConvertCornerRep( BaseArray )

end

