function N= countCM( X,Arr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% X
X=X(X>0) ;
Xa=[X; max(Arr) ] ;

QQ=accumarray(Xa,1) ;
QQ(end) = QQ(end)-1 ;
N= QQ(Arr) ;
% m=Arr ;

end

