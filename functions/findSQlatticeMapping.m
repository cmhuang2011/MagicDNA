function [ SQMapping ] = findSQlatticeMapping( r ,XYbound)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% tic

% x=-XYbound(1):r*sqrt(3):XYbound(1);
% x= 0:r*sqrt(3):XYbound(1) ;
x=sort(unique( [ 0:r*2:XYbound(1),  0:-2*r:-XYbound(1)] )) ; 



% y=-XYbound(2):3*r:XYbound(2);
% y=sort(unique( [ 0:3*r:XYbound(2),  0:-3*r:-XYbound(2)] ) ) ;
y=sort(unique( [ 0:r*2:XYbound(2),  0:-2*r:-XYbound(2)] )) ; 

% nrow=length(y);
% nColumn=length(x);
% SQcenter=zeros(nColumn*nrow,4);
% for i=1:length(x)
%     for j=1:length(y)
% 
%     cx=x(i);
%     cy=y(j) ;
%     k=(i-1)*nrow+j;
%     SQcenter(k,:)=[cx,cy,i,j ];
%     end
% end


[X,Y] = meshgrid(1:length(x),1:length(y)); tt=[X(:) ,Y(:)] ;
t1 =x(tt(:,1) )' ;
t2 =y(tt(:,2) )'  ;
SQcenter=[t1,t2,tt] ;




% HClattice.HCcenter=HCcenter;
% HClattice.nrow=nrow;
% HClattice.nColumn=nColumn;

SQMapping= SQcenter ;
% toc
end

