function [ HCMapping ] = findHClatticeMapping( r ,XYbound)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% tic

% x=-XYbound(1):r*sqrt(3):XYbound(1);
% x= 0:r*sqrt(3):XYbound(1) ;
x=sort(unique( [ 0:r*sqrt(3):XYbound(1),  0:-r*sqrt(3):-XYbound(1)] )) ; 



% y=-XYbound(2):3*r:XYbound(2);
y=sort(unique( [ 0:3*r:XYbound(2),  0:-3*r:-XYbound(2)] ) ) ;

% clc
% tic

% nColumn=length(x);
% nrow=length(y);
% HCcenter=zeros(nColumn*nrow,4);
% for i=1:length(x)
%     for j=1:length(y)
%       
%      if mod(i+j,2)==0
%          cy=y(j)+0.5*r;
%      else
%          cy=y(j)-0.5*r;
%      end
%      cx=x(i);
% %     scatter(cx,cy,'.');
%     k=(i-1)*nrow+j;
%     HCcenter(k,:)=[cx,cy,i,j ];
%     end
% end

% toc

[X,Y] = meshgrid(1:length(x),1:length(y)); tt=[X(:) ,Y(:)] ;
ind2D=mod(tt(:,1) +tt(:,2) ,2)==0;
t1 =x(tt(:,1) )' ;
t2 =y(tt(:,2) )' + 0.5*ind2D -0.5*(~ind2D) ;
HCcenter=[t1,t2,tt] ;





% HClattice.HCcenter=HCcenter;
% HClattice.nrow=nrow;
% HClattice.nColumn=nColumn;

HCMapping= HCcenter ;
% toc
end

