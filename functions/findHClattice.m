function [ HClattice ] = findHClattice( r ,XYbound)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


x=0:r*sqrt(3):XYbound(1);
nColumn=length(x);


y=0:3*r:XYbound(2);
nrow=length(y);

HCcenter=zeros(nColumn*nrow,2);

for i=1:length(x)
    for j=1:length(y)
        
     if mod(i+j,2)==0
         cy=y(j)+0.5*r;
     else
         cy=y(j)-0.5*r;
     end
     cx=x(i);
%     scatter(cx,cy,'.');
    
    k=(i-1)*nrow+j;
    HCcenter(k,:)=[cx,cy];
    

    end
end

HClattice.HCcenter=HCcenter;
HClattice.nrow=nrow;
HClattice.nColumn=nColumn;


end

