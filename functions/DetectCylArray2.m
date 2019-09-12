function [ OutP ,in,on,Nxy] = DetectCylArray2( OuterVertexs) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Vertexs=OuterVertexs(1:end-1,:);
d=2 ; %nm


Boundary=[min(Vertexs(:,1)) max(Vertexs(:,1)) min(Vertexs(:,2)) max(Vertexs(:,2)) ] ;
% nx=[flip(0:-d:Boundary(1)-d) d:d:Boundary(2)+d];
% ny=[flip(0:-d:Boundary(3)-d) d:d:Boundary(4)+d];
nx=[floor(Boundary(1)/d)*d :d: ceil(Boundary(2)/d)*d];
ny=[floor(Boundary(3)/d)*d :d: ceil(Boundary(4)/d)*d];

k=0;
TestPointsInG=zeros(length(nx)*length(ny),3);
for i=1:length(nx)
    for j=1:length(ny)
        k=k+1;
        TestPointsInG(k,:)=[ nx(i) ny(j) 0];        
    end
end

OutP=TestPointsInG(:,1:2);
[in,on] = inpolygon(OutP(:,1)',OutP(:,2)',OuterVertexs(:,1)',OuterVertexs(:,2)');
OutP=OutP';
Nxy=[length(nx) length(ny)];



end

