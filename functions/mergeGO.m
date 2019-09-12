function outputHandle = mergeGO( pHcell , targetAxe )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%---allocate
n=0;
for k=1:length(pHcell)
   n=n+ length(pHcell{k}.XData) ;
   n=n+1 ; 
end
XYZ =zeros(n,3) ;
n2=1;
for k=1:length(pHcell)
XYZ(n2:n2+length(pHcell{k}.XData)-1,:  ) = [pHcell{k}.XData' , pHcell{k}.YData' , pHcell{k}.ZData'] ;
XYZ(n2+length(pHcell{k}.XData),:)     = [NaN ,NaN,NaN] ;
n2=n2+length(pHcell{k}.XData)+1 ;
end
 XYZ(end,:)=[];
outputHandle=plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3)) ;

outputHandle.LineStyle =pHcell{1}.LineStyle ; 
outputHandle.LineWidth =pHcell{1}.LineWidth ;
outputHandle.Color= pHcell{1}.Color ;

outputHandle.Parent = targetAxe ;
% sdfsf=3
end

