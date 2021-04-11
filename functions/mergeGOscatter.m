function outputHandle = mergeGOscatter( pHcell , targetAxe )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%---allocate
n=0;
for k=1:length(pHcell)
   n=n+ length(pHcell{k}.XData) ;
   n=n+1 ; 
end
XYZ =zeros(n,3) ; 
if isfield(pHcell{1}.UserData,'Seq')
seqallocate =  char(ones(1,n )) ;
end

n2=1;
for k=1:length(pHcell)
    XYZ(n2:n2+length(pHcell{k}.XData)-1,:  ) = [pHcell{k}.XData' , pHcell{k}.YData' , pHcell{k}.ZData'] ;
    XYZ(n2+length(pHcell{k}.XData),:)     = [NaN ,NaN,NaN] ;
    if isfield(pHcell{1}.UserData,'Seq')
        seqallocate(n2:n2+length(pHcell{k}.XData)-1) = pHcell{k}.UserData.Seq ;
        seqallocate(n2+length(pHcell{k}.XData)) = '?' ;
    end
    n2=n2+length(pHcell{k}.XData)+1 ;

end
XYZ(end,:)=[];
    if isfield(pHcell{1}.UserData,'Seq')
    seqallocate(end)='';
    end
outputHandle=scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3)) ;

% outputHandle.CData= pHcell{1}.CData ;
% outputHandle.Marker = pHcell{1}.Marker ;
    if isfield(pHcell{1}.UserData,'Seq')
    outputHandle.UserData.Seq= seqallocate;
    end

outputHandle.Parent = targetAxe ;
% sdfsf=3
end

