function [ Res ] = showstapins( NewStapList )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

InitialLength=zeros(size(NewStapList,1),1);
IniSegLength=cell(size(NewStapList,1),1);
TT=[];


for i=1:length(InitialLength)
    
    Strand=NewStapList{i};
    LL=0;
    SegL=zeros(size(Strand,2)-1,1);
        for j=1:size(Strand,1)-1        
           if mod(j,2)==1
               Add=1+abs(Strand(j,2)-Strand(j+1,2));
           else
              Add=0 ;        
           end
           SegL(j)=Add;
           LL=LL+Add;      
        end    
    InitialLength(i)=LL;
    SegL(SegL==0)=[];
    IniSegLength{i}=SegL;
    
    
end


Res.InitialLength=InitialLength;
Res.IniSegLength=IniSegLength;

end



