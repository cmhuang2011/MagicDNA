function [ StapleCell ] = CalibStapDir( StapleCell,RVec)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here


for j=1:size(StapleCell,1)
    
    StappStrand2=  StapleCell{j};
    FirstCyl=StappStrand2(1,1);
   
    
%     if  FirstCyl>28
%         sdfsf=234
%         
%     end
    NANOCylnum=RVec(FirstCyl);
    
     if mod(NANOCylnum,2)==1 && StappStrand2(2,2)>StappStrand2(1,2)                 
     StappStrand=StappStrand2;
     elseif mod(NANOCylnum,2)==0 && StappStrand2(2,2)<StappStrand2(1,2)   
     StappStrand=StappStrand2;
     elseif mod(NANOCylnum,2)==0 && StappStrand2(2,2)>StappStrand2(1,2)   
     StappStrand=flip(StappStrand2);
     elseif mod(NANOCylnum,2)==1 && StappStrand2(2,2)<StappStrand2(1,2)  
     StappStrand=flip(StappStrand2);
     else
        StappStrand=StappStrand2;  %--------
     end    
    
%     show=StapleCell{j} ;
     
     
    StapleCell{j}=StappStrand;
    
end


end

