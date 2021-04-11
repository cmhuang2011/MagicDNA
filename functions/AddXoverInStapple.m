function [whichstrands, NStrands,NoOfBreakStrand] = AddXoverInStapple(StappStrands,StapXover )
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here

CrossCyls=union(StapXover(:,1),StapXover(:,1));
ZPosition=mean(union(StapXover(:,2),StapXover(:,2)));
whichstrands=[]    ;          %find out Xover between which two strands
CorrLoation=[];
% whichstrands2=[]    ;          %find out Xover between which two strands
% CorrLoation2=[];


NStrands=StappStrands;
for cellindex=1:size(StappStrands,1)
   TemMat=StappStrands{cellindex};
%    for ScanInMat=1:size(TemMat,1)-1
%        LastP=TemMat(ScanInMat,:);
%        NextP=TemMat(ScanInMat+1,:)   ;    
%        if LastP(1)==NextP(1) && ismember(LastP(1),CrossCyls) && (LastP(2)-ZPosition)*(NextP(2)-ZPosition)<0          
%            whichstrands=[whichstrands cellindex]    ;
%            CorrLoation=[CorrLoation ;ScanInMat];
%        end      
%    end   
   
   S1 = TemMat(1:end-1,1) ==TemMat(2:end,1)  ;
   S2 = ismember( TemMat(1:end-1,1),CrossCyls)   ;
   S3 =  ( TemMat(1:end-1,2)-ZPosition).*(TemMat(2:end,2)-ZPosition)<0         ;
   
   if sum(and(and(S1,S2),S3) )>0
%        sdfsf=3 ;
       whichstrands=[whichstrands cellindex]    ;
       CorrLoation=[CorrLoation; find(and(and(S1,S2),S3))];
   end
       
 
end
% whichstrands
% whichstrands2
% CorrLoation
% CorrLoation2
% if sum(whichstrands~=whichstrands2)>0  || sum(CorrLoation~=CorrLoation2)>0 
% sdsf=2
% end

NoOfBreakStrand=length(intersect(whichstrands,whichstrands));

switch NoOfBreakStrand
    case 0
        NStrands=StappStrands   ; 
        return
    case 1
      
        StrandT=StappStrands{whichstrands(1)}(1:CorrLoation(1),:);
         if length(CorrLoation)==1
             NStrands=StappStrands;
             return           
         end
         
         
        StrandB=StappStrands{whichstrands(1)}(CorrLoation(2)+1:end,:);       
        RemainStrand=StappStrands{whichstrands(1)}(CorrLoation(1)+1:CorrLoation(2),:);
        StrandT=[StrandT ;StrandT(end,1) ZPosition];
        StrandB=[StrandB(1,1) ZPosition ;StrandB];
        RemainStrand=[ RemainStrand(1,1) ZPosition; RemainStrand;  RemainStrand(end,1) ZPosition];
        
        if StrandT(end,2)<StrandT(end-1,2)
            StrandT(:,2)=ceil(StrandT(:,2));
            StrandB(:,2)=ceil(StrandB(:,2));
            NStrands{whichstrands(1)}=[StrandT ;StrandB];
            RemainStrand(:,2)=floor(RemainStrand(:,2));          
            NStrands{size(StappStrands,1)+1}=RemainStrand;
        else
            StrandT(:,2)=floor(StrandT(:,2));
            StrandB(:,2)=floor(StrandB(:,2));
            NStrands{whichstrands(1)}=[StrandT; StrandB];
            RemainStrand(:,2)=ceil(RemainStrand(:,2));            
            NStrands{size(StappStrands,1)+1}=RemainStrand;
        end
        
   
        
    case 2

        NewStrand1T=StappStrands{whichstrands(1)}(1:CorrLoation(1),:);
        NewStrand2T=StappStrands{whichstrands(2)}(1:CorrLoation(2),:);
        NewStrand1B=StappStrands{whichstrands(1)}(CorrLoation(1)+1:end,:);
        NewStrand2B=StappStrands{whichstrands(2)}(CorrLoation(2)+1:end,:);

        NewStrand1T=[NewStrand1T ; NewStrand1T(end,1) ZPosition];
        NewStrand2T=[NewStrand2T ; NewStrand2T(end,1) ZPosition];
        NewStrand1B=[NewStrand1B(1,1) ZPosition;NewStrand1B];
        NewStrand2B=[NewStrand2B(1,1) ZPosition;NewStrand2B];

        if  (NewStrand1T(end,2)-NewStrand1T(end-1,2))*(NewStrand2T(end,2)-NewStrand2T(end-1,2))>0
           NNS1= [NewStrand1T; flip(NewStrand2T)];
           NNS2=[flip(NewStrand1B); NewStrand2B]   ;
        else
             NNS1= [NewStrand1T;  NewStrand2B];
             NNS2=[NewStrand2T; NewStrand1B];
        end
     zzzz=NNS1(find(NNS1(:,2)==ZPosition)-1,2);
        if    zzzz(1)<ZPosition
            NNS1(:,2)=floor(NNS1(:,2));
            NNS2(:,2)=ceil(NNS2(:,2));
        else
            NNS1(:,2)=ceil(NNS1(:,2));
            NNS2(:,2)=floor(NNS2(:,2));
        end

        NStrands{whichstrands(1)}=NNS1;
        NStrands{whichstrands(2)}=NNS2;
end

% a=3333

end

