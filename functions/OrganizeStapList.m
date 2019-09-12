function [ Output ] = OrganizeStapList( StappList,Part ,tt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
StapBreakPoint=Part.stapBP;
IJTable=zeros(1000,3);
n=0;
OStappL=StappList;
for i=1:length(StappList)
    MasterStrand=StappList{i};
    MasterSP=MasterStrand(1,:);
    MasterEP=MasterStrand(end,:);
   for j=i+1:length(StappList)
     SlaveStrand=StappList{j};
     SlaveSP=SlaveStrand(1,:);
     SlaveEP=SlaveStrand(end,:);    
      
    Mat=[MasterSP ;MasterEP;SlaveSP;SlaveEP];
    Trigger=0;
    if abs(Mat(1,2)-Mat(3,2))<=1  && Trigger==0
         MasIndex=  intersect(find(Mat(1,1)==StapBreakPoint(:,1)), find(Mat(1,2)==StapBreakPoint(:,2)) );
         SlaIndex=  intersect(find(Mat(3,1)==StapBreakPoint(:,1)), find(Mat(3,2)==StapBreakPoint(:,2)) );
         if ~isempty(SlaIndex)&& ~isempty(MasIndex)  && ceil(MasIndex/4)==ceil(SlaIndex/4)
             n=n+1;
             IJTable(n,1:2)=[i j];
             IJTable(n,3)=1;
             Trigger=1;
         end
    end
    if abs(Mat(1,2)-Mat(4,2))<=1 && Trigger==0
         MasIndex=  intersect(find(Mat(1,1)==StapBreakPoint(:,1)), find(Mat(1,2)==StapBreakPoint(:,2)) );
         SlaIndex=   intersect(find(Mat(4,1)==StapBreakPoint(:,1)), find(Mat(4,2)==StapBreakPoint(:,2)) );
         if ~isempty(SlaIndex)&& ~isempty(MasIndex) && ceil(MasIndex/4)==ceil(SlaIndex/4)
             n=n+1;
             IJTable(n,1:2)=[i j];
             IJTable(n,3)=2;
             Trigger=1;
         end
    end
    if abs(Mat(2,2)-Mat(3,2))<=1 && Trigger==0
         MasIndex=  intersect(find(Mat(2,1)==StapBreakPoint(:,1)), find(Mat(2,2)==StapBreakPoint(:,2)) );
         SlaIndex=  intersect(find(Mat(3,1)==StapBreakPoint(:,1)), find(Mat(3,2)==StapBreakPoint(:,2)) );

         if ~isempty(SlaIndex)&& ~isempty(MasIndex)  && ceil(MasIndex/4)==ceil(SlaIndex/4)
           
              n=n+1;
              IJTable(n,1:2)=[i j];
               IJTable(n,3)=3;
               Trigger=1;
         end
    end
    
    if abs(Mat(2,2)-Mat(4,2))<=1 && Trigger==0
         MasIndex=  intersect(find(Mat(2,1)==StapBreakPoint(:,1)), find(Mat(2,2)==StapBreakPoint(:,2)) );
         SlaIndex=   intersect(find(Mat(4,1)==StapBreakPoint(:,1)), find(Mat(4,2)==StapBreakPoint(:,2)) );
         if ~isempty(SlaIndex) && ~isempty(MasIndex)  && ceil(MasIndex/4)==ceil(SlaIndex/4)
           
            n=n+1;
            IJTable(n,1:2)=[i j];
            IJTable(n,3)=4;
            Trigger=1;
         end
   
     end
   
   end
end
IJTable(find(IJTable(:,1)==0),:)=[];
AppearStrands=union(IJTable(:,1:2),[]);
SSS=IJTable(:,1:2);
RepeatTable=zeros(length(AppearStrands),2);
RepeatTable(:,1)=AppearStrands;
for i=1:length(SSS(:))
RepeatTable(find(RepeatTable(:,1)==SSS(i)),2)=RepeatTable(find(RepeatTable(:,1)==SSS(i)),2)+1;
end
RemoveCell=cell(size(IJTable,1),1);
for t=1:size(IJTable,1)
    RemoveCell{t}=union(IJTable(t,1),IJTable(t,2));
end

List=1:size(RemoveCell,1);
for y=1:size(RemoveCell,1)
    Elements=RemoveCell{y};
    for x=1:length(Elements)
        DList=setdiff(List,y);
        for z=1:length(DList)
           if  ismember(Elements(x),RemoveCell{ DList(z)}  )
               RemoveCell{y}=union(RemoveCell{y}, RemoveCell{ DList(z)} );
           end            
        end        
    end    
end
% 
if nnz(IJTable)==0   %means no strand should be connected
Output=StappList;
    return   
end
[CateIJTable, CaseCell]=CateElements(IJTable(:,1:3));
% DoingSeq=size(IJTable,1):-1:1

UltimateStrand=[];
for s=1:size(CateIJTable,1)
    SingleDoingList=CateIJTable{s};
    CorrespondingInstruct=CaseCell{s};
    for w=1:length(CorrespondingInstruct)
        CaseFlag=IJTable(CorrespondingInstruct(w),3);
%         Seg1=OStappL{SingleDoingList(w)}
%         Seg2=OStappL{SingleDoingList(w+1)}
        Seg1=OStappL{  IJTable(CorrespondingInstruct(w),1) };
        Seg2=OStappL{  IJTable(CorrespondingInstruct(w),2)  };
        
        switch CaseFlag
            case 1 
            NewStrand=[flip(Seg2,1) ;Seg1]    ;    
            case 2
            NewStrand=[  Seg2   ;Seg1];
            case 3
            NewStrand=[  Seg1   ;Seg2]   ;
            case 4
            NewStrand=[  Seg1   ;flip(Seg2,1)];

        end
        OStappL{IJTable(CorrespondingInstruct(w),1)}=NewStrand;
        OStappL{IJTable(CorrespondingInstruct(w),2)}=NewStrand;
    end
  UltimateStrand=union(UltimateStrand,SingleDoingList(w+1));
  
end
DeleteStrands=setdiff(AppearStrands,UltimateStrand);
for h=1:length(DeleteStrands)
 OStappL{   DeleteStrands(h)}=[];
end
Output=OStappL(~cellfun('isempty',OStappL))  ;
end

