function [Output,CaseCell ] = CateElements( InputMatrix )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

A=cell(max(max(InputMatrix)),1);
A{1}=union(A{1},InputMatrix(1,1:2));
UsedCell=1;
SelectedElements=union([],InputMatrix(1,1:2));
for i=2:size(InputMatrix,1)
     E1=   InputMatrix(i,1);
     E2=   InputMatrix(i,2);
    %  if i==5
    %      asdfsaf=122414
    %  end 
     if ismember(E1,SelectedElements) ||  ismember(E2,SelectedElements)
         Gindex=[];
          for j=1:UsedCell
            if ismember(E1,A{j}) || ismember(E2,A{j})
            Gindex=union(j,Gindex);
            end
          end
    %       A{:}
    %       [E1 E2]
    %       Gindex
    %       i
        if length(Gindex)==1
          A{Gindex}=union(A{Gindex},[E1 E2]);
          SelectedElements=union(SelectedElements,[E1 E2]);

        else
%             Gindex
%             A{Gindex(1)}
%             A{Gindex(2)}
            A{Gindex(1)}=union(A{Gindex(1)},A{Gindex(2)} );
            A{Gindex(1)}=union(A{Gindex(1)},[E1 E2]);
            A{Gindex(2)}=[];
            SelectedElements=union(SelectedElements,[E1 E2]);

        end
     else

         UsedCell=UsedCell+1;
         UsedCell ;
         A ;
         A{UsedCell}=union(A{UsedCell},[E1 E2]);
         SelectedElements=union(SelectedElements,[E1 E2]);  
     end

end

UsedIndex=1:size(InputMatrix,1);
Output=A(~cellfun('isempty',A))  ;

for k=1:size(Output,1)
   ArrayinCell= Output{k};
   Count=zeros(size(ArrayinCell));
   for  CC=1: size(InputMatrix,1)
       if ismember(InputMatrix(CC,1),ArrayinCell)
       Count( find(ArrayinCell==InputMatrix(CC,1)))= Count( find(ArrayinCell==InputMatrix(CC,1)))+1;
       end
       
       if ismember(InputMatrix(CC,2),ArrayinCell)
       Count( find(ArrayinCell==InputMatrix(CC,2)))= Count( find(ArrayinCell==InputMatrix(CC,2)))+1;
       end
       
   end  
   
   if length(ArrayinCell)~=2
       NArr=zeros(size(ArrayinCell));
       First=find(Count==1);
       if isempty(First)
       First=find(Count==2);
       First=First(1);
       else
       First=First(1)  ;    
       end
       NArr(1)=ArrayinCell(First);       
       Second=union( find(InputMatrix(:,1)==NArr(1)) , find(InputMatrix(:,2)==NArr(1)));
        
        Second=Second(1);
        NArr(2)=setdiff(InputMatrix(Second,1:2), NArr(1));
       UsedIndex=setdiff(UsedIndex,Second);
       
       for q=3:length(ArrayinCell)
       TargetIndex=union( find(InputMatrix(:,1)==NArr(q-1)) , find(InputMatrix(:,2)==NArr(q-1)));
       TargetIndex=intersect(TargetIndex,UsedIndex);
%        IM=InputMatrix(TargetIndex,1:2);
%        Na=NArr(q-1);
%        ArrayinCell;
%        Count;
%        TargetIndex;
%        NArr(q-1);
%        InputMatrix;
       NArr(q)=setdiff(InputMatrix(TargetIndex,1:2), NArr(q-1));
       UsedIndex=setdiff(UsedIndex,TargetIndex)    ;        
       end
   Output{k}=NArr ;      
   end
end

 CaseCell=cell(size(Output));
for j2=1:size(CaseCell,1)
    Vec=Output{j2};
    InstrucVec=zeros(size(Vec,1)-1,1);
    for k2=1:size(Vec,1)-1
        Edge=[Vec(k2) Vec(k2+1)];
        AAA=0;
        for Scan=1:size(InputMatrix,1)
            Tar=[InputMatrix(Scan,1) InputMatrix(Scan,2)];
            if length(union(Edge,Tar))==2
                AAA=Scan ;               
            end
            
            
        end
        
        if AAA~=0
%         InstrucVec(k2)=InputMatrix(AAA,3)
         InstrucVec(k2)=AAA;
        end
    end
   CaseCell{j2}= InstrucVec;
    
end
 
end

