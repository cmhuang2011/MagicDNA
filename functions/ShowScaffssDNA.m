function  ShowScaffssDNA( GetHyperB )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fH= gcf ;
QQ=GetHyperB.ScafUnUsed ;

if ~isempty(QQ) % if have saved overhangs sequences , July 3 2019
    dQQ =diff(QQ)~=1;
    Head=[QQ(1),QQ(find(dQQ)+1 )] ;
    Tails=   [QQ(find(dQQ) ) , QQ(end)];
    HeadTail=[Head;Tails] ;
    UnUsedSeq= cell(size(HeadTail,2) ,1) ;
    for k=1: length(UnUsedSeq)
        UnUsedSeq{k}=GetHyperB.pSeqAll(HeadTail(1,k):HeadTail(2,k)) ;
    end
    %     GetHyperB.ScafUnUsedSeq = UnUsedSeq ;
end


f96 = figure(96); clf ; hold on ;
axCurr= gca;
% [ScafHelixXYZ,ScafBasesCenter  ]=plotScaf2_Helix_V2_noGraphics( GetHyperB,{GetHyperB.scafC5},Isstap ,[0,0,1] ,TM ) ;     % get scaf strands coordinate

% ScafssStrands = cell(1, size(HeadTail,2)) ;  % C5
% cc=1;
% for scafi = 1: length(GetHyperB.ScafAllBase)
% for k = 1:length(ScafssStrands)
%     %    if HeadTail(1,k)~= HeadTail(2,k)
%     ScafssStrands{cc} = GetHyperB.ScafAllBase{scafi}(HeadTail(1,k):HeadTail(2,k),: ) ;
%     %    end
%     cc=cc+1 ;
% end
% end
Isstap= 0 ;  TM=1 ;
[ScafHelixXYZ,ScafBasesCenter  ]=plotScaf2_Helix_V2_noGraphics( GetHyperB,GetHyperB.scafC5,Isstap ,[0,0,1] ,TM ) ;     % get scaf strands coordinate

GatherAll= [];
for k=1: length(ScafHelixXYZ)
    GatherAll=[GatherAll; ScafHelixXYZ{k}];
end

pH=cell(size(HeadTail,2),1) ; tH= cell(size(HeadTail,2),1) ;
for k = 1: size(HeadTail,2)
    IndStartEnd =  HeadTail(:,k);
    XYZ= GatherAll(IndStartEnd(1):IndStartEnd(2) , :) ;
    
    pH{k} = plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3) ,'-ok','MarkerFaceColor','k' ,'MarkerSize',3) ;
    mXYZ= mean(XYZ ,1) ;
    tH{k} = text(mXYZ(1), mXYZ(2), mXYZ(3)  ,  strcat('\leftarrow', num2str(abs(diff(IndStartEnd) )+1  )),'Fontsize',12  ,'HitTest', 'off', 'clipping', 'on') ;
end

for k = 1: size(HeadTail,2)
pH{k}.ButtonDownFcn=@(src,evn)HighLightSelect(src,evn,pH,tH) ;
    
end

AssemblyMain=findobj(fH,'Tag','AssemblyMain');
h_patch = findobj(AssemblyMain,'Type','Patch');
new_PatchH = copyobj(h_patch,gca);
for patchi=1:length(new_PatchH)
    new_PatchH(patchi).FaceAlpha=0.2;
    new_PatchH(patchi).HitTest='off';
    new_PatchH(patchi).PickableParts='none' ;
end


%
%     new_handle_patch{iBun}.FaceAlpha=0.1;
%     new_handle_patch{iBun}.PickableParts ='none';


axis equal




% sdsf=3

set(f96,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) ) 
xlabel('X') ; ylabel('Y') ; zlabel('Z') ; %box on
set(axCurr,'View',AssemblyMain.View);



end

function HighLightSelect(src,evn,pH,tH) 

for k =1: length(pH)
 if  isequal(pH{k},src )
     Select = k  ;
 end
 pH{k}.Color= [0,0,0] ; tH{k}.Color= [0,0,0] ;
end

 pH{Select}.Color= [1,0,0] ; tH{Select}.Color= [1,0,0] ;

% sdfsf=3

end


