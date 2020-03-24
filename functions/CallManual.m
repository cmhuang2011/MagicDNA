function  CallManual(UpperPB_H,evn,hyperB,AssFigH,patchH,SelectEdge,Pop_NumberConn,lineIndex )
%UNTITLED2 manually assign bridging connection between two bundles
%   Detailed explanation goes here
% dfgdg=3
MainScreeSize = get(0,'screensize');
fLocal= figure('name','Manually Assign','Position',MainScreeSize,'numbertitle','off');
fLocal.Units='normalized';fLocal.OuterPosition=[0 0 1 1];  clf;
drawnow;warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFig = get(handle(fLocal), 'JavaFrame'); jFig.setMaximized(true); drawnow;


clf; hold on;  axis equal;grid on;
ax=gca;  ax.Position(3)=0.65; ax.Position(1)=0.1;
ax.Parent=fLocal;
str2=strcat('Number of Connection =  ',{' '} ,num2str(Pop_NumberConn.Value ));
set(gcf,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) )
title(str2);  ax.UserData.OriTitle= str2 ;

CellSaveScatter=cell(2,1);
new_handle_Scat=cell(2,1);
new_handle_patch=cell(2,1) ;
for twoBundlei=1:2
    XYZ=[];
    ScatH= patchH{3,SelectEdge(twoBundlei)};
    XYZ(:,3)= ScatH.ZData;
    XYZ(:,1)= ScatH.XData;
    XYZ(:,2)= ScatH.YData;
    
    CellSaveScatter{twoBundlei}=XYZ;
    %    scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'MarkerFaceColor',ScatH.CData);
    CopyScatterh=patchH{3,SelectEdge(twoBundlei)};   %Copy patch 3-scatter
    new_handle_Scat{twoBundlei} = copyobj(CopyScatterh,ax);
    %      new_handle_Scat{twoBundlei}.SizeData
    %     new_handle_Scat{twoBundlei}.SizeData=50;
    %    new_handle_Scat{twoBundlei}.UserData=twoBundlei;
    new_handle_Scat{twoBundlei}.UserData.OriScatColor=mean(new_handle_Scat{twoBundlei}.CData,1);
    new_handle_Scat{twoBundlei}.MarkerFaceColor=new_handle_Scat{twoBundlei}.CData ;
    new_handle_Scat{twoBundlei}.HitTest='on';
    new_handle_Scat{twoBundlei}.UserData.SelectEdge=SelectEdge;
    new_handle_Scat{twoBundlei}.UserData.ij=twoBundlei;
    
    fLocal.UserData.ScatH(twoBundlei)=new_handle_Scat{twoBundlei};
    
    CopyVolh=patchH{1,SelectEdge(twoBundlei)};     %Copy patch 1-Volume
    new_handle_patch{twoBundlei} = copyobj(CopyVolh,ax);
    new_handle_patch{twoBundlei}.FaceAlpha=0.1;
    new_handle_patch{twoBundlei}.PickableParts ='none';
    if twoBundlei==1
        new_handle_patch{twoBundlei}.FaceColor= [0.4,1,0.4 ] ;
        new_handle_Scat{twoBundlei}.CData= ones(length(new_handle_Scat{twoBundlei}.XData) ,1)*[0,0.8,0.2 ] ;
        new_handle_Scat{twoBundlei}.MarkerFaceColor= [0,0.7,0 ] ;
        %        new_handle_Scat{twoBundlei}.MarkerEdgeColor= [0,0.8,0.2  ] ;
        new_handle_Scat{twoBundlei}.UserData.OriScatColor=[0,0.8,0.2 ] ;
    else
        new_handle_patch{twoBundlei}.FaceColor= [0.4,0.4,1 ] ;
        new_handle_Scat{twoBundlei}.CData= ones(length(new_handle_Scat{twoBundlei}.XData) ,1)*[0,0.2,0.8 ];
        new_handle_Scat{twoBundlei}.MarkerFaceColor= [0,0,0.7 ] ;
        %        new_handle_Scat{twoBundlei}.MarkerEdgeColor= [0,0.2,0.8 ] ;
        new_handle_Scat{twoBundlei}.UserData.OriScatColor=[0,0.2,0.8 ] ;
    end
    new_handle_Scat{twoBundlei}.MarkerFaceColor='flat';
%     new_handle_Scat{twoBundlei}.MarkerFaceColor=[0.5 0.5 0.5] ;
    new_handle_Scat{twoBundlei}.MarkerEdgeColor='flat';
end



AllCombinationDistTable=zeros(size(CellSaveScatter{1},1)*size(CellSaveScatter{2},1),3);
for iPoint=1:size(CellSaveScatter{1},1)
    for jPoint=1:size(CellSaveScatter{2},1)
        dd=norm([CellSaveScatter{1}(iPoint,:)-CellSaveScatter{2}(jPoint,:)]);
        AllCombinationDistTable( (iPoint-1)*size(CellSaveScatter{2},1)+jPoint,:)=[dd,iPoint,jPoint];
    end
end

SortRes=sortrows(AllCombinationDistTable,1);
RefinBi=SortRes(1,2);
RefinBj=SortRes(1,3);
Pdata=[CellSaveScatter{1}(RefinBi,:);CellSaveScatter{2}(RefinBj,:)];

pRefH=plot3(Pdata(:,1),Pdata(:,2),Pdata(:,3),':','HitTest','off');
XYZ2(:,3)= new_handle_Scat{1}.ZData;
XYZ2(:,1)= new_handle_Scat{1}.XData;
XYZ2(:,2)= new_handle_Scat{1}.YData;

RefXYZ=XYZ2(RefinBi,:);

dM=XYZ2-ones(size(XYZ2,1),1)* RefXYZ  ;
dd1=zeros(size(dM,1),1);
for k=1:length(dd1)
    dd1(k)=norm( dM(k,:));
end

sld1 = uicontrol('Parent',fLocal,'Style', 'slider','Units','normalized',...
    'Position', [0.8 0.88 0.15 0.03]);
sld1.Callback=@(src,evn)movesld1(src,evn,new_handle_Scat{1},RefinBi,dd1) ;
sld1.Value=1; sld1.Max= length(new_handle_Scat{1}.XData); sld1.Min=1;

XYZ3(:,1)= new_handle_Scat{2}.XData;
XYZ3(:,2)= new_handle_Scat{2}.YData;
XYZ3(:,3)= new_handle_Scat{2}.ZData;

RefXYZ=XYZ3(RefinBj,:);
%    scatter3(RefXYZ(1),RefXYZ(2),RefXYZ(3),'s')

dM=  XYZ3 -ones(size(XYZ3,1),1)*RefXYZ  ;
dd2=zeros(size(dM,1),1);
for k=1:length(dd2)
    dd2(k)=norm( dM(k,:));
end

sld2 = uicontrol('Parent',fLocal,'Style', 'slider','Units','normalized',...
    'Position', [0.8 0.73 0.15 0.03]);
sld2.Callback=@(src,evn)movesld2(src,evn,new_handle_Scat{2},RefinBj,dd2) ;
sld2.Value=1; sld2.Max= length(new_handle_Scat{2}.XData); sld2.Min=1;

NRequiredConn=Pop_NumberConn.Value;
fLocal.UserData.XoverInd=zeros(NRequiredConn,2);

set(fLocal,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) ) ;
% set(fLocal,'KeyPressFcn',{@(src,evn)keySelect(src,evn)}) ;


PopAssign = uicontrol('Parent',fLocal,'Units','normalized','style','popup','position',[0.8 0.58 0.1 0.1],...
    'string',num2cell(1:NRequiredConn));
sld3 = uicontrol('Parent',fLocal,'Style', 'slider','Units','normalized','Value',0.1, 'Min',0.1, 'Max',1,...
    'Position', [0.8 0.55 0.15 0.03]);
sld3.Callback=@(src,evn)movesld3(src,evn,new_handle_patch) ;


btn_Yes = uicontrol('Parent',fLocal,'Units','normalized',....
    'Position',[0.8 0.4 0.08 0.1],...
    'String','OK', 'Callback',{@(src,evn)ChooseThisPair(src,evn,NRequiredConn,PopAssign) });
btn_Yes.UserData.plotLinkH=cell(1,Pop_NumberConn.Value);

btn_Close = uicontrol('Parent',fLocal,'Units','normalized',....
    'Position',[0.9 0.4 0.08 0.1],...
    'String','Close', 'Callback',{@(src,evn)closeMA(src,evn,fLocal.UserData.XoverInd,UpperPB_H,hyperB,lineIndex ) });
% fLocal.CloseRequestFcn=@(src,evn)showNotFinished(src,evn,fLocal.UserData.XoverInd,UpperPB_H) ;

movesld1(sld1,[],new_handle_Scat{1},[],dd1) ;
movesld2(sld2,[],new_handle_Scat{2},[],dd2) ;
ax.View=patchH{3,SelectEdge(twoBundlei)}.Parent.View ;
str1=strcat('For Bundle', num2str(SelectEdge(1)), ', slide to show more' ) ;
txtH1 = uicontrol(fLocal,'Style','text','Unit','normalized','FontSize',14,'Position', [0.8 0.93 0.2 0.03], 'String',str1,'HorizontalAlignment','left');
str2=strcat('For Bundle', num2str(SelectEdge(2)),', slide to show more' ) ;
txtH2 = uicontrol(fLocal,'Style','text','Unit','normalized','FontSize',14,'Position', [0.8 0.78 0.2 0.03], 'String',str2,'HorizontalAlignment','left');
xlabel('X') ;ylabel('Y') ;zlabel('Z') ;
%--------------------
ax2 = axes('Position',[0.8 0.05 0.15 0.15],'Tag','AxForBundlei'); hold on ; % show base position in projective plane 02282020
Bundlei=hyperB.containBundle{SelectEdge(1)} ;
CylCenters= Bundlei.CylInplanePosition ;r=1;
theta=[(0:1:359)' ;nan]; np=length(theta) ;
dxy= r*[cosd(theta),sind(theta) ];
AllPoints= zeros( size(CylCenters,1)*np ,2) ;
pt_h=cell(size(CylCenters,1),1) ;
for k=1:size(CylCenters,1)
    AllPoints( np*(k-1)+1:np*k , :) = ones(np,1)*CylCenters(k,:) + dxy ;
    QQ=ones(np,1)*CylCenters(k,:) + dxy ;
    if ismember( k,Bundlei.AGroup)
        pt_h{k} = patch(QQ(1:end-1,1),QQ(1:end-1,2) ,[0.8,0.8,0] ,'FaceAlpha',0.2) ;
    else
        pt_h{k} = patch(QQ(1:end-1,1),QQ(1:end-1,2) ,[0.2,0.7,0.7] ,'FaceAlpha',0.2) ;
    end
    pt_h{k}.UserData.Ind= k ;
    pt_h{k}.UserData.State= 0 ;
    
end
plot(AllPoints(:,1),AllPoints(:,2) ,'k','LineWidth',2 ,'HitTest','off'  ) ;
text(CylCenters(:,1),CylCenters(:,2) ,num2str((1:k)') ,'FontSize',16 ,'HorizontalAlignment','center','HitTest','off'  ,'clipping', 'on') ;
axis equal ;axis ij ;   axis off;


%-----------------------
ax3 = axes('Position',[0.8 0.22 0.15 0.15],'Tag','AxForBundlej'); hold on ;  % ax2.UserData.NoZ=true ;
Bundlei=hyperB.containBundle{SelectEdge(2)} ;
CylCenters= Bundlei.CylInplanePosition ;r=1;
theta=[(0:1:359)' ;nan]; np=length(theta) ;
dxy= r*[cosd(theta),sind(theta) ];
AllPoints= zeros( size(CylCenters,1)*np ,2) ;
pt_h=cell(size(CylCenters,1),1) ;
for k=1:size(CylCenters,1)
    AllPoints( np*(k-1)+1:np*k , :) = ones(np,1)*CylCenters(k,:) + dxy ;
    QQ=ones(np,1)*CylCenters(k,:) + dxy ;
    if ismember( k,Bundlei.AGroup)
        pt_h{k} = patch(QQ(1:end-1,1),QQ(1:end-1,2) ,[0.8,0.8,0] ,'FaceAlpha',0.2) ;
    else
        pt_h{k} = patch(QQ(1:end-1,1),QQ(1:end-1,2) ,[0.2,0.7,0.7] ,'FaceAlpha',0.2) ;
    end
    pt_h{k}.UserData.Ind= k ;
    pt_h{k}.UserData.State= 0 ;
    
end
plot(AllPoints(:,1),AllPoints(:,2) ,'k','LineWidth',2 ,'HitTest','off'  ) ;
text(CylCenters(:,1),CylCenters(:,2) ,num2str((1:k)') ,'FontSize',16 ,'HorizontalAlignment','center','HitTest','off' ,'clipping', 'on') ;
axis equal ;axis ij ;  axis off;
%---------------------------
new_handle_Scat{1}.ButtonDownFcn=@(src,evn)ClickScatter(src,evn,hyperB,ax2,ax);
new_handle_Scat{2}.ButtonDownFcn=@(src,evn)ClickScatter(src,evn,hyperB,ax3,ax);

axes(ax); uistack(pRefH,'bottom') ;
AssignIcon( btn_Yes,'ManualConnect.jpg' );
AssignIcon( btn_Close,'FinishManual.jpg' );

uiwait ;
end

% function showNotFinished(src,evn,MAForceList,UpperPB_H,hyperB)
%     if nnz(MAForceList)==2*size(MAForceList,1)
%     L1fig=UpperPB_H.Parent.Parent;
%
%     UpperPB_H.UserData.MAResult=MAForceList;
%     % closereq;
%      delete(gcf)
%      delete(L1fig)
%     else
%      warndlg(' Not Finished  ','!! Warning !!')
%     end
% end

function closeMA(src,evn,MAForceList,UpperPB_H,hyperB,lineIndex)
if nnz(MAForceList)==2*size(MAForceList,1)
    L1fig=UpperPB_H.Parent.Parent;
    % UpperPB_H.UserData.MAResult=MAForceList;
    % closereq;
    hyperB.choice.Bundle1=4;
    hyperB.choice.Bundle2=4;
    hyperB.choice.NumberOfConnect=size(MAForceList,1);
    hyperB.choice.UseMethod='Manual';
    hyperB.choice.ManualRes=MAForceList;
    hyperB.TakeSeveralV3{lineIndex}=MAForceList;
    delete(gcf)
    delete(L1fig)
end
end




function ChooseThisPair(src,evn,NRequiredConn,PopAssign)
fH=src.Parent;ax=gca;
[~,XoverIndin1]=ismember([1,0,0],fH.UserData.ScatH(1).CData,'rows');
pCoor1=[fH.UserData.ScatH(1).XData(XoverIndin1),...
    fH.UserData.ScatH(1).YData(XoverIndin1),....
    fH.UserData.ScatH(1).ZData(XoverIndin1)];

[~,XoverIndin2]=ismember([1,0,0],fH.UserData.ScatH(2).CData,'rows');

pCoor2=[fH.UserData.ScatH(2).XData(XoverIndin2),...
    fH.UserData.ScatH(2).YData(XoverIndin2),....
    fH.UserData.ScatH(2).ZData(XoverIndin2)];

Coor=[pCoor1;pCoor2] ;



if ~ismember(XoverIndin1,fH.UserData.XoverInd(:,1)) && ~ismember(XoverIndin2,fH.UserData.XoverInd(:,2))
    
    if ~isempty(src.UserData.plotLinkH{PopAssign.Value})
        delete(src.UserData.plotLinkH{PopAssign.Value})
    end
    
    fH.UserData.XoverInd(PopAssign.Value,:)=[XoverIndin1,XoverIndin2];
    src.UserData.plotLinkH{PopAssign.Value}= plot3(Coor(:,1), Coor(:,2), Coor(:,3),'k');
    % elseif ismember(XoverIndin1,fH.UserData.XoverInd(:,1)) && ismember(XoverIndin2,fH.UserData.XoverInd(:,2))
    
else
    warndlg('Node Repeated','!! Warning !!')
end

nC=0;
for iConn=1:NRequiredConn
    if ~isempty(src.UserData.plotLinkH{iConn})
        nC=nC+1;
    end
end
nMore= NRequiredConn-nC;
str2=strcat('Required Connection =   ',{' '} ,num2str(nMore ));
title(str2);

xlabel('X') ;ylabel('Y') ;zlabel('Z') ;
% HelloPB=1;

end


% function keySelect(src,evn)
% right=0;
% switch    evn.Character
%
%     case 'q'    %bundle 1 scatter move up
%         right=1;
%         scattH=src.UserData.ScatH(1);
%
%
%         Avail=find(src.UserData.ScatH(1).SizeData==150);
%         [~,Select]=ismember([1,0,0],src.UserData.ScatH(1).CData,'rows');
%         if isempty(Select)||Select==0; Select=  Avail(randi(length(Avail)));  end
%         if Select~=Avail(1)
%             Become=Avail(find(Avail==Select)-1);
%             if  ~ismember(Select,Avail);Become=Avail(randi(length(Avail)));end
%         else
%             Become=Avail(end);
%         end
%
%     case 'a'    %bundle 1 scatter move down
%         right=1;
%         scattH=src.UserData.ScatH(1);
%         Avail=find(src.UserData.ScatH(1).SizeData==150);
%         [~,Select]=ismember([1,0,0],src.UserData.ScatH(1).CData,'rows');
%         if isempty(Select)||Select==0 ; Select=  Avail(randi(length(Avail)));  end
%         if Select~=Avail(end)
%             Become=Avail(find(Avail==Select)+1);
%             if  ~ismember(Select,Avail);Become=Avail(randi(length(Avail)));end
%         else
%             Become=Avail(1);
%         end
%
%     case 'e'   %bundle 2 scatter move up
%         right=1;
%         scattH=src.UserData.ScatH(2);
%
%         Avail=find(scattH.SizeData==150);
%         [~,Select]=ismember([1,0,0],scattH.CData,'rows');
%         if isempty(Select)||Select==0; Select=  Avail(randi(length(Avail)));  end
%         if Select~=Avail(1)
%             Become=Avail(find(Avail==Select)-1);
%             if  ~ismember(Select,Avail);Become=Avail(randi(length(Avail)));end
%         else
%             Become=Avail(end);
%         end
%
%
%     case 'd'    %bundle 2 scatter move down
%         right=1;
%         scattH=src.UserData.ScatH(2);
%
%         Avail=find(scattH.SizeData==150);
%         [~,Select]=ismember([1,0,0],scattH.CData,'rows');
%         if isempty(Select)||Select==0; Select=  Avail(randi(length(Avail)));  end
%         if Select~=Avail(end)
%             Become=Avail(find(Avail==Select)+1);
%             if  ~ismember(Select,Avail);Become=Avail(randi(length(Avail)));end
%         else
%             Become=Avail(1);
%         end
%
% end
% if right==1
%     scattH.CData=scattH.UserData.OriScatColor;
%     scattH.CData=repmat(scattH.CData,size( scattH.XData'));
%     scattH.CData(Become,:)=[1,0,0];
%
%
% end
% end


function ClickScatter(src,evn,hyperB,axInplane,axMain)
%
XYZ(:,3)= src.ZData;
XYZ(:,1)= src.XData;
XYZ(:,2)= src.YData;
distances = sqrt(sum(bsxfun(@minus, XYZ, evn.IntersectionPoint).^2,2));
ind=distances==min(distances);
src.CData=src.UserData.OriScatColor;
src.CData=repmat(src.CData,size( src.XData'));
src.CData(ind,:)=[1,0,0];
% src.CData

src.UserData;

BundleInd =  src.UserData.SelectEdge( src.UserData.ij);
%  BundleInd =  src.UserData.SelectEdge( src.UserData.ij);
ExternalXoverAsFB = hyperB.containBundle{BundleInd}.ExternalXoverAsFB ;
NodeInfo = ExternalXoverAsFB(ind,:)  ;
%  XYZ (ind,:)
%  src.Pa
if  NodeInfo(1)>0  % side nodes
    CylinderInd = NodeInfo(1) ;
    NicksBase= NodeInfo([2 4]) ;
    strAdd=  sprintf('Bundle: %i, Cyl: %i, Nick:[%i %i]',BundleInd, CylinderInd,NicksBase )  ;
    %     StrAdd = strcat('Cyl-' , num2str(CylinderInd) , '  ','Nick-[',NicksBase
    title([ src.Parent.UserData.OriTitle strAdd] ) ;
else
    CylindersInd = -NodeInfo([1 3] ) ;
    EndBase =  NodeInfo([2 4]) ;
    strAdd=  sprintf('Bundle: %i, Cyl: [%i %i], EndBases:[%i %i]',BundleInd, CylindersInd,EndBase )  ;
    title([ src.Parent.UserData.OriTitle strAdd] ) ;
end
% sdfsdf=3;

Position3D = [src.XData(ind) ; src.YData(ind) ;src.ZData(ind) ] ;
TM = hyperB.containBundle{BundleInd}.TransformMatrix2     ;

% A=TM(1:3,1:3) ;
% NodeInfo
% inv(TM(1:3,1:3))
% XYZinplane = inv(TM(1:3,1:3))*( NodeInfo(5:7)'-TM(1:3,4) );
XYZinplane = inv(TM(1:3,1:3))*( Position3D-TM(1:3,4) );

XY = ( XYZinplane(1:2) )  ;

if isempty(findobj(axInplane,'type','scatter'))
    sH =scatter(axInplane,XY(1),  XY(2),86, 'r','filled' ) ;
else
    sH=findobj(axInplane,'type','scatter')  ;
    sH.XData =XY(1);  sH.YData =XY(2);
end


if isequal( findobj(gcf,'Tag','AxForBundlei'), axInplane )
    sH.CData = [ 0,1, 0] ;
else
    sH.CData = [ 0,0, 1] ;
end

%  mean( hyperB.containBundle{BundleInd}.CylInplanePosition )

axes(axMain) ;




end

function movesld1(src,evn,ScatterH,RefinBi,dd)
src.Value=round(src.Value);

[A,B]=sort(dd);
scatsize=zeros(1,length(A));
scatsize(  B(1:src.Value) )=36;
scatsize( B(src.Value+1:end) )=5;
ScatterH.SizeData=scatsize';

end

function movesld2(src,evn,ScatterH,RefinBj,dd)
src.Value=round(src.Value);
[A,B]=sort(dd);
scatsize=zeros(1,length(A));
Big=B(1:src.Value) ;
Small=setdiff(1:length(dd),Big);

scatsize(Big )=36;
scatsize( Small)=5;
ScatterH.SizeData=scatsize';

end

function movesld3(src,evn,new_handle_patch)

new_handle_patch{1}.FaceAlpha=src.Value ;
new_handle_patch{2}.FaceAlpha=src.Value ;

end
