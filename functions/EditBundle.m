function EditBundle(src,evn,GetHyperB,fH,patchH,popupH )
%UNTITLED Summary of this function goes here
%   for editting cylinder models of assigned bundle
clc

% decision.pair=2;
AskWhichBundle(src,evn,GetHyperB,fH,patchH ,popupH.Value(1))  ;

SelectBundle= src.UserData.SelectBundle ;
TotalBundle = length(GetHyperB.containBundle ) ;

% MainScreeSize = get(0,'screensize');% MainScreeSize(3)=round(0.5*MainScreeSize(3)) ;
% fLocal= figure(124) ;
% % fLocal.name='Edit Bundle dSDNA' ;
% fLocal.Position=MainScreeSize ;  fLocal.numbertitle='off';
fLocal= figure('name','Edit Bundle dSDNA','numbertitle','on' );
fLocal.Units='normalized';fLocal.OuterPosition=[0 0 1 1];  clf;

drawnow;warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFig = get(handle(fLocal), 'JavaFrame'); jFig.setMaximized(true); drawnow;


fLocal.UserData.InvShift =1 ;




clf; hold on;  axis equal;grid on;
% ax1 = axes('Position',[0.1 0.1 0.7 0.8]);
ax1=gca;  ax1.Position(3)=0.55; ax1.Position(1)=0.05;
ax1.Parent=fLocal;
str2=strcat('For Bundle ',{' '} ,num2str(SelectBundle ));
title(str2)

newLines =  cell(TotalBundle,1);
new_handle_patch=  cell(TotalBundle,1);
% newNodes =  cell(TotalBundle,1);
for iBun=1:TotalBundle
    
    CopyLineh=patchH{2,iBun};     %Copy patch 1-Volume
    newLines{iBun} = copyobj(CopyLineh,ax1);
    %      newNodes{iBun} =cell(1,length( newLines{iBun})) ;
    %----------------------
    CopyVolh=patchH{1,iBun};     %Copy patch 1-Volume
    new_handle_patch{iBun} = copyobj(CopyVolh,ax1);
    new_handle_patch{iBun}.FaceAlpha=0.1;
    new_handle_patch{iBun}.PickableParts ='none';
    %------------------
    if iBun==SelectBundle
        new_handle_patch{iBun}.FaceColor= [0.8,0.2,0.2 ] ;
        %         newLines{iBun}.Color =
        SaveNodes = zeros(2*length( newLines{iBun}) ,3) ;
        for k=1:length( newLines{iBun})
            newLines{iBun}(k).Color =0.2*ones(1,3)  ;
            newLines{iBun}(k).LineWidth=4 ;
            newLines{iBun}(k).HitTest='off' ;
            SaveNodes(2*k-1:2*k,:) = [newLines{iBun}(k).XData ;newLines{iBun}(k).YData ;newLines{iBun}(k).ZData ]';
            newLines{iBun}(k).UserData=[] ;
            newLines{iBun}(k).UserData.OriData = [newLines{iBun}(k).XData ; newLines{iBun}(k).YData ; newLines{iBun}(k).ZData ] ;
        end
        %           sdfs=3
        newNodes_3D= scatter3(SaveNodes(:,1),SaveNodes(:,2),SaveNodes(:,3) ,128,'sb','filled' ) ;
        newNodes_3D.UserData.OriData=SaveNodes ;
    else
        
        new_handle_patch{iBun}.FaceColor= [0,0.2,0.8 ] ;
        for k=1:length( newLines{iBun})
            newLines{iBun}(k).Color =0.8*ones(1,3)  ;
            newLines{iBun}(k).LineWidth=2 ;
            newLines{iBun}(k).HitTest='off' ;
        end
    end
    % -----------------------------
end
newNodes_3D.CData= repmat(newNodes_3D.CData ,length(newNodes_3D.XData),1) ;


ax1.View=patchH{1,1}.Parent.View ;
ax1.XLim=patchH{1,1}.Parent.XLim ;
ax1.YLim=patchH{1,1}.Parent.YLim ;
ax1.ZLim=patchH{1,1}.Parent.ZLim ;

xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
%------------------2D projection

ax2 = axes('Position',[0.63 0.65 0.35 0.28]); hold on ; ax2.UserData.NoZ=true ;
Bundle =GetHyperB.containBundle{SelectBundle}  ;
PlanarNodes = zeros(2*length(Bundle.Zbase1) ,2) ;
h_2D_lines = cell(length(Bundle.Zbase1) ,1) ;

for k=1: length(Bundle.Zbase1)
    PlanarNodes(2*k-1:2*k,:) =[ [Bundle.Zbase1(k),Bundle.Zbase2(k)]' , [k;k] ] ;
    h_2D_lines{k}=plot([Bundle.Zbase1(k),Bundle.Zbase2(k)], [k,k] ,'HitTest','off' ,'LineWidth',2  ) ;
    h_2D_lines{k}.UserData.OriData = [Bundle.Zbase1(k),Bundle.Zbase2(k); [k,k] ] ;
end
h_2DNodes = scatter(PlanarNodes(:,1), PlanarNodes(:,2) ,'ob','filled') ; h_2DNodes.SizeData=64;
h_2DNodes.CData= repmat(h_2DNodes.CData ,size(PlanarNodes,1),1) ;
h_2DNodes.ButtonDownFcn=@(src,evn)HitScatter(src,evn,newNodes_3D) ;
h_2DNodes.UserData.OriData=PlanarNodes ;

newNodes_3D.ButtonDownFcn=@(src,evn)HitScatter2(src,evn,h_2DNodes) ;

%--------change to paired colors
PairListForThisBundle = GetHyperB.SavePremPair.CellPairList{SelectBundle} ;
DefColors = get(groot,'defaultAxesColorOrder') ;
for k=1:size(PairListForThisBundle,1)
    if k<= size(DefColors,1)
        h_2D_lines{PairListForThisBundle(k,1) }.Color = DefColors(k,:) ;
        h_2D_lines{PairListForThisBundle(k,2) }.Color = DefColors(k,:) ;
        %    hello=1
    else
        QQ =mod(k,size(DefColors,1)) ;
        if QQ==0
            QQ=size(DefColors,1);
        end
        OtherColor =DefColors(QQ ,:) + 0*(rand(1,3)-0.5) ;
        %      OtherColor=OtherColor/max(OtherColor)
        h_2D_lines{PairListForThisBundle(k,1) }.Color = OtherColor  ;
        h_2D_lines{PairListForThisBundle(k,2) }.Color = OtherColor ;        
    end
end


%---------


%  plot([Bundle.Zbase1()   )

CXLim=ax2.XLim ;
% ax2.XLim=[0, CXLim(2) ] ;
ylabel('Cylinders in this bundle')
xlabel('Cylinder length (bp)')
axis ij ; grid on ;

btn_UpdateBundle = uicontrol(fLocal,'Style', 'pushbutton', 'String', 'Update','Unit','normalized', 'Position', [0.83 0.33 0.08 0.1] );
btn_UpdateBundle.Callback=@(src,evn)UpdateBundle(src,evn,GetHyperB,fH,patchH,h_2DNodes,SelectBundle,newLines{SelectBundle} ) ;

btn_SaveBundle = uicontrol(fLocal,'Style', 'pushbutton', 'String', 'Save this bundle','Unit','normalized', 'Position', [0.88 0.33 0.08 0.1] );
btn_SaveBundle.Callback=@(src,evn)SaveBundle(src,evn,GetHyperB,fH,patchH,h_2DNodes,SelectBundle,newLines{SelectBundle} ) ;



btn_SelectLeftAll = uicontrol(fLocal,'Style', 'pushbutton', 'String', 'SelectLeftAll','Unit','normalized', 'Position', [0.63 0.45 0.04 0.1] );
btn_SelectLeftAll.Callback=@(src,evn)SelectLeftAll(src,evn,h_2DNodes,newNodes_3D) ;

btn_SelectRightAll = uicontrol(fLocal,'Style', 'pushbutton', 'String', 'SelectRightAll','Unit','normalized', 'Position', [0.75 0.45 0.04 0.1] );
btn_SelectRightAll.Callback=@(src,evn)SelectRightAll(src,evn,h_2DNodes,newNodes_3D) ;

btn_CheckAll = uicontrol(fLocal,'Style', 'pushbutton', 'String', 'CheckAll','Unit','normalized', 'Position', [0.8 0.45 0.04 0.1] );
btn_CheckAll.Callback=@(src,evn)CheckAll(src,evn,h_2DNodes,newNodes_3D) ;

btn_UuCheckAll = uicontrol(fLocal,'Style', 'pushbutton', 'String', 'UuCheckAll','Unit','normalized', 'Position', [0.79 0.45 0.04 0.1] );
btn_UuCheckAll.Callback=@(src,evn)UnCheckAll(src,evn,h_2DNodes,newNodes_3D) ;

h_text = text(PlanarNodes(:,1), PlanarNodes(:,2) , strcat('\leftarrow',num2str(PlanarNodes(:,1)) ),'HitTest','off' ,'Clipping','on' ) ;
checkH_ShowNT = uicontrol(fLocal, 'Style', 'checkbox','String', 'Show/Hide','Unit','normalized','Position', [0.9 0.58 0.04 0.03]);
checkH_ShowNT.Callback=@(src,evn)showNT(src,evn,h_2DNodes,h_text)  ;

middle = 0.5*(PlanarNodes(1:2:end,:) +   PlanarNodes(2:2:end,:) ) ;
Diff_nm =  0.34* (PlanarNodes(2:2:end,1) -  PlanarNodes(1:2:end,1)) ;
h_lengthLabel=text(middle(:,1), middle(:,2) , strcat('\downarrow',num2str(Diff_nm,3),' nm' ),'HitTest','off' ,'VerticalAlignment','bottom','Clipping','on') ;
% for k=1:length(h_lengthLabel)
%     h_lengthLabel(k).
% end
% h_lengthLabel(1).UserData.htext =h_text;

title(ax2,'Step size =1. Use control to change ');



AssignIcon( btn_SelectLeftAll,'EditLeft.jpg' ) ; btn_SelectLeftAll.TooltipString='Select all nodes on the left side(Z1).' ;
AssignIcon( btn_SelectRightAll,'EditRight.jpg' ) ;btn_SelectRightAll.TooltipString='Select all nodes on the right side(Z2).' ;
AssignIcon( btn_CheckAll,'EditAll.jpg' ) ; btn_CheckAll.TooltipString='Select all nodes.' ; 
AssignIcon( btn_UuCheckAll,'EditNone.jpg' ) ; btn_UuCheckAll.TooltipString='Unselect all nodes.' ; 
% align([btn_SelectLeftAll btn_SelectRightAll btn_CheckAll btn_UuCheckAll],'distribute','bottom');
align([btn_SelectLeftAll btn_SelectRightAll btn_CheckAll btn_UuCheckAll btn_UpdateBundle btn_SaveBundle],'distribute','top');


AssignIcon( btn_UpdateBundle,'EditUpdate.jpg' ); btn_UpdateBundle.TooltipString='Temporary change the cylinder model. After editing all bundles, it needs to update the mechanism in Assebmly by exporting all bundle again. ' ;
AssignIcon( btn_SaveBundle,'EditSaveBundle.jpg' ); btn_SaveBundle.TooltipString='Save this cylinder model into the library for futre use as Bottom-up approach.' ;
% align([btn_UpdateBundle btn_SaveBundle ],'distribute','bottom');

ax2.YTick=1:length(Bundle.Zbase1) ;
ylim([0, length(Bundle.Zbase1)+1] ) ;

uistack(h_2DNodes,'top') ; % July 17 2019
set(ax2, 'XGrid', 'on', 'YGrid', 'off') ;
%-------------


ax3 = axes('Position',[0.63 0.08 0.35 0.35]); hold on ; ax3.UserData.NoZ=true ;
CylCenters= Bundle.CylInplanePosition ;
r=1;
theta=[(0:1:359)' ;nan]; np=length(theta) ;
dxy= r*[cosd(theta),sind(theta) ];
AllPoints= zeros( size(CylCenters,1)*np ,2) ;
pt_h=cell(size(CylCenters,1),1) ;
for k=1:size(CylCenters,1)
    AllPoints( np*(k-1)+1:np*k , :) = ones(np,1)*CylCenters(k,:) + dxy ;
    QQ=ones(np,1)*CylCenters(k,:) + dxy ;
    if ismember( k,Bundle.AGroup)
        pt_h{k} = patch(QQ(1:end-1,1),QQ(1:end-1,2) ,[0.8,0.8,0] ,'FaceAlpha',0.2) ;
    else
        pt_h{k} = patch(QQ(1:end-1,1),QQ(1:end-1,2) ,[0.2,0.7,0.7] ,'FaceAlpha',0.2) ;
    end
    pt_h{k}.UserData.Ind= k ;
    pt_h{k}.UserData.State= 0 ;
    
end
% AllPoints
plot(AllPoints(:,1),AllPoints(:,2) ,'k','LineWidth',2 ,'HitTest','off' ) ;
text(CylCenters(:,1),CylCenters(:,2) ,num2str((1:k)') ,'FontSize',16 ,'HorizontalAlignment','center','HitTest','off' ,'Visible','on') ;

MarkDistance1 = plot(AllPoints(1,1),AllPoints(1,2),'Visible','off' ,'LineWidth',2.5,'Color','m') ;
MarkDistance2 = text(AllPoints(1,1),AllPoints(1,2),'','Visible','off','FontSize',18,'Color','m') ;

PairListForThisBundle;
SwapAGroup = ismember(PairListForThisBundle(:,1) , Bundle.AGroup ) ; 
PairListForThisBundle(SwapAGroup,:) = flip(PairListForThisBundle(SwapAGroup,:) ,2) ;

plotPairingInax3 = cell(size(PairListForThisBundle,1) ,1) ;
% SaveXY = nan*ones(3*size(PairListForThisBundle,1)-1 ,2) ;
for k= 1: size(PairListForThisBundle,1)
% plotPairingInax3{k} = plot()
   SaveXY(1,:) =  [ mean(pt_h{PairListForThisBundle(k,1)}.Vertices(:,1)),mean(pt_h{PairListForThisBundle(k,1)}.Vertices(:,2))] ;
   SaveXY(2 ,:) =  [ mean(pt_h{PairListForThisBundle(k,2)}.Vertices(:,1)),mean(pt_h{PairListForThisBundle(k,2)}.Vertices(:,2))] ;
   plotPairingInax3{k} = plot(SaveXY(:,1) , SaveXY(:,2),'Color','b' ,'LineWidth',3 , 'Visible','off') ; 
    plotPairingInax3{k}.UserData.Ind= k ;
    plotPairingInax3{k}.UserData.State= 0 ;
end
ax3.UserData.OriPair= PairListForThisBundle ;
ax3.UserData.CurrentPair= PairListForThisBundle ;   

% plotPairingInax3=plot(SaveXY(:,1) , SaveXY(:,2),'Color','r' ,'LineWidth',2) ; 
for k= 1: length(PairListForThisBundle)
plotPairingInax3{k}.ButtonDownFcn= @(src,evn)ChangePair(src,evn,plotPairingInax3) ;
end
%         MarkDistance1.XData = [ mean( pt_h{States==1}.Vertices(:,1) ) ,  mean( pt_h{States==2}.Vertices(:,1)) ] ;
%         MarkDistance1.YData = [ mean( pt_h{States==1}.Vertices(:,2)) ,  mean( pt_h{States==2}.Vertices(:,2)) ] ;


for k=1:size(CylCenters,1)
    uistack(pt_h{k},'top')
    pt_h{k}.ButtonDownFcn= @(src,evn)clickPatch(src,evn, pt_h,MarkDistance1,MarkDistance2) ;
end
drawnow; 
for k= 1: size(PairListForThisBundle ,1)
 uistack(plotPairingInax3{k},'top')
end

try
axis equal ;
catch
end
axis ij ;
axis off;


axes(ax1) ;

%                Global_XYZ=Bundle.HelixRegardlessCylinder(1.06,  0 ,InplaneXY,BasesArr,ThefSimilarDirCylder) ;   %  already assign as scaffold domain (rr,isstap,......)
EndsPosition = zeros( length(Bundle.Zbase1) , 6  );   % [x y z ]@Z1 [x y z]@Ze
EndsOrientation = zeros( length(Bundle.Zbase1) , 6  );   % [x y z ]@Z1 [x y z]@Ze
for k= 1 :size(  EndsPosition ,1)
    rr= 0 ;
    XYZ = Bundle.HelixRegardlessCylinder(rr,  0 ,Bundle.CylInplanePosition(k,:) ,[Bundle.Zbase1(k),Bundle.Zbase2(k)],k) ;
    EndsPosition(k, :) = [  XYZ(1,1:3) ,XYZ(2,1:3)] ;
    rr= 1 ;
    XYZ2 = Bundle.HelixRegardlessCylinder(rr,  0 ,Bundle.CylInplanePosition(k,:) ,[Bundle.Zbase1(k),Bundle.Zbase2(k)],k) ;    
    EndsOrientation(k, :) = [  XYZ2(1,1:3)-XYZ(1,1:3) , XYZ2(2,1:3)-XYZ(2,1:3)] ;    
end
QuiverM =   [[EndsPosition(:,1:3); EndsPosition(:,4:6)] ,[EndsOrientation(:,1:3); EndsOrientation(:,4:6)]] ;


h_quiver = quiver3(QuiverM(:,1),QuiverM(:,2),QuiverM(:,3),QuiverM(:,4),QuiverM(:,5),QuiverM(:,6) ,'AutoScale','off','HitTest','off','Color',[0, 0.7, 0] )  ;
h_quiver.LineWidth =1.5 ;

set(gcf,'KeyPressFcn',@(src,evn)KeyPressMultiple(src,evn,ax1,ax2 ,newLines{SelectBundle},newNodes_3D,h_2D_lines,h_2DNodes ,h_text ,h_lengthLabel ,plotPairingInax3,ax3 ,GetHyperB,SelectBundle,h_quiver   ) )

% sdfsdf=3;
uistack(newNodes_3D,'top') ; % July 17 2019


%---------add legend as instructions 
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            ForLegend2=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            ForLegend3=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');    %ForLegend3=scatter3(mean(ax.XLim),mean(ax.YLim),mean(ax.ZLim));
            hLg= legend([ForLegend,ForLegend2,ForLegend3],'x','y','z' ,'Location','northwest' ) ;
            hLg.String={'[Q][A] =x';'[W][S] =y';'[E][D] =z'};
            hLg.Interpreter='tex';  %latex
            hLg.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            ForLegend2.Marker='.' ; ForLegend2.Marker='none';
            ForLegend3.Marker='.' ; ForLegend3.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_EditBundle( src,evn,ax1 );
            %             fHH.UserData.LgIndication.Title.Interpreter='latex';
            
            %             drawnow;
            hLg.Units='normalized'; hLg.AutoUpdate ='off';
            hLg.Title.String='Click me for instructions' ;
            hLg.Position=[0.0063 0.9528 0.1569 0.0387];
%---------------

popup3.TooltipString = 'Assign single-stranded lengths to the selected connections';




uiwait(fLocal);

end

function ChangePair(src,evn,plotPairingInax3 )
WhichLine =[];
 States = zeros(size(plotPairingInax3)) ;
for k= 1: length(plotPairingInax3)
      States(k) =  plotPairingInax3{k}.UserData.State ;
    if isequal(plotPairingInax3{k} ,src )
        WhichLine=k ;
    end
%     plotPairingInax3{k}.Color =    [0,0,1] ;
end
    if max(States) ==0
        States( src.UserData.Ind ) =1 ;
    elseif  max(States) ==1
        States(States==1) = 2 ;
        States( src.UserData.Ind ) =1 ;
    elseif  max(States) ==2
        if States( src.UserData.Ind )~=1
            States(States==2) = 0 ;
            States(States==1) = 2 ;
            States( src.UserData.Ind ) =1 ;
        end
    end

    for k=1:length(plotPairingInax3)
        plotPairingInax3{k}.UserData.State= States(k) ;
        if States(k)>0
             plotPairingInax3{k}.Color =    [1,0,0] ;
        else
             plotPairingInax3{k}.Color =    [0,0,1] ;
        end
    end

%     changeColor = 

States;
% WhichLine

% if evn.Button==3  % left click
%  plotPairingInax3{WhichLine}.Color =    [1,0,0] ;
%     
% end


ax3=gca;
% ax3.UserData.OriPair
% QQQQsdfssdfsgf=3
end


function  clickPatch(src,evn, pt_h,MarkDistance1,MarkDistance2)
    if evn.Button==1

    States = zeros(size(pt_h)) ;
    for k=1:length(pt_h)
        States(k) =  pt_h{k}.UserData.State ;
    end
    if max(States) ==0
        States( src.UserData.Ind ) =1 ;
    elseif  max(States) ==1
        States(States==1) = 2 ;
        States( src.UserData.Ind ) =1 ;
        if max(States)==1; return;end
        MarkDistance1.XData = [ mean( pt_h{States==1}.Vertices(:,1) ) ,  mean( pt_h{States==2}.Vertices(:,1)) ] ;
        MarkDistance1.YData = [ mean( pt_h{States==1}.Vertices(:,2)) ,  mean( pt_h{States==2}.Vertices(:,2)) ] ;
        MarkDistance2.Position=[ mean( MarkDistance1.XData)+0.2 ,  mean( MarkDistance1.YData) ]  ;
        MarkDistance2.String = strcat(num2str(sqrt(diff(MarkDistance1.XData ).^2 +  diff(MarkDistance1.YData ).^2 ),3),' nm' ) ;

        MarkDistance1.Visible='on';
        MarkDistance2.Visible='on';
          title(MarkDistance2.String) ;
    elseif  max(States) ==2
        if States( src.UserData.Ind )~=1
            States(States==2) = 0 ;
            States(States==1) = 2 ;
            States( src.UserData.Ind ) =1 ;
            MarkDistance1.XData = [ mean( pt_h{States==1}.Vertices(:,1) ) ,  mean( pt_h{States==2}.Vertices(:,1)) ] ;
            MarkDistance1.YData = [ mean( pt_h{States==1}.Vertices(:,2)) ,  mean( pt_h{States==2}.Vertices(:,2)) ] ;
            MarkDistance2.Position=[ mean( MarkDistance1.XData)+0.2 ,  mean( MarkDistance1.YData) ]  ;
            MarkDistance2.String = strcat(num2str(sqrt(diff(MarkDistance1.XData ).^2 +  diff(MarkDistance1.YData ).^2 ),3),' nm' ) ;       
            
            title(MarkDistance2.String) ;
        end 
    end
    for k=1:length(pt_h)
        pt_h{k}.UserData.State= States(k) ;
    end

    end

% States

% sdfsdf=3
end


function showNT(src,evn,h_2DNodes ,h_text)

if  src.Value==1
    for k=1: length(h_text)
        h_text(k).Visible='on' ;
    end
else
    for k=1: length(h_text)
        h_text(k).Visible='off' ;
    end
end

end

function UpdateText(h_2DNodes ,h_text ,h_lengthLabel)

for k=1: length(h_text)
    h_text(k).Position(1)= h_2DNodes.XData(k) ;
    h_text(k).String = strcat( '\leftarrow', num2str(h_2DNodes.XData(k))  )  ;
end

PlanarNodes= [h_2DNodes.XData' , h_2DNodes.YData' ] ;
middle = 0.5*(PlanarNodes(1:2:end,:) +   PlanarNodes(2:2:end,:) ) ;
Diff_nm =  0.34* (PlanarNodes(2:2:end,1) -  PlanarNodes(1:2:end,1)) ;
% h_lengthLabel=text(middle(:,1), middle(:,2) , strcat('\downarrow',num2str(Diff_nm),' nm' ),'HitTest','off' ,'VerticalAlignment','bottom') ;
for k = 1:length(h_lengthLabel)
    h_lengthLabel(k).Position = [middle(k,:)] ;
    h_lengthLabel(k).String =  strcat('\downarrow',num2str(Diff_nm(k),3),' nm' ) ;
end


% sdfsf=3 ;
end

function SaveBundle(src,evn,GetHyperB,fH,patchH,h_2DNodes,SelectBundle,new3DLine )



if  isa(GetHyperB.containBundle{SelectBundle} ,'BundleCylinderSQ')      %if SQ
    Part1=BundleCylinderSQ(1,[],h_2DNodes.XData(1:2:end) ,h_2DNodes.XData(2:2:end) ,GetHyperB.containBundle{SelectBundle}.CylInplanePosition)   ;
    Part1.TransformMatrix2 = GetHyperB.containBundle{SelectBundle}.TransformMatrix2 ;
elseif  isa(GetHyperB.containBundle{SelectBundle} ,'BundleCylinderHC')       %if HC
    Part1=BundleCylinderHC(1,[],h_2DNodes.XData(1:2:end) ,h_2DNodes.XData(2:2:end) ,GetHyperB.containBundle{SelectBundle}.CylInplanePosition)   ;
    Part1.TransformMatrix2 = GetHyperB.containBundle{SelectBundle}.TransformMatrix2 ;
end

InportHyperBundle=hyperbundle(1,Part1);  % first bundle
trial=0;
while trial<=1000
    [Doable,CellPairList,CellUnpair]=checkpairableH(InportHyperBundle,1);
    if  (Doable==1   ||  trial>=200)
        break;
    end
    trial=trial+1;
end
Doable ;
if  Doable==1
    Psave.GUISavePart=Part1;
    Psave.transM = GetHyperB.containBundle{SelectBundle}.TransformMatrix2 ;
    uisave({'Psave'},'Partx');
else
    opts = struct('WindowStyle','modal','Interpreter','tex');
    errordlg('\fontsize{18} Error! Cylinders can''t be paired. Consider changing postions of nodes or property: Tolerance in BundleCylinder class. Data is not saved','File Error',opts) ;
    
end

% sdf3=3

end



function UpdateBundle(src,evn,GetHyperB,fH,patchH,h_2DNodes,SelectBundle, new3DLine)


if  isa(GetHyperB.containBundle{SelectBundle} ,'BundleCylinderSQ')      %if SQ
    Part1=BundleCylinderSQ(1,[],h_2DNodes.XData(1:2:end) ,h_2DNodes.XData(2:2:end) ,GetHyperB.containBundle{SelectBundle}.CylInplanePosition)   ;
elseif  isa(GetHyperB.containBundle{SelectBundle} ,'BundleCylinderHC')       %if HC
    Part1=BundleCylinderHC(1,[],h_2DNodes.XData(1:2:end) ,h_2DNodes.XData(2:2:end) ,GetHyperB.containBundle{SelectBundle}.CylInplanePosition)   ;
end

InportHyperBundle=hyperbundle(1,Part1);  % first bundle
trial=0;
while trial<=1000
    [Doable,CellPairList,CellUnpair]=checkpairableH(InportHyperBundle,1);
    if  (Doable==1   ||  trial>=200)
        break;
    end
    trial=trial+1;
end
Doable ;
if  Doable==0
    opts = struct('WindowStyle','modal','Interpreter','tex');
%     errordlg('\fontsize{18} Cylinders can''t be paired. Consider changing postions of nodes or change propert:Tol in BundleCylinder class. Data is updated','File Error',opts) ;
    errordlg('\fontsize{18} Error! Cylinders can''t be paired. Consider changing postions of nodes or property: Tolerance in BundleCylinder class. Data is not saved','File Error',opts) ;
    
    return
end
%-----------




newZ1 =h_2DNodes.XData(1:2:end) ;
newZ2 =h_2DNodes.XData(2:2:end) ;

if strcmp(GetHyperB.containBundle{SelectBundle}.type, 'SQ') && min(newZ1)<32
    Shift =  -floor(min(newZ1)/32)*32  ;
elseif strcmp(GetHyperB.containBundle{SelectBundle}.type, 'HC') && min(newZ1)<21
    Shift =  -floor(min(newZ1)/21)*21  ;
else
    Shift=0 ;
end


GetHyperB.containBundle{SelectBundle}.Zbase1 =newZ1 + Shift ;
GetHyperB.containBundle{SelectBundle}.Zbase2 =newZ2 + Shift;

for k= 1:length(new3DLine)
    patchH{2,SelectBundle}(k).XData =new3DLine(k).XData  ;
    patchH{2,SelectBundle}(k).YData =new3DLine(k).YData  ;
    patchH{2,SelectBundle}(k).ZData =new3DLine(k).ZData  ;
end

close(src.Parent)

end



function HitScatter(src,evn ,newNodes_3D)
% d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2 ;
% Ind = d==min(d) ;
%  src.CData(Ind',:)=[1,0,0] ;
%  src.CData(~Ind',:)=repmat([0,0,1]  ,length(src.XData)-1,1) ;


d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2 ;
Ind = find(d==min(d)) ;
if evn.Button==1
    src.CData(Ind,:) =[1,0,0] ;
    newNodes_3D.CData(Ind,:) =[1,0,0] ;
elseif evn.Button==3
    src.CData(Ind,:) =[0,0,1] ;
    newNodes_3D.CData(Ind,:) =[0,0,1] ;
end
end

function SelectLeftAll(src,evn,h_2DNodes,newNodes_3D)
h_2DNodes.CData= repmat([1,0,0 ; 0, 0 1 ],length(h_2DNodes.XData )/2,1) ;
newNodes_3D.CData=repmat([1,0,0 ; 0, 0 1 ],length(h_2DNodes.XData )/2,1) ;
end
function SelectRightAll(src,evn,h_2DNodes,newNodes_3D)
h_2DNodes.CData= repmat([0,0,1 ; 1, 0, 0 ],length(h_2DNodes.XData )/2,1) ;
newNodes_3D.CData=repmat([0,0,1 ; 1, 0, 0 ],length(h_2DNodes.XData )/2,1) ;
end


function CheckAll(src,evn,h_2DNodes,newNodes_3D)
h_2DNodes.CData = ones(length(h_2DNodes.XData) ,1)*[1,0,0]  ;
newNodes_3D.CData=ones(length(h_2DNodes.XData), 1)*[1,0,0] ;
end

function UnCheckAll(src,evn,h_2DNodes,newNodes_3D)
h_2DNodes.CData = ones(length(h_2DNodes.XData) ,1)*[0,0,1]  ;
newNodes_3D.CData=ones(length(h_2DNodes.XData), 1)*[0,0,1] ;
end

function HitScatter2(src,evn ,h_2DNodes)
% d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2 ;
% Ind = d==min(d) ;
%  src.CData(Ind',:)=[1,0,0] ;
%  src.CData(~Ind',:)=repmat([0,0,1]  ,length(src.XData)-1,1) ;


d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2  + (src.ZData-evn.IntersectionPoint(3)).^2 ;
Ind = find(d==min(d))  ;
if evn.Button==1
    src.CData(Ind,:) =[1,0,0] ;
    h_2DNodes.CData(Ind,:) =[1,0,0] ;
elseif evn.Button==3
    src.CData(Ind,:) =[0,0,1] ;
    h_2DNodes.CData(Ind,:) =[0,0,1] ;
end

end


function AskWhichBundle(BtnEdit,evn,GetHyperB,fH,patchH,popValue)
% decision.pair=2;
decision.whichBundle=1;

fh = dialog('units','pixels','position',[300 300 300 150],...
    'menubar','none','name','Edit Bundle',...
    'numbertitle','off', 'resize','off');
movegui(fh,'center')
str=num2cell(1:length(GetHyperB.containBundle));

pop2 = uicontrol('style','pop','Parent',fh,...
    'units','pixels','position',[170 75 100 30],    'string',str);
pop2.Value=popValue ;

txtH2 = uicontrol('Style','text','FontSize',12,'Position', [140 110 150 20],....
    'String','Select Bundle');

btn = uicontrol('style','pushbutton','Parent',fh,...
    'unit','pix',    'position',[50 75 100 30],...
    'string','OK',  'callback',@(src,evn)export(src,evn,pop2,BtnEdit)    );
uiwait(fh);

% decision= pop2.UserData ;
% ChangeAxesLimit
end

function export(src,evn,pop2,BtnEdit)
% decision.pair=pop.Value;
BtnEdit.UserData.SelectBundle =pop2.Value;
decision.whichBundle=pop2.Value;
% decision

% pop2.UserData= pop2.Value;
close(gcf)
end

function KeyPressMultiple(src,evn ,ax1,ax2 ,newLines,newNodes,h_2D_lines,h_2DNodes,h_text ,h_lengthLabel ,plotPairingInax3,ax3 ,GetHyperB,SelectBundle,h_quiver )

% %     set(gcf,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) )
Bundle= GetHyperB.containBundle{SelectBundle} ;
ax=gca;
interval =2; R=1.1;
% Inv
evn.Character ;
evn.Key  ; % for arrows
switch evn.Key
    case 'control'
        src.UserData.InvShift=5 ;
        title(ax2,'Step size =5. Use shift to change ')
    case 'shift'
        src.UserData.InvShift=1 ;
        title(ax2,'Step size =1. Use control to change ')
        
    case 'leftarrow'
        Inv =src.UserData.InvShift;
        [~, NodeInd] = ismember(h_2DNodes.CData,[1,0,0],'rows') ;
        NodeInd=find(NodeInd) ;
        h_2DNodes.XData(NodeInd)  = h_2DNodes.XData(NodeInd) -Inv ;    % 2D scatter
        %-------------
        Odds = mod(NodeInd,2)==1 ;
        LinesInds = (NodeInd(Odds)+1)/2 ;
        for k=1: length(LinesInds)
            h_2D_lines{LinesInds(k)}.XData(1) =  h_2D_lines{LinesInds(k)}.XData(1) -Inv ;  % 2D Lines
        end
        Evvs = mod(NodeInd,2)==0 ;
        LinesIndsEvv = (NodeInd(Evvs))/2 ;
        for k=1: length(LinesIndsEvv)
            h_2D_lines{LinesIndsEvv(k)}.XData(2) =  h_2D_lines{LinesIndsEvv(k)}.XData(2) -Inv ; % 2D Lines
        end
        %--------------- 3D
        for k=1: length(LinesIndsEvv)
            Ref3D= newLines(LinesIndsEvv(k)).UserData.OriData ;  %[ X1 X2 ;Y1 Y2 ;Z1 Z2]
            Ref2D= h_2D_lines{LinesIndsEvv(k)}.UserData.OriData ;  %[ X1 X2 ;Y1 Y2 ;Z1 Z2]
            Current2D =    h_2D_lines{LinesIndsEvv(k)}.XData(2) ;
            
            vq = interp1(Ref2D(1,:)',Ref3D',Current2D,'linear' ,'extrap') ;
            newLines(LinesIndsEvv(k)).XData(2) = vq(1) ;
            newLines(LinesIndsEvv(k)).YData(2) = vq(2) ;
            newLines(LinesIndsEvv(k)).ZData(2) = vq(3) ;
        end
        
        for k=1: length(LinesInds)
            Ref3D= newLines(LinesInds(k)).UserData.OriData ;  %[ X1 X2 ;Y1 Y2 ;Z1 Z2]
            Ref2D= h_2D_lines{LinesInds(k)}.UserData.OriData ;  %[ X1 X2 ;Y1 Y2 ;Z1 Z2]
            Current2D =    h_2D_lines{LinesInds(k)}.XData(1) ;
            vq = interp1(Ref2D(1,:)',Ref3D',Current2D,'linear' ,'extrap') ;
            newLines(LinesInds(k)).XData(1) = vq(1) ;
            newLines(LinesInds(k)).YData(1) = vq(2) ;
            newLines(LinesInds(k)).ZData(1) = vq(3) ;
        end
        
        for k=1:size(newLines,1)
            Nodes = [ newLines(k).XData ;  newLines(k).YData ;  newLines(k).ZData] ;
            newNodes.XData(2*k-1:2*k) =  Nodes(1,1:2) ;
            newNodes.YData(2*k-1:2*k) =  Nodes(2,1:2) ;
            newNodes.ZData(2*k-1:2*k) =  Nodes(3,1:2) ;
        end
        %       axes(ax1) ;      axis equal ;
        UpdateText(h_2DNodes ,h_text ,h_lengthLabel) ;
        %-------update quiver, 07092019
        EndsPosition = zeros( length(Bundle.Zbase1) , 6  );   % [x y z ]@Z1 [x y z]@Ze
        EndsOrientation = zeros( length(Bundle.Zbase1) , 6  );   % [x y z ]@Z1 [x y z]@Ze
        for k= 1 :size(  EndsPosition ,1)
            rr= 0 ;
            XYZ = Bundle.HelixRegardlessCylinder(rr,  0 ,Bundle.CylInplanePosition(k,:) ,h_2DNodes.XData(2*k-1:2*k) ,k) ;
            EndsPosition(k, :) = [  XYZ(1,1:3) ,XYZ(2,1:3)] ;
            rr= 1 ;
            XYZ2 = Bundle.HelixRegardlessCylinder(rr,  0 ,Bundle.CylInplanePosition(k,:) ,h_2DNodes.XData(2*k-1:2*k),k) ;
            EndsOrientation(k, :) = [  XYZ2(1,1:3)-XYZ(1,1:3) , XYZ2(2,1:3)-XYZ(2,1:3)] ;
        end
        QuiverM =   [[EndsPosition(:,1:3); EndsPosition(:,4:6)] ,[EndsOrientation(:,1:3); EndsOrientation(:,4:6)]] ;
        h_quiver.XData = QuiverM(:,1) ;
        h_quiver.YData = QuiverM(:,2) ;
        h_quiver.ZData = QuiverM(:,3) ;
        h_quiver.UData = QuiverM(:,4) ;
        h_quiver.VData = QuiverM(:,5) ;
        h_quiver.WData = QuiverM(:,6) ;
        %--------
%         sdfsf=3
        
        
    case 'rightarrow'
        Inv =src.UserData.InvShift;
        [~, NodeInd] = ismember(h_2DNodes.CData,[1,0,0],'rows') ;
        NodeInd=find(NodeInd) ;
        h_2DNodes.XData(NodeInd)  = h_2DNodes.XData(NodeInd) +Inv ;      % 2D scatter
        
        Odds = mod(NodeInd,2)==1 ;
        LinesInds = (NodeInd(Odds)+1)/2 ;
        for k=1: length(LinesInds)
            h_2D_lines{LinesInds(k)}.XData(1) =  h_2D_lines{LinesInds(k)}.XData(1) +Inv ;  % 2D Lines
        end
        Evvs = mod(NodeInd,2)==0 ;
        LinesIndsEvv = (NodeInd(Evvs))/2 ;
        for k=1: length(LinesIndsEvv)
            h_2D_lines{LinesIndsEvv(k)}.XData(2) =  h_2D_lines{LinesIndsEvv(k)}.XData(2) +Inv ; % 2D Lines
        end
        %------------
        for k=1: length(LinesIndsEvv)
            Ref3D= newLines(LinesIndsEvv(k)).UserData.OriData ;  %[ X1 X2 ;Y1 Y2 ;Z1 Z2]
            Ref2D= h_2D_lines{LinesIndsEvv(k)}.UserData.OriData ;  %[ X1 X2 ;Y1 Y2 ;Z1 Z2]
            Current2D =    h_2D_lines{LinesIndsEvv(k)}.XData(2) ;
            vq = interp1(Ref2D(1,:)',Ref3D',Current2D,'linear' ,'extrap') ;
            newLines(LinesIndsEvv(k)).XData(2) = vq(1) ;
            newLines(LinesIndsEvv(k)).YData(2) = vq(2) ;
            newLines(LinesIndsEvv(k)).ZData(2) = vq(3) ;
        end
        
        for k=1: length(LinesInds)
            Ref3D= newLines(LinesInds(k)).UserData.OriData ;  %[ X1 X2 ;Y1 Y2 ;Z1 Z2]
            Ref2D= h_2D_lines{LinesInds(k)}.UserData.OriData ;  %[ X1 X2 ;Y1 Y2 ;Z1 Z2]
            Current2D =    h_2D_lines{LinesInds(k)}.XData(1) ;
            vq = interp1(Ref2D(1,:)',Ref3D',Current2D,'linear' ,'extrap') ;
            newLines(LinesInds(k)).XData(1) = vq(1) ;
            newLines(LinesInds(k)).YData(1) = vq(2) ;
            newLines(LinesInds(k)).ZData(1) = vq(3) ;
        end
        for k=1:size(newLines,1)
            Nodes = [ newLines(k).XData ;  newLines(k).YData ;  newLines(k).ZData] ;
            newNodes.XData(2*k-1:2*k) =  Nodes(1,1:2) ;
            newNodes.YData(2*k-1:2*k) =  Nodes(2,1:2) ;
            newNodes.ZData(2*k-1:2*k) =  Nodes(3,1:2) ;
        end
        UpdateText(h_2DNodes ,h_text ,h_lengthLabel) ;
        %       axes(ax1); axis auto ;
                %-------update quiver, 07092019
        EndsPosition = zeros( length(Bundle.Zbase1) , 6  );   % [x y z ]@Z1 [x y z]@Ze
        EndsOrientation = zeros( length(Bundle.Zbase1) , 6  );   % [x y z ]@Z1 [x y z]@Ze
        for k= 1 :size(  EndsPosition ,1)
            rr= 0 ;
            XYZ = Bundle.HelixRegardlessCylinder(rr,  0 ,Bundle.CylInplanePosition(k,:) ,h_2DNodes.XData(2*k-1:2*k) ,k) ;
            EndsPosition(k, :) = [  XYZ(1,1:3) ,XYZ(2,1:3)] ;
            rr= 1 ;
            XYZ2 = Bundle.HelixRegardlessCylinder(rr,  0 ,Bundle.CylInplanePosition(k,:) ,h_2DNodes.XData(2*k-1:2*k),k) ;
            EndsOrientation(k, :) = [  XYZ2(1,1:3)-XYZ(1,1:3) , XYZ2(2,1:3)-XYZ(2,1:3)] ;
        end
        QuiverM =   [[EndsPosition(:,1:3); EndsPosition(:,4:6)] ,[EndsOrientation(:,1:3); EndsOrientation(:,4:6)]] ;
        h_quiver.XData = QuiverM(:,1) ;
        h_quiver.YData = QuiverM(:,2) ;
        h_quiver.ZData = QuiverM(:,3) ;
        h_quiver.UData = QuiverM(:,4) ;
        h_quiver.VData = QuiverM(:,5) ;
        h_quiver.WData = QuiverM(:,6) ;
        %--------
    case 'return'  % for pairing 
        States = zeros(size(plotPairingInax3)) ;
        for k= 1: length(plotPairingInax3)
            States(k) =  plotPairingInax3{k}.UserData.State ;
        end
        
        IndLines = find(States>0 ) ;
        if length(IndLines)<2
            return
        end
        
        CollectXY=zeros(4,2) ; %[x1 x2 ;y1 y2 ;x3 x4 ;y3 y4 ]
        for k=1:length(IndLines)
            CollectXY(2*k-1 ,:) =   plotPairingInax3{IndLines(k)}.XData ;
            CollectXY(2*k ,:) =   plotPairingInax3{IndLines(k)}.YData ;
        end
        

        CylinderM = ax3.UserData.CurrentPair(IndLines,1:2) ;
        checkZ1 =abs(h_2DNodes.XData(2*CylinderM(1,1)-1)-h_2DNodes.XData(2*CylinderM(2,2)-1 )) ; % indexing = 2k-1
        checkZ2 =abs(h_2DNodes.XData(2*CylinderM(1,1))-h_2DNodes.XData(2*CylinderM(2,2) )) ;
        Check2 = 0 ;
        if checkZ1>15 || checkZ2>15   % also check the cylinders ends
            Check2 = 1 ;
%             return
        end
        
        QQ= ax3.UserData.CurrentPair ;
        QQ(IndLines(1) ,2) =ax3.UserData.CurrentPair(IndLines(2) ,2) ;
        QQ(IndLines(2) ,2) =ax3.UserData.CurrentPair(IndLines(1) ,2) ;
        ax3.UserData.CurrentPair =QQ  ;
        
        plotPairingInax3{IndLines(1)}.XData=[CollectXY(1,1) ,CollectXY(3,2)] ;
        plotPairingInax3{IndLines(1)}.YData=[CollectXY(2,1) ,CollectXY(4,2)] ;
        plotPairingInax3{IndLines(2)}.XData=[CollectXY(3,1) ,CollectXY(1,2)] ;
        plotPairingInax3{IndLines(2)}.YData=[CollectXY(4,1) ,CollectXY(2,2)] ;
        
        CurrentPair =ax3.UserData.CurrentPair ;
        
      
        CylInplane =GetHyperB.containBundle{SelectBundle}.CylInplanePosition ;
        Check = 1 ;
        for k = 1 :size(CurrentPair)
            d= CylInplane(CurrentPair(k,1) ,:) -CylInplane(CurrentPair(k,2) ,:) ;
            d2 = sqrt(d(1).^2+d(2).^2) ;
           checkZ1 =abs(h_2DNodes.XData(2*CurrentPair(k,1)-1)-h_2DNodes.XData(2*CurrentPair(k,2)-1 )) ; % indexing = 2k-1
            checkZ2 =abs(h_2DNodes.XData(2*CurrentPair(k,1))-h_2DNodes.XData(2*CurrentPair(k,2) )) ;
            
            if round(d2)~=2  ||  checkZ1>15 || checkZ2>15  % check the new pair still neighbor
                  Check = 0 ;
%                   k
            end       
                        
        end
        
        
        if Check ==1
         if isempty(GetHyperB.SavePremPairManual)
            GetHyperB.SavePremPairManual = GetHyperB.SavePremPair ;
            GetHyperB.SavePremPairManual.CellPairList{SelectBundle}  =ax3.UserData.CurrentPair ;
         else
            GetHyperB.SavePremPairManual.CellPairList{SelectBundle}  =ax3.UserData.CurrentPair ;
         end
        else
            fw=msgbox('Temporary change. Data not update. Check cylinder neighbors or ends  !! ');
         
        end
         
%         sdfsf=3
end

 ax = gca ;
 
switch  evn.Character
   
    case 'q'
        ax.XLim= ax.XLim + interval ;
    case 'Q'
        ax.XLim= R*(ax.XLim -mean(ax.XLim))+ mean(ax.XLim)   ;
    case 'a'
        ax.XLim= ax.XLim - interval ;
    case 'A'
        ax.XLim= (ax.XLim -mean(ax.XLim))/R + mean(ax.XLim)   ;
        %-----------
    case 'w'
        ax.YLim= ax.YLim + interval ;
    case 's'
        ax.YLim= ax.YLim - interval ;
    case 'W'
        ax.YLim= R*(ax.YLim -mean(ax.YLim))+ mean(ax.YLim)   ;
    case 'S'
        ax.YLim= (ax.YLim -mean(ax.YLim))/R + mean(ax.YLim)   ;
        %------------
    case 'e'
        if isfield(ax.UserData, 'NoZ');return ;end
        ax.ZLim= ax.ZLim + interval ;
    case 'd'
        if isfield(ax.UserData, 'NoZ');return ;end
        ax.ZLim= ax.ZLim - interval ;
    case 'E'
        if isfield(ax.UserData, 'NoZ');return ;end
        ax.ZLim= R*(ax.ZLim -mean(ax.ZLim))+ mean(ax.ZLim)   ;
    case 'D'
        if isfield(ax.UserData, 'NoZ');return ;end
        ax.ZLim= (ax.ZLim -mean(ax.ZLim))/R + mean(ax.ZLim)   ;
        
        %---------------
    case 'x'
        axis equal;
        fprintf('set axis equal \n')
    case 'X'
        axis auto
        fprintf('set axis auto \n')
    case 'P'
        fprintf('Capturing..... \n') ;
        print('-r100',gcf,'SnapShotFromEditBundle_r100','-dpng')  ;
        fprintf('Captured current figure. \n') ;
    case 'p'
        fprintf('Capturing..... \n') ;
        
        print('-r100',gcf,'SnapShotFromEditBundle_r100','-dpng')  ;
        fprintf('Captured current figure. \n') ;
    case 'n'
        axis normal
        fprintf('set axis normal \n')
    case 'N'
        axis normal
        fprintf('set axis normal \n')
        
        %--------------------
        
end




end