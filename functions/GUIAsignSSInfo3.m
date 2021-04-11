function ssOption = GUIAsignSSInfo3( HyperB,fH,Prev_ssDNA,assemblyAxs_View )
% This is the GUI for assign details about ss connetion and tuning the
% orientation
%-
%   Detailed explanation goes here

% fH= figure('name','single-strand options','Position',[300 300 800 600],'numbertitle','off');

%     d = dialog('Position',[300 300 900 250],'Name','sing
% le-strand options');

% sdf=3
% Prev_ssDNA.BundleShiftTwoSide
tic

ssOption=[];
if ~isempty(Prev_ssDNA) ; ssOption= Prev_ssDNA ;   end

d = figure('Position',[300 300 900 250],'Name','single-strand options');
d.Units='normalized';d.OuterPosition=[0 0 1 1];  clf;
aH=gca;
aH.Position=[0.1,0.1,0.75,0.85]; axis equal;


drawnow;
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFig = get(handle(d), 'JavaFrame');
jFig.setMaximized(true);
drawnow;

str1='Length (unit: nt) of extended scaffold ends to form loops for preventing base stacking ';
txt1 = uicontrol('Parent',d,'Style','text','Units','normalized',...
    'Position',[0.85 0.9 0.1 0.05], 'String',str1 ,'Enable','Inactive','HorizontalAlignment','left');

str = cellfun(@num2str,num2cell(3:10),'uniformoutput',0);  % single strand options
popup1 = uicontrol('Parent',d , 'Style','popup','Units','normalized',...
    'Position',[0.85 0.85 0.1 0.05], 'String',str);


if ~isempty(Prev_ssDNA) ;  popup1.Value= Prev_ssDNA.OverRoute-2 ;   end

FBtwosideAsEnd=zeros(size(HyperB.ForcedConnectList,1),2);
for col=1:2
    for kc1=1:size(FBtwosideAsEnd,1)
        BS=HyperB.ForcedConnectList(kc1,3*col-2:3*col);
        BS(2)=-BS(2);
%         ComList=HyperB.containBundle{BS(1)}.ExternalXoverAsFB(:,1:2);
        ComList=[HyperB.containBundle{BS(1)}.ExternalXoverAsFB(:,1:2);HyperB.containBundle{BS(1)}.ExternalXoverAsFB(:,3:4) ];
        if ismember(BS(2:3),ComList,'rows')
            FBtwosideAsEnd(kc1,col)=1;
        end
    end
end

SortFB=zeros(sum(sum(FBtwosideAsEnd)),4);
[indu,indv]=find(FBtwosideAsEnd);

for endFBi=1:sum(sum(FBtwosideAsEnd))
    SortFB(endFBi,1:3)=HyperB.ForcedConnectList(indu(endFBi),3*indv(endFBi)-2:3*indv(endFBi));
end

SortFB=sortrows(SortFB);
Datat=mat2cell(SortFB, ones(1,size(SortFB,1) ) , [3,1]  );

%         Datat
for rev=1:size(Datat,1)
    Datat{rev,1}= num2str(Datat{rev,1});
end
str = cellfun(@num2str,num2cell(0:10),'uniformoutput',0);
t = uitable(d,'Units','normalized','Position',[0.8 0.2 0.2 0.5],'ColumnFormat',({[]  str}),...
    'ColumnEditable', true,'Data',Datat,'ColumnName',{'Bundle-Cyl-Base'; 'Clearance'},'Visible','off');

oldData=1 ;  % prevent data structure is old and not compatible.
if ~isempty(Prev_ssDNA)
    if ~isempty(Prev_ssDNA.individual)
        if length( t.Data(:,2))==length(Prev_ssDNA.individual(:,4))   % 02132019
            t.UserData.EndsPM=  Prev_ssDNA.individual(:,1:3) ;
            t.Data(:,2)= num2cell(Prev_ssDNA.individual(:,4));
            oldData= 0 ;
        end
    end
end

btn1 = uicontrol('Parent',d,'Units','normalized', 'Position',[0.9 0.05 0.08 0.1],   'String','OK');

str22 = cellfun(@num2str,num2cell([0:10,12,14,16]),'uniformoutput',0);  % single strand options
popup3 = uicontrol('Parent',d , 'Style','popup','Units','normalized','Position',[0.8 0.15 0.1 0.05], 'String',str22);



btn2 = uicontrol('Parent',d,'Units','normalized','Position',[0.8 0.05 0.08 0.1],  'String','Change selected ss lengths');

str2= cellfun(@num2str,num2cell(0:5),'uniformoutput',0);
popup2 = uicontrol('Parent',d ,  'Style','popup','Units','normalized','Position',[0.85 0.65 0.1 0.1],    'String',str2);
popup2.Callback=@(src,evn)SetAll(src,evn,t) ;




AllPatchH=findobj(fH,'type','patch');
new_PatchH = copyobj(AllPatchH,aH);
%     new_PatchH(1).HitTest
for patchi=1:length(new_PatchH)
    new_PatchH(patchi).FaceAlpha=0.2;
    new_PatchH(patchi).HitTest='off';
     new_PatchH(patchi).PickableParts='none' ;
end
AlllineH=findobj(fH,'type','Line');
%     TwoAxiesinFH=findobj(fH.Children,'type','Axes') ;

checkH_View = uicontrol('Style', 'checkbox','String', 'Axis auto/equal','Unit','normalized','Position', [0.85 0.75 0.1 0.1]);
checkH_View.Callback=@(src,evn)checkFcn2(src,evn,aH)  ;
new_handle=[];
% new_handle = copyobj(AlllineH,aH);
% for ks1=1:length(AlllineH)
%     if strcmp( AlllineH(ks1).Parent.UserData,'ax2')
%         new_handle(ks1).Visible='off';
%     elseif  sum( AlllineH(ks1).Color==[0.2 0.4 0.2])==3   %cylinder lines
%     else
%         new_handle(ks1).Color=[0,0,0];
%         new_handle(ks1).Visible='off';
%     end
%     new_handle(ks1).HitTest='off';
% end
checkH_Length = uicontrol('Style', 'checkbox','String', 'Show ssDNA lengths','Unit','normalized','Position', [0.85 0.75 0.1 0.1]);

%            [0.2 0.4 0.2]
hold on;
SaveGHelix=cell(1,length(HyperB.containBundle));
for Bundlei=1:length(HyperB.containBundle)
    QQWithPM10Bases =HyperB.containBundle{Bundlei}.HelixXYZG;
    SaveGHelix{Bundlei}= QQWithPM10Bases;
end
XXX=t.Data(:,1);
YYY=zeros(size(XXX,1),3);
ZZZ=zeros(size(XXX,1),3);
for kk=1:size(YYY,1)
    QQ= str2num( XXX{kk});  Q0=QQ;
    if HyperB.containBundle{QQ(1)}.Zbase1(QQ(2))==QQ(3)+2; QQ(3)=QQ(3)+2; end
    if HyperB.containBundle{QQ(1)}.Zbase2(QQ(2))==QQ(3)-2; QQ(3)=QQ(3)-2; end
    YYY(kk,:)= QQ;
    ZZZ(kk,:)= Q0;
    t.Data{kk,1}=num2str(QQ);
end
t.UserData.Endnodes=YYY; t.UserData.EndsPM2=ZZZ;

title('\color[rgb]{1 0 0}Selected(Left)  \color[rgb]{0 0 1}Non-Selected(Right)  \color[rgb]{0.2 0.2 0.2}non-changeable');

ppH=cell(1,size(HyperB.ForcedConnectList,1));    % connection plot handle
for fconni=1:size(HyperB.ForcedConnectList,1)
    OneConnec=HyperB.ForcedConnectList(fconni,:); CC=OneConnec ;
    if HyperB.containBundle{CC(1)}.Zbase1(CC(2))==CC(3)+2; OneConnec(3)=CC(3)+2; end
    if HyperB.containBundle{CC(1)}.Zbase2(CC(2))==CC(3)-2; OneConnec(3)=CC(3)-2; end
    if HyperB.containBundle{CC(4)}.Zbase1(CC(5))==CC(6)+2; OneConnec(6)=CC(6)+2; end
    if HyperB.containBundle{CC(4)}.Zbase2(CC(5))==CC(6)-2; OneConnec(6)=CC(6)-2; end
    % revise +- 2bp for variable OneConnec
    PACell=HyperB.containBundle{OneConnec(1)}.findHelixQ(OneConnec(2),OneConnec(3));
    PA=PACell{1};
    PBCell=HyperB.containBundle{OneConnec(4)}.findHelixQ(OneConnec(5),OneConnec(6));
    PB=PBCell{1};
    
    ppH{fconni}=plot3([PA(1);PB(1)], [PA(2);PB(2)], [PA(3);PB(3)] ,'b','LineWidth',1) ;
    ppH{fconni}.UserData.BCB0=OneConnec;
    ppH{fconni}.UserData.BCBCur=OneConnec;
    if ~ismember(OneConnec(1:3),YYY,'rows') && ~ismember(OneConnec(4:6),YYY,'rows')  %Both End nodes not selectable
        ppH{fconni}.HitTest= 'off'  ;
        ppH{fconni}.Color=[0.2,0.2,0.2];
    end
    
    ppH{fconni}.UserData.labelL=text(mean( ppH{fconni}.XData ),mean( ppH{fconni}.YData ),mean( ppH{fconni}.ZData ),'\leftarrow 0', 'clipping', 'on','FontSize',16 ) ;
    ppH{fconni}.UserData.labelL.HitTest= 'off'  ;
end
for kss = 1:length(ppH)
    ppH{kss}.ButtonDownFcn=@(src,evn)selectConnection(src,evn,ppH);
end
    
% plotHelixH=cell( length(HyperB.containBundle) , 20);
SaveXYZforPlot3= zeros(10000,3) ; nS =1 ;
for iB=1:length(HyperB.containBundle)
    for jCyl=1:length(SaveGHelix{iB})
        Qs= SaveGHelix{iB}{jCyl}; QQ=[Qs(11:end-10,1:3);nan,nan,nan] ;
%         plotHelixH{iB,jCyl}=  plot3(Qs(11:end-10,1) ,Qs(11:end-10,2), Qs(11:end-10,3),'.-k' ,'HitTest','off','PickableParts','none'  );
%         plotHelixH{iB,jCyl}.UserData.PXYZ=[Qs(11:end-10,1),Qs(11:end-10,2),Qs(11:end-10,3)];
%         plotHelixH{iB,jCyl}.UserData.TwoEndBP=[HyperB.containBundle{iB}.Zbase1(jCyl),HyperB.containBundle{iB}.Zbase2(jCyl)];
        SaveXYZforPlot3(nS:nS+size(QQ,1)-1 , : ) = QQ; nS=nS+size(QQ,1) ;
    end
end
SaveXYZforPlot3=SaveXYZforPlot3(1:nS-1,:) ;
plotHelixH=  plot3( SaveXYZforPlot3(:,1) ,SaveXYZforPlot3(:,2), SaveXYZforPlot3(:,3),'.-k' ,'HitTest','off','PickableParts','none'  ) ;
 
btn2.Callback=@(src,evn)UpdateSSConnection(src,evn,t,popup3,ppH);

Str={'--'}; for k2=1:length(HyperB.containBundle); Str{k2}=strcat('Bundle  ',num2str(k2));     end

popup_bundle = uicontrol('Parent',d , 'Style','popup','Units','normalized',...
    'Position',[0.8 0.15 0.1 0.05], 'String',Str);
sld1 = uicontrol('Style', 'slider','Parent',d,'Units','normalized',....
    'Min',0,'Max',20,'Value',0,'Position', [0.7 0.65 0.25 0.05],'Visible','off' );    %  blow
sld2 = uicontrol('Style', 'slider','Parent',d,'Units','normalized',....
    'Min',0,'Max',20,'Value',1,'Position', [0.7 0.65 0.25 0.05] ,'Visible','off');    %  blow

%-------July 17, add descriptions
str5='Extend single stranded scaffold at all ends. ';
txt5 = uicontrol('Parent',d,'Style','text','Units','normalized',...
    'Position',[0.88 0.46 0.1 0.05], 'String',str5,'Enable','Inactive','HorizontalAlignment','left');

str6='Specify single stranded scaffold lengths on selected connections';
txt6 = uicontrol('Parent',d,'Style','text','Units','normalized',...
    'Position',[0.88 0.34 0.1 0.05], 'String',str6,'Enable','Inactive','HorizontalAlignment','left');

%---



btn1.Position=[0.88 0.05 0.09 0.1] ;  % ok 
btn2.Position=[0.88 0.17 0.09 0.1] ;   % change 
t.Position=[   0.88 0.2 0.18 0.3] ;   % hide 

popup1.Position=[ 0.88 0.85 0.09 0.05] ;   % scaffold loop
txt1.Position=[0.88 0.9 0.09 0.05];      % description loop
popup2.Position=[ 0.88 0.4 0.09 0.05] ;    % apply all
popup3.Position=[ 0.88 0.28 0.09 0.05] ;  % apply selected connections

checkH_Length.Position=[ 0.88 0.6 0.09 0.05] ;checkH_Length.Value = 1 ;
checkH_View.Position= [0.88 0.7 0.09 0.05] ;

popup_bundle.Position=[0.88 0.72 0.09 0.05];  % show bundle 
sld1.Position=[0.8 0.7 0.1 0.03] ;  % hide
sld2.Position=[0.8 0.6 0.1 0.03] ;  % hide
popup_bundle.Callback=@(src,evn)SelectBundlePop(src,evn,HyperB,new_PatchH,sld1,sld2)  ;

% sld1.Callback= @(src,evn)slds1Fcn(src,evn,HyperB,ppH,plotHelixH,new_handle,popup_bundle,sld2,t) ;
% sld2.Callback= @(src,evn)slds2Fcn(src,evn,HyperB,ppH,plotHelixH,new_handle,popup_bundle,sld1,t) ;
% slider function has been deactive. Plot the helixes with single plot to
% improve graphical performance

txt1.ButtonDownFcn=@(src,evn)uisetfont(txt1);
txt5.ButtonDownFcn=@(src,evn)uisetfont(txt5);
txt6.ButtonDownFcn=@(src,evn)uisetfont(txt6);



d.UserData=zeros(length(HyperB.containBundle),2);
if ~isempty(Prev_ssDNA)
    d.UserData=  Prev_ssDNA.BundleShiftTwoSide;
end


btn1.Callback=@(src,evn)Exportresult(src,evn,t,popup1,ppH) ;

checkH_Length.Callback=@(src,evn)CheckFcnshowssDNAL(src,evn,t,ppH) ;
t.CellEditCallback=@(src,evn)tableChange(src,evn,ppH)  ;
% toc1 =toc

if ~isempty(Prev_ssDNA) && oldData==0
    t.UserData.EndsPM=  Prev_ssDNA.individual(:,1:3) ;
    t.Data(:,2)= num2cell(Prev_ssDNA.individual(:,4));
end


for seleB= 1: length(HyperB.containBundle)
    popup_bundle.Value=seleB ;
    SelectBundlePop(popup_bundle,[],HyperB,new_PatchH,sld1,sld2)  ;
%     slds1Fcn(sld1,[],HyperB,ppH,plotHelixH,new_handle,popup_bundle,sld2,t) ;
%     slds2Fcn(sld2,[],HyperB,ppH,plotHelixH,new_handle,popup_bundle,sld1,t)
end

% for kk=1:length(ppH)
%  uistack(ppH{kk},'top');
% end
% toc2 =toc
% 
% Cstr    = cell(1, length(ppH));
% Cstr(:) = {'top'};
% 
% toc3 =toc
% 
% %  uistack(ppH{:},'top');
% cellfun(@uistack,ppH,Cstr) ;
% toc4 =toc

tableChange(t,[],ppH) ;
d.KeyPressFcn=@(src,evn)ChangeAxesLimit(src,evn) ;
xlabel('X') ; ylabel('Y') ; zlabel('Z') ; %box on
set(gca,'View',assemblyAxs_View);

AssignIcon( btn2,'ChangeSelectedSsDNA.jpg' );
AssignIcon( btn1,'OKssDNA.jpg' );
popup3.TooltipString = 'Assign single-stranded lengths to the selected connections';
popup2.TooltipString = 'Apply to all ends';
popup1.TooltipString = 'Extension as Scaffold loop ';
popup_bundle.TooltipString= 'Use to indicate bundles ' ;
btn2.TooltipString='Update the selected connections(red) into designated lengths of ssDNA. ' ;
btn1.TooltipString='If all connections are good, complete this step of ssDNA and ready to do scaffold routing.' ;

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
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_ssScaf( src,evn,aH );
            %             fHH.UserData.LgIndication.Title.Interpreter='latex';
            
            %             drawnow;
            hLg.Units='normalized'; hLg.AutoUpdate ='off';
            hLg.Title.String='Click me for instructions' ;
            hLg.Position=[0.0063 0.9528 0.1569 0.0387];
%---------------

uiwait(d);

% sdf=3
%-------------------------------------------------------End main

    function ChangeAxesLimit(src,evn)
        ax=gca;
        interval =5; R=1.1;
        switch  evn.Character
            case 'h'
            implay('ssScaf.mp4');
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
                ax.ZLim= ax.ZLim + interval ;
            case 'd'
                ax.ZLim= ax.ZLim - interval ;
            case 'E'
                ax.ZLim= R*(ax.ZLim -mean(ax.ZLim))+ mean(ax.ZLim)   ;
            case 'D'
                ax.ZLim= (ax.ZLim -mean(ax.ZLim))/R + mean(ax.ZLim)   ;
            case 'P'
                fprintf('Capturing current figure \n') ;
                print('-r100',gcf,'SnapShotFromSSDNA_r100','-djpeg')  ;
            case 'p'
                fprintf('Capturing current figure \n') ;
                print('-r100',gcf,'SnapShotFromSSDNA_r100','-djpeg')  ;
            case 'x'
                axis equal;
                fprintf('set axis equal \n')
            case 'X'
                axis auto
                fprintf('set axis auto \n')
            case 'n'
                axis normal
                fprintf('set axis normal \n')
            case 'N'
                axis normal
                fprintf('set axis normal \n')
                
        end
        
    end


    function tableChange(src,evn,ppH)
        for k_ss=1: length(ppH)
            BCB0 =ppH{k_ss}.UserData.BCB0  ;
            if ismember(BCB0(1:3),t.UserData.Endnodes, 'rows')
                [~,indA]=ismember(BCB0(1:3) ,t.UserData.Endnodes, 'rows');
                L_A= t.Data{indA,2} ;
            else
                L_A=0;
            end
            if ismember(BCB0(4:6),t.UserData.Endnodes, 'rows')
                [~,indB]=ismember(BCB0(4:6) ,t.UserData.Endnodes, 'rows');
                L_B= t.Data{indB,2} ;
            else
                L_B=0;
                
            end
            %              if isvarname(L_A) &&   isvarname(L_B)
            ppH{k_ss}.UserData.labelL.String=strcat('\leftarrow ',num2str( L_A+ L_B) ) ;
            %              end
        end
    end  % end fcn tableChange

    function CheckFcnshowssDNAL(src,evn,t,ppH)
        %         L_ssConn= zeros(length(ppH) ,1) ;
        %         for k_ss=1: length(L_ssConn)
        %              BCB0 =ppH{k_ss}.UserData.BCB0  ;
        %              if ismember(BCB0(1:3),t.UserData.Endnodes, 'rows')
        %                 [~,indA]=ismember(BCB0(1:3) ,t.UserData.Endnodes, 'rows');
        %                 L_A= t.Data{indA,2} ;
        %              end
        %              if ismember(BCB0(4:6),t.UserData.Endnodes, 'rows')
        %                 [~,indB]=ismember(BCB0(4:6) ,t.UserData.Endnodes, 'rows');
        %                 L_B= t.Data{indB,2} ;
        %              end
        %              L_ssConn(k_ss)= L_A+ L_B;
        %             ppH{k_ss}.UserData.labelL.String=strcat('\leftarrow ',num2str( L_A+ L_B) ) ;
        %         end
        if src.Value ==1
            for k_ss=1:length(ppH)
                ppH{k_ss}.UserData.labelL.Visible='on';
            end
        else
            for k_ss=1:length(ppH)
                ppH{k_ss}.UserData.labelL.Visible='off';
            end
        end
        
    end  % end fcn CheckFcnshowssDNAL

    function slds1Fcn(src,evn,HyperB,ppH,plotHelixH,new_handle,popup_bundle,sld2,t)
        d=gcf;
        sBundle2=popup_bundle.Value;
        Part=HyperB.containBundle{sBundle2};
        medZbase= round(src.Value);
        medZbaseCC= round(sld2.Value);
        d.UserData(sBundle2,:)=[medZbase ,medZbaseCC];
        for Cylplot=1:length(new_handle)   %for cylinder plot
            if   strcmp(new_handle(Cylplot).Visible,'on')&&   new_handle(Cylplot).UserData(1)==sBundle2
                PlotXYZ2=Part.CylinderXYZGlobal;
                PlotXYZ=[PlotXYZ2(:,1),PlotXYZ2(:,4),PlotXYZ2(:,2),PlotXYZ2(:,5),PlotXYZ2(:,3),PlotXYZ2(:,6)];
                CylinBundle=new_handle(Cylplot).UserData(2);
                extrac=PlotXYZ(CylinBundle,:);
                OldInterval=[Part.Zbase1(CylinBundle) , Part.Zbase2(CylinBundle)];
                NewInt=OldInterval +[medZbase,0];
                v=[extrac(1:2:5);extrac(2:2:6)] ;
                vq1 = interp1(OldInterval,v,NewInt(1));
                new_handle(Cylplot).XData(1)=vq1(1);
                new_handle(Cylplot).YData(1)=vq1(2);
                new_handle(Cylplot).ZData(1)=vq1(3);
            end
        end
        for helx=1:length(Part.Zbase1)
            plotHelixH{sBundle2,helx}.XData =   plotHelixH{sBundle2,helx}.UserData.PXYZ(1+medZbase:end-medZbaseCC,1);
            plotHelixH{sBundle2,helx}.YData =   plotHelixH{sBundle2,helx}.UserData.PXYZ(1+medZbase:end-medZbaseCC,2);
            plotHelixH{sBundle2,helx}.ZData =   plotHelixH{sBundle2,helx}.UserData.PXYZ(1+medZbase:end-medZbaseCC,3);
        end
        
        for conni=1:length(ppH)
            BCB0= ppH{conni}.UserData.BCB0 ;
            LP= BCB0(1:3)  ;   RP= BCB0(4:6) ;
            if LP(1)==sBundle2 &&  abs(Part.Zbase1(LP(2))-LP(3))<20
                BCBCur=  ppH{conni}.UserData.BCBCur ;
                BCBCur(3)= BCB0(3)+medZbase;
                PACell=HyperB.containBundle{BCBCur(1)}.findHelixQ(BCBCur(2),BCBCur(3) ); PA=PACell{1};
                PBCell=HyperB.containBundle{BCBCur(4)}.findHelixQ(BCBCur(5),BCBCur(6));PB=PBCell{1};
                ppH{conni}.XData=[   PA(1),PB(1)];    ppH{conni}.UserData.labelL.Position(1)= mean([   PA(1),PB(1)]) ;
                ppH{conni}.YData=[   PA(2),PB(2)];   ppH{conni}.UserData.labelL.Position(2)= mean([   PA(2),PB(2)]) ;
                ppH{conni}.ZData=[   PA(3),PB(3)];     ppH{conni}.UserData.labelL.Position(3)= mean([   PA(3),PB(3)]) ;
                ppH{conni}.UserData.BCBCur =  BCBCur  ;
            end
            if RP(1)==sBundle2 && abs( Part.Zbase1(RP(2))-RP(3) )<20
                BCBCur=  ppH{conni}.UserData.BCBCur ;
                BCBCur(6)= BCB0(6)+medZbase;
                
                PACell=HyperB.containBundle{BCBCur(1)}.findHelixQ(BCBCur(2),BCBCur(3)); PA=PACell{1};
                PBCell=HyperB.containBundle{BCBCur(4)}.findHelixQ(BCBCur(5),BCBCur(6)); PB=PBCell{1};
                ppH{conni}.XData=[   PA(1),PB(1)];
                ppH{conni}.YData=[   PA(2),PB(2)];
                ppH{conni}.ZData=[   PA(3),PB(3)];
                ppH{conni}.UserData.BCBCur =  BCBCur  ;
                ppH{conni}.UserData.labelL.Position(1)= mean([   PA(1),PB(1)]) ;
                ppH{conni}.UserData.labelL.Position(2)= mean([   PA(2),PB(2)]) ;
                ppH{conni}.UserData.labelL.Position(3)= mean([   PA(3),PB(3)]) ;
                
            end
        end
        updateT(ppH,t)
    end % end fcn slds1Fcn

    function slds2Fcn(src,evn,HyperB,ppH,plotHelixH,new_handle,popup_bundle,sld1,t)
        d=gcf;
        sBundle2=popup_bundle.Value;
        Part=HyperB.containBundle{sBundle2};
        medZbase= round(src.Value);
        medZbaseCC= round(sld1.Value);
        d.UserData(sBundle2,:)=[medZbaseCC ,medZbase];
        for Cylplot=1:length(new_handle)   %for cylinder plot
            if   strcmp(new_handle(Cylplot).Visible,'on')&&   new_handle(Cylplot).UserData(1)==sBundle2
                PlotXYZ2=Part.CylinderXYZGlobal;
                PlotXYZ=[PlotXYZ2(:,1),PlotXYZ2(:,4),PlotXYZ2(:,2),PlotXYZ2(:,5),PlotXYZ2(:,3),PlotXYZ2(:,6)];
                CylinBundle=new_handle(Cylplot).UserData(2);
                extrac=PlotXYZ(CylinBundle,:);
                OldInterval=[Part.Zbase1(CylinBundle) , Part.Zbase2(CylinBundle)];
                NewInt=OldInterval -[0,medZbase];
                v=[extrac(1:2:5);extrac(2:2:6)] ;
                vq1 = interp1(OldInterval,v,NewInt(2));
                new_handle(Cylplot).XData(2)=vq1(1);
                new_handle(Cylplot).YData(2)=vq1(2);
                new_handle(Cylplot).ZData(2)=vq1(3);
            end
        end
        for helx=1:length(Part.Zbase1)
            plotHelixH{sBundle2,helx}.XData =   plotHelixH{sBundle2,helx}.UserData.PXYZ(1+medZbaseCC:end-medZbase,1);
            plotHelixH{sBundle2,helx}.YData =   plotHelixH{sBundle2,helx}.UserData.PXYZ(1+medZbaseCC:end-medZbase,2);
            plotHelixH{sBundle2,helx}.ZData =   plotHelixH{sBundle2,helx}.UserData.PXYZ(1+medZbaseCC:end-medZbase,3);
        end
        
        for conni=1:length(ppH)
            BCB0= ppH{conni}.UserData.BCB0 ;
            LP= BCB0(1:3)  ;   RP= BCB0(4:6) ;
            if LP(1)==sBundle2 && abs(Part.Zbase2(LP(2))-LP(3))<20
                BCBCur=  ppH{conni}.UserData.BCBCur ;
                BCBCur(3)= BCB0(3)-medZbase;
                PACell=HyperB.containBundle{BCBCur(1)}.findHelixQ(BCBCur(2),BCBCur(3) ); PA=PACell{1};
                PBCell=HyperB.containBundle{BCBCur(4)}.findHelixQ(BCBCur(5),BCBCur(6));PB=PBCell{1};
                ppH{conni}.XData=[   PA(1),PB(1)];    ppH{conni}.UserData.labelL.Position(1)= mean([   PA(1),PB(1)]) ;
                ppH{conni}.YData=[   PA(2),PB(2)];   ppH{conni}.UserData.labelL.Position(2)= mean([   PA(2),PB(2)]) ;
                ppH{conni}.ZData=[   PA(3),PB(3)];     ppH{conni}.UserData.labelL.Position(3)= mean([   PA(3),PB(3)]) ;
                ppH{conni}.UserData.BCBCur =  BCBCur  ;
            end
            if RP(1)==sBundle2 && abs(Part.Zbase2(RP(2))-RP(3))<20
                BCBCur=  ppH{conni}.UserData.BCBCur ;
                BCBCur(6)= BCB0(6)-medZbase;
                PACell=HyperB.containBundle{BCBCur(1)}.findHelixQ(BCBCur(2),BCBCur(3)); PA=PACell{1};
                PBCell=HyperB.containBundle{BCBCur(4)}.findHelixQ(BCBCur(5),BCBCur(6)); PB=PBCell{1};
                ppH{conni}.XData=[   PA(1),PB(1)];    ppH{conni}.UserData.labelL.Position(1)= mean([   PA(1),PB(1)]) ;
                ppH{conni}.YData=[   PA(2),PB(2)];   ppH{conni}.UserData.labelL.Position(2)= mean([   PA(2),PB(2)]) ;
                ppH{conni}.ZData=[   PA(3),PB(3)];     ppH{conni}.UserData.labelL.Position(3)= mean([   PA(3),PB(3)]) ;
                ppH{conni}.UserData.BCBCur =  BCBCur  ;
            end
        end
        updateT(ppH,t)
    end % end fcn slds2Fcn

    function updateT(ppH,t)
        for ppj2=1:length(ppH)
            if ismember(ppH{ppj2}.UserData.BCB0(1:3),t.UserData.Endnodes, 'rows')
                [~,ind]=ismember(ppH{ppj2}.UserData.BCB0(1:3),t.UserData.Endnodes, 'rows');
                t.Data{ind,1}=num2str(ppH{ppj2}.UserData.BCBCur(1:3));
            end
            
            if ismember(ppH{ppj2}.UserData.BCB0(4:6),t.UserData.Endnodes, 'rows')
                [~,ind]=ismember(ppH{ppj2}.UserData.BCB0(4:6),t.UserData.Endnodes, 'rows');
                t.Data{ind,1}=num2str(ppH{ppj2}.UserData.BCBCur(4:6));
            end
        end
    end % end fcn updateT



    function SelectBundlePop(src,evn,HyperB,new_Patch,sld1,sld2)
        d=gcf; sBundle=src.Value ;
        for k=1:length(new_Patch)
            if  strcmp(  new_Patch(k).Tag,num2str(sBundle) )
                new_Patch(k).FaceColor=[1,0,0];
                new_Patch(k).FaceAlpha=0.4;
            else
                new_Patch(k).FaceColor= [0.5,0.5,0.5];
                new_Patch(k).FaceAlpha=0.2;
            end
        end
        sld1.Value=d.UserData(sBundle,1);   % load selected bundle shift into slider
        sld2.Value=d.UserData(sBundle,2);
    end  % end fcn SelectBundlePop


    function UpdateSSConnection(src,evn,t,popup3,ppH)
        RedConn=[] ;
        for k=1:length(ppH)
            if sum( ppH{k}.Color==[1,0,0])==3
                RedConn=union(RedConn,k);
            end
        end
        ssBP= str2double(popup3.String{popup3.Value});
        for redi=1:length(RedConn)
            Connection=ppH{RedConn(redi)}.UserData.BCB0 ;
            if ismember(  Connection(1:3),t.UserData.Endnodes,'rows') && ismember(  Connection(4:6),t.UserData.Endnodes,'rows')
                [~,indA]=ismember(  Connection(1:3),t.UserData.Endnodes,'rows') ;
                [~,indB]=ismember(  Connection(4:6),t.UserData.Endnodes,'rows') ;
                t.Data{indA,2}= round(ssBP/2);
                t.Data{indB,2}= ssBP-round(ssBP/2) ;
            elseif ~ismember(  Connection(1:3),t.UserData.Endnodes,'rows') && ismember(  Connection(4:6),t.UserData.Endnodes,'rows')
                [~,indB]=ismember(  Connection(4:6),t.UserData.Endnodes,'rows') ;
                t.Data{indB,2}= ssBP  ;
            elseif ismember(  Connection(1:3),t.UserData.Endnodes,'rows') && ~ismember(  Connection(4:6),t.UserData.Endnodes,'rows')
                [~,indA]=ismember(  Connection(1:3),t.UserData.Endnodes,'rows') ;
                t.Data{indA,2}= ssBP;
            end
            ppH{RedConn(redi)}.Color=[0,0,1];
            ppH{RedConn(redi)}.UserData.labelL.String=strcat('\leftarrow',num2str(ssBP) ) ;
        end
    end   % end fcn UpdateSSConnection


    function selectConnection(src,evn,ppH)
        if  evn.Button==1
            src.Color=[1,0,0];
        elseif evn.Button==3
            src.Color=[0,0,1];
        end
%         ppH
        nR =0 ; nB=0 ;
        for k = 1 :length(ppH)
            if ~isempty(ppH{k})
                
                if sum(ppH{k}.Color ==[1,0,0] ) ==3
                    nR =nR+1 ;
                elseif sum(ppH{k}.Color ==[0,0,1] ) ==3
                    nB=nB+1 ;
                end
            end
        end
%             [nR , nB]
     strTT=strcat('\color[rgb]{1 0 0}Selected(Left) ', num2str(nR),' \color[rgb]{0 0 1}Non-Selected(Right)', num2str(nB),'  \color[rgb]{0.2 0.2 0.2}non-changeable' ) ;
     
     title(strTT);
%       title('\color[rgb]{1 0 0}Selected(Left)  \color[rgb]{0 0 1}Non-Selected(Right)  \color[rgb]{0.2 0.2 0.2}non-changeable');
  
    end  % end fcn UpdateSSConnection


    function SetAll(src,~,t)
        answer = questdlg('Are you sure to change all? This will apply to all connections.', ...
	'Change all ends', ...
	'No','Yes','Yes');
        switch answer
            case 'No'
                return
        end
        
        for k22=1:size(t.Data,1)
            t.Data{k22,2}=str2double(src.String{src.Value});
        end
        tableChange(t,[],ppH) ;
    end  % end fcn SetAll

    function Exportresult(src,~,t,popup1,ppH)
        ssOption.OverRoute=str2double(popup1.String{popup1.Value});
        %       yy= cell2mat( t.Data(:,2));
        ssOption.individual=[t.UserData.EndsPM2, cell2mat( t.Data(:,2))];
        ssOption.BundleShiftTwoSide=d.UserData;
        FCBCB=cell(size(ppH));         %force connection Bundle Cylinder Base Rep, before and after tuning
        for kkkk=1:length(FCBCB)
            [EndorSideA,indA]=  ismember(ppH{kkkk}.UserData.BCB0(1:3),t.UserData.Endnodes,'rows');
            [EndorSideB,indB]=  ismember(ppH{kkkk}.UserData.BCB0(4:6),t.UserData.Endnodes,'rows') ;
            FCBCB{  kkkk}=ppH{kkkk}.UserData;
            BCBPM2a=zeros(1,6);
            if EndorSideA==1
                BCBPM2a( 1:3)=t.UserData.EndsPM2(indA,:) ;
            else
                BCBPM2a( 1:3)=ppH{kkkk}.UserData.BCB0(1:3) ;
            end
            if EndorSideB==1
                BCBPM2a( 4:6)=t.UserData.EndsPM2(indB,:) ;
            else
                BCBPM2a( 4:6)=ppH{kkkk}.UserData.BCB0(4:6) ;
            end
            FCBCB{  kkkk}.BCBPM2=BCBPM2a;
        end
        ssOption.ForceConnUpdate=FCBCB;
        delete(gcf)
    end % end fcn Exportresult






end  %end of all


