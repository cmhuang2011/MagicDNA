function readPoints_2(src,evn,allGraphics,ss_STEP,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nargin

% clc ;
%---------------
if nargin==4
opts.Interpreter = 'tex';
opts.Default = 'Cartesian grids';
quest = '\fontsize{16} Call the Spline sktech GUI (freeform) or use Cartesian/cylindrical grids';
answer = questdlg(quest,'Importing points',...
    'Free-form','Cartesian','Cylindrical',opts) ;
else
answer='ComingFromfreeformcurve';
    
end

% answer = questdlg(quest,'Importing points',...
%     'XYZ txt','Cartesian','Cylindrical','Free-form',opts) ;

%----------
if strcmp(answer,'XYZ txt')
    [FileName,PathName] = uigetfile('*.txt','Select the file with XYZ points');
    fffid=fopen(strcat(PathName,FileName));
    %         Oneline = fgetl(fffid);
    count=0;
    %     str2num(Oneline)
    
    while 1
        Oneline = fgetl(fffid) ;
        count=count+1 ;
        if isempty(Oneline)  % ||strcmp(Oneline, '-1')  ||count>20
            break
        end
        if isa(Oneline,'double')
            if Oneline== -1
                break
            end
        end
        
    end
    fclose(fffid);
    %-------------------------
    XYZ = zeros( count-1,3) ;
    fffid2=fopen(strcat(PathName,FileName)); count=0 ;
    while 1
        Oneline = fgetl(fffid2) ;
        count=count+1 ;
        if ~isa(Oneline,'double')
            XYZ(count ,: ) = str2num( Oneline) ;
        end
        if isempty(Oneline)  % ||strcmp(Oneline, '-1')  ||count>20
            break
        end
        if isa(Oneline,'double')
            if Oneline== -1
                break
            end
        end
    end
    fclose(fffid2);
    %------------------------------------------------
elseif strcmp(answer,'Cartesian')
    prompt = {'Spacing length in X','# of Instances in X','Spacing length in Y','# of Instances in Y','Spacing length in Z','# of Instances in Z'};
    dlgtitle = 'Input';
    dims = ones(1,6);
    definput = {'10','8','10','1','10','8'};opt.FontSize=10;
    answer = inputdlg(prompt,dlgtitle,dims,definput,opt) ;
    XYZSet =str2double(answer) ;
    dX= linspace(0, XYZSet(1)*(XYZSet(2)-1), XYZSet(2) ) ;
    dY= linspace(0, XYZSet(3)*(XYZSet(4)-1), XYZSet(4) ) ;
    dZ= linspace(0, XYZSet(5)*(XYZSet(6)-1), XYZSet(6) ) ;
    [X,Y,Z] = meshgrid(dX,dY,dZ) ;
    XYZ=[reshape(X,numel(X),1),reshape(Y,numel(Y),1),reshape(Z,numel(Z),1)] ;
elseif strcmp(answer,'Cylindrical')
    prompt = {'R spacing','R N Instance','Theta N Instance','Z spacing','Z N Instance'};
    dlgtitle = 'Input';
    dims = ones(1,5);
    definput = {'10','2','12','10','2'};opt.FontSize=10;
    answer = inputdlg(prompt,dlgtitle,dims,definput,opt) ;
    Cyl_paras =str2double(answer) ;
    dR= linspace(Cyl_paras(1),Cyl_paras(1)*(Cyl_paras(2)), Cyl_paras(2) ) ;
    dTheta= linspace(0,360-360/Cyl_paras(3),Cyl_paras(3))   ;
    dZ= linspace(0, Cyl_paras(4)*(Cyl_paras(5)-1), Cyl_paras(5) ) ;
    
    dX=dR'*cosd(dTheta) ; dX=dX(:) ;
    dY=dR'*sind(dTheta) ;  dY=dY(:) ;
    
    [X,Y,Z] = meshgrid(dX,dY,dZ) ;
    XYZ= [repelem([dX,dY],length(dZ),1) ,   repmat(dZ',size(dX,1),1) ] ;
    
    XYZ=[[repelem([0,0],length(dZ),1),dZ'] ; XYZ];
    %      XYZ=[reshape(X,numel(X),1),reshape(Y,numel(Y),1),reshape(Z,numel(Z),1)] ;
    
    % sdfsdf=3
elseif strcmp(answer,'Free-form')
    if ~isfield(ss_STEP.UserData ,'UseFFCbefore')
           
    xyz = [ 0 0 0  ; 1 5 10 ] ;
    FFC  = freeformcurve(1, xyz,ss_STEP)   ;
    FFC.VisualizeFFC ;
    ss_STEP.UserData.UseFFCbefore=1 ;
    else
     FFC=   ss_STEP.UserData.FFC ;
     FFC.ss_STEP =ss_STEP ;
     FFC.VisualizeFFC ;   
    end
    
    
    uiwait(FFC.fH);
    
    return
    
end


% XYZ  ;



% MainScreeSize = get(0,'screensize');% MainScreeSize(3)=round(0.5*MainScreeSize(3)) ;
% fAskLIne= figure('name','Assign lines','Position',MainScreeSize,'numbertitle','on'); clf;
fAskLIne= figure('name','Assign lines','numbertitle','on'); clf;
drawnow;warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFig = get(handle(fAskLIne), 'JavaFrame'); jFig.setMaximized(true); drawnow;
% fAskLIne_keypress
fAskLIne.KeyPressFcn=@(src,evn)fAskLIne_keypress(src,evn);
% fAskLIne=figure(3235) ; clf ;




sH = scatter3(    XYZ(:,1), XYZ(:,2), XYZ(:,3),126 ,'o', 'filled' ) ; hold on ;
sH.CData = repmat([0,0,1 ], length(sH.XData) ,1) ;
sH.ButtonDownFcn=@(src,evn)HitScatter2(src,evn) ;
sH.UserData.XYZ =XYZ ;
text( XYZ(:,1), XYZ(:,2), XYZ(:,3) ,strcat( '\leftarrow', num2str( (1:size(XYZ ,1))' )) ,'HitTest','off' ,'Clipping','on' ,'Visible','off' ) ;
ax=gca; ax.Position(3) =0.6 ;
axis auto ;
drawnow;
axis equal ;

mXYZ=mean(XYZ) ;
pH = plot3(mXYZ(1),mXYZ(2),mXYZ(3),'LineWidth',2,'HitTest','off' )  ; pH.UserData.Initial =1 ;

title('\color[rgb]{1 0 0}LeftClick  \color[rgb]{0 1 0}RightClick  ');

btn_AddLine = uicontrol(fAskLIne,'Style', 'pushbutton', 'String', 'Add Line','Unit','normalized', 'Position', [0.85 0.45 0.1 0.1] );
btn_AddLine.Callback=@(src,evn)Add_line(src,evn,sH,pH) ;

btn_FinishLine = uicontrol(fAskLIne,'Style', 'pushbutton', 'String', 'Send to Sketch','Unit','normalized', 'Position', [0.85 0.3 0.1 0.1] );
btn_FinishLine.Callback=@(src,evn)Export(src,evn,sH,pH,ss_STEP ,fAskLIne) ;
axLim= [ax.XLim , ax.YLim ,ax.ZLim];
axLim=axLim*1.1 ;
ax.XLim=axLim(1:2) ; ax.YLim=axLim(3:4) ; ax.ZLim=axLim(5:6) ;
uistack(sH,'top');

%---------add legend as instructions, single row
ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
hLg= legend(ForLegend,'Click me for instructions','Location','northwest' ) ;

hLg.String={'\bf{Click me for instructions}'};
hLg.Interpreter='tex';  %latex
hLg.Orientation='horizontal';
ForLegend.Marker='.' ; ForLegend.Marker='none';
hLg.ButtonDownFcn=@(src,evn)LegendBoxing_ReadPoints( src,evn,ax );

hLg.Units='normalized'; hLg.AutoUpdate ='off';
hLg.Position=[0.0063 0.9728 0.1569 0.025];
%------------------

uiwait(fAskLIne);

end

function fAskLIne_keypress(src,evn)
switch evn.Key
    case 'h'
        implay('PointToLines.mp4');
end
end

function Export(src,evn,sH,pH,ss_STEP ,fAskLIne)

ss_STEP.UserData.UsePoints.pXYZ=sH.UserData.XYZ ;
Inds =pH.UserData.Inds  ;
IIFlip=Inds(:,2)>Inds(:,1) ;
Inds(IIFlip,:) =flip(Inds(IIFlip,:) ,2) ;

ss_STEP.UserData.UsePoints.Edges= Inds ;

ss_STEP.UserData.UsePoints.Edges ;
close(fAskLIne) ;
end


function Add_line(src,evn,sH, pH )
OriRed =find ( and(and(sH.CData(:,1)==1 ,  sH.CData(:,2)==0) , sH.CData(:,3)==0)) ;
OriGreen =find ( and(and(sH.CData(:,1)==0 ,  sH.CData(:,2)==1) , sH.CData(:,3)==0)) ;
XYZ = sH.UserData.XYZ ;
if ~isempty(OriRed)  &&   ~isempty(OriGreen)
    if  pH.UserData.Initial ==1
        QQ= XYZ([OriRed ,OriGreen], :) ;  mQQ =mean(QQ) ;
        pH.XData=QQ(:,1);  pH.YData=QQ(:,2) ;pH.ZData=QQ(:,3)  ;
        pH.UserData.Inds = [OriRed , OriGreen] ;
        pH.UserData.Initial =0 ;
        sH.CData = repmat([0,0,1 ], length(sH.XData) ,1) ;
        str=strcat('L_{',num2str(OriRed) ,num2str(OriGreen),'}') ;
        %         text(mQQ(1) ,mQQ(2),mQQ(3) , str ,'HitTest','off','Clipping','on' ) ;
        
    else
        if  ismember([OriGreen,OriRed] , pH.UserData.Inds,'rows') ||  ismember([OriRed,OriGreen] , pH.UserData.Inds,'rows')
            fprintf('Edge Repeat.  no action. \n' )
            return
        end
        QQ= XYZ([OriRed ,OriGreen], :) ;  mQQ =mean(QQ) ;
        pH.XData(end+1:end+3)=[nan QQ(:,1)'];  pH.YData(end+1:end+3)=[nan QQ(:,2)'] ;pH.ZData(end+1:end+3)=[nan QQ(:,3)']  ;
        pH.UserData.Inds(end+1,:)  = [OriRed , OriGreen] ;
        sH.CData = repmat([0,0,1 ], length(sH.XData) ,1) ;
        str=strcat('L_{',num2str(OriRed) ,num2str(OriGreen),'}') ;
        %         text(mQQ(1) ,mQQ(2),mQQ(3) , str ,'HitTest','off','Clipping','on' ) ;
    end
else
    fprintf('Select Red and Green Nodes with mouse.\n')
    
end
end


function HitScatter2(src,evn )
d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2  + (src.ZData-evn.IntersectionPoint(3)).^2 ;
Ind = find(d==min(d))  ;
if evn.Button==1
    OriRed =find ( and(and(src.CData(:,1)==1 ,  src.CData(:,2)==0) , src.CData(:,3)==0)) ;
    if ~isempty(OriRed)
        src.CData(OriRed,:) =[0,0,1 ] ;
    end
    src.CData(Ind,:) =[1,0,0] ;
elseif evn.Button==3
    OriGreen =find ( and(and(src.CData(:,1)==0 ,  src.CData(:,2)==1) , src.CData(:,3)==0)) ;
    if ~isempty(OriGreen)
        src.CData(OriGreen,:) =[0,0,1 ] ;
    end
    src.CData(Ind,:) =[0,1,0] ;
end
end