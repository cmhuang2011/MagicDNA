f = figure; 
% fg=uifigure;
clf;
% f.Pointer ='fleur';

f.Name='MagicDNA' ;
TStrO='<html><b><i><font size=5>Square</h1></html>' ;  % ALL tab use html to set font size, graphics are possible to use in future.
% replace(TStrO,"Square","Hello")

drawnow;warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFig = get(handle(f), 'JavaFrame'); jFig.setMaximized(true); drawnow; % this line makes figure(gcf) always maximize the window for MATLAB 2018

tgroup = uitabgroup('Parent', f);   drawnow;
tab_part = uitab('Parent', tgroup, 'Title', replace(TStrO,"Square","Part"));
tab_mechanism = uitab('Parent', tgroup, 'Title', replace(TStrO,"Square","Mechanism"),'Tag','s_Mech');
tgroup.SelectedTab =tab_mechanism ;

s = sprintf('Old way to create bundles. \n Recommend to use Custom cross-sections and "Edit Bundle" to create bundle library.');
tab_part.TooltipString = s;


Sub_group_part = uitabgroup('Parent', tab_part);
ss_HC = uitab('Parent', Sub_group_part, 'Title', replace(TStrO,"Square","HoneyComb"));
ss_SQ = uitab('Parent', Sub_group_part, 'Title', TStrO);

Sub_group_mech = uitabgroup('Parent', tab_mechanism);
ss_STEP = uitab('Parent', Sub_group_mech, 'Title', replace(TStrO,"Square","Sketch"),'Tag','ss_STEP');
ss_Assembly = uitab('Parent', Sub_group_mech, 'Title',replace(TStrO,"Square","Assembly") ,'Tag','ss_Assembly');
ss_Assembly.ButtonDownFcn=@(src,evn)RecoverKeyMove(src,evn,1) ; 


ss_STEP.ButtonDownFcn =@(src,evn)ssSTEP_Button(src,evn);

%--------merged as Routing tab
ss_routing=uitab('Parent', Sub_group_mech, 'Title',  replace(TStrO,"Square","Routing"));
RoutingTab = uitabgroup('Parent', ss_routing);
% ss_RoutingParameter = uitab('Parent', RoutingTab, 'Title', replace(TStrO,"Square","parameter"),'Tag','ss_RoutingParameter' );




ss_Scaffold = uitab('Parent', RoutingTab, 'Title', replace(TStrO,"Square","scaffold"),'Tag','ss_Scaffold' );
ss_Scaffold.ButtonDownFcn =@(src,evn)ssSTEP_Button(src,evn);

ss_Staple = uitab('Parent', RoutingTab, 'Title', replace(TStrO,"Square","staple"),'Tag' , 'ss_Staple' );
ss_Staple.ButtonDownFcn =@(src,evn)ssSTEP_Button(src,evn);

ss_OH = uitab('Parent', RoutingTab, 'Title', replace(TStrO,"Square","overhang"));
ss_OH.ButtonDownFcn =@(src,evn)ssOH_Button(src,evn);

ss_json = uitab('Parent', RoutingTab, 'Title', replace(TStrO,"Square","seq & diagram"),'Tag','ss_json');
% ss_json.ButtonDownFcn =@(src,evn)ssJSON_Button(src,evn);
ss_json.ButtonDownFcn=@(src,evn)ssJSON_Button(src,evn) ; 


%-------

% ss_oxDNA = uitab('Parent', Sub_group_mech, 'Title', 'oxDNA');
ss_oxDNA = uitab('Parent', Sub_group_mech, 'Title',  replace(TStrO,"Square","Simulation"));

SubSub_group_oxDNA = uitabgroup('Parent', ss_oxDNA);
sss_Confforce= uitab('Parent', SubSub_group_oxDNA, 'Title',  replace(TStrO,"Square","Topogoly and Conf"),'Tag' , 'sss_Confforce' );
% sss_pattern= uitab('Parent', SubSub_group_oxDNA, 'Title',  replace(TStrO,"Square","pattern"),'Tag','sss_pattern' );

sss_Confforce.ButtonDownFcn =@(src,evn)ssSTEP_Button(src,evn);


sss_Traj= uitab('Parent', SubSub_group_oxDNA, 'Title',  replace(TStrO,"Square","Trajectory") ,'Tag' , 'sss_Traj' );
sss_Traj.ButtonDownFcn =@(src,evn)ssSTEP_Button(src,evn);



% ss_oxDNA = uitab('Parent', Sub_group_mech, 'Title', '<html><font size=8><img src="file:///C:\Users\huang.2011\Dropbox\ShareFolderDNA\OOP_caddom\single.jpg">oxDNA</h1></html>');

% <img src="file:///C:\Users\huang.2011\Dropbox\ShareFolderDNA\OOP_caddom\single.jpg">
%  t.ColumnName{3}='<html><font size=8>Enable</h1></html>' ;

%-------------

% test=tabdlg 
allaxes.HC_Cross=axes(ss_HC,'Position',[0.05 0.45 0.5 0.5 ]);title('CrossSection Editor'  );
allaxes.HC_Extrude=axes(ss_HC,'Position',[0.05 0.05 0.3 0.3 ],'Color',[0.8,0.8,0.8]);title('Extrude CrossSection'  );
allaxes.HC_CylModel=axes(ss_HC,'Position',[0.6 0.55 0.3 0.4 ]);title('Cylinder Model'  );
allaxes.HC_Section=axes(ss_HC,'Position',[0.8 0.05 0.15 0.25 ]);title('Section Check'  );

dy=0.05;
allbtn.HC_Show3D= uicontrol('Style','pushbutton','Parent',ss_HC,...
    'Units','normalized','Position',[0.37 0.26+dy 0.1  0.05],'String','Show3D');
allbtn.HC_Adjust= uicontrol('Style','pushbutton','Parent',ss_HC,...
    'Units','normalized','Position',[0.37 0.16+dy 0.1  0.05],'String','Adjust');
allbtn.HC_SaveSection= uicontrol('Style','pushbutton','Parent',ss_HC,...
    'Units','normalized','Position',[0.37 0.06 0.1  0.05],'String','SaveSection');
allbtn.HC_Clear= uicontrol('Style','pushbutton','Parent',ss_HC,...
    'Units','normalized','Position',[0.6 0.4 0.1  0.05],'String','Clear');
allbtn.HC_MakeColseLoop= uicontrol('Style','pushbutton','Parent',ss_HC,...
    'Units','normalized','Position',[0.7 0.4 0.1  0.05],'String','MakeColseLoop');
allbtn.HC_SavePart= uicontrol('Style','pushbutton','Parent',ss_HC,...
    'Units','normalized','Position',[0.83 0.36 0.1  0.05],'String','SavePart');
allbtn.HC_LoadPart= uicontrol('Style','pushbutton','Parent',ss_HC,...
    'Units','normalized','Position',[0.83 0.41 0.1  0.05],'String','LoadPart');
allbtn.HC_ClearPart= uicontrol('Style','pushbutton','Parent',ss_HC,...
    'Units','normalized','Position',[0.83 0.46 0.1  0.05],'String','ClearPart');

alledit.HC_show3Dlower= uicontrol('Style','edit','Parent',ss_HC,...
    'Units','normalized','Position',[0.37 0.21+dy 0.05 0.05 ],'UserData','show3Dlower','String',50);
alledit.HC_show3Dupper= uicontrol('Style','edit','Parent',ss_HC,...
    'Units','normalized','Position',[0.43 0.21+dy 0.05 0.05 ],'UserData','show3Dupper','String',100);
alledit.HC_Adjustlower= uicontrol('Style','edit','Parent',ss_HC,...
    'Units','normalized','Position',[0.37 0.11+dy 0.05 0.05],'UserData','Adjustlower','String',-5);
alledit.HC_Adjustupper= uicontrol('Style','edit','Parent',ss_HC,...
    'Units','normalized','Position',[0.43 0.11+dy 0.05 0.05 ],'UserData','Adjustupper','String',5);

allpop.HC_SaveSection = uicontrol('Style','popupmenu','Parent',ss_HC,...
    'Units','normalized','Position',[0.37 0 0.1  0.05],'Value',1,'String','New');
%--------
allbg.HC_Loop = uibuttongroup(ss_HC,'Position',[0.5 0.22 0.1  0.15],'Title','Loop') ;
allbg.HC_Unit = uibuttongroup(ss_HC,'Position',[0.5 0.05 0.1  0.15],'Title','Unit') ;
allbg.HC_Cylinder = uibuttongroup(ss_HC,'Position',[0.62 0.05 0.13  0.32],'Title','Cylinder') ;

allrb.HC_Loop_Outer = uicontrol('Style','radiobutton','Parent',allbg.HC_Loop,...
    'Units','normalized','Position',[0.05 0.3 0.8 0.2],'String','Outer Loop'  );
allrb.HC_Loop_Inner = uicontrol('Style','radiobutton','Parent',allbg.HC_Loop,...
    'Units','normalized','Position',[0.05 0.8 0.8 0.2],'String','Inner Loop'  );

allrb.HC_Unit_nm = uicontrol('Style','radiobutton','Parent',allbg.HC_Unit,...
    'Units','normalized','Position',[0.05 0.3 0.8 0.2],'String','nm'  );
allrb.HC_Unit_base = uicontrol('Style','radiobutton','Parent',allbg.HC_Unit,...
    'Units','normalized','Position',[0.05 0.8 0.8 0.2],'String','base'  );

allrb.HC_Cyl_Draw = uicontrol('Style','radiobutton','Parent',allbg.HC_Cylinder,...
    'Units','normalized','Position',[0.05 0.8 0.8 0.2],'String','Draw Outline'  );
allrb.HC_Cyl_Select = uicontrol('Style','radiobutton','Parent',allbg.HC_Cylinder,...
    'Units','normalized','Position',[0.05 0.5 0.8 0.2],'String','Select Cylinder'  );
allrb.HC_Cyl_AddOrDel = uicontrol('Style','radiobutton','Parent',allbg.HC_Cylinder,...
    'Units','normalized','Position',[0.05 0.2 0.8 0.2],'String','Add or Delete'  );


allaxes.SQ_Cross=axes(ss_SQ,'Position',[0.05 0.45 0.5 0.5 ]);title('CrossSection Editor'  );
allaxes.SQ_Extrude=axes(ss_SQ,'Position',[0.05 0.05 0.3 0.3 ],'Color',[0.8,0.8,0.8]);title('Extrude CrossSection'  ); axis equal; view(45,30) ;
allaxes.SQ_CylModel=axes(ss_SQ,'Position',[0.6 0.55 0.3 0.4 ]);title('Cylinder Model'  );
allaxes.SQ_Section=axes(ss_SQ,'Position',[0.8 0.05 0.15 0.25 ]);title('Section Check'  );

New_UIControl = copyobj(findobj(ss_HC,'Type','UIControl','-depth',1),ss_SQ) ;
allbtnO=allbtn;
for k=1:length(New_UIControl)
   if  strcmp(New_UIControl(k).Style,'pushbutton')
       currentDate = strcat('SQ_',New_UIControl(k).String);
       allbtn.(currentDate) = New_UIControl(k) ;    
   elseif strcmp(New_UIControl(k).Style,'edit')
       currentDate = strcat('SQ_',New_UIControl(k).UserData);
       alledit.(currentDate) = New_UIControl(k) ;    
   elseif strcmp(New_UIControl(k).Style,'popupmenu')
       currentDate = strcat('SQ_','SaveSection');
       allpop.(currentDate) = New_UIControl(k) ;    
   end    
end

allbg.SQ_Loop = uibuttongroup(ss_SQ,'Position',[0.5 0.22 0.1  0.15],'Title','Loop') ;
allbg.SQ_Unit = uibuttongroup(ss_SQ,'Position',[0.5 0.05 0.1  0.15],'Title','Unit') ;
allbg.SQ_Cylinder = uibuttongroup(ss_SQ,'Position',[0.62 0.05 0.13  0.32],'Title','Cylinder') ;

allrb.SQ_Loop_Outer = uicontrol('Style','radiobutton','Parent',allbg.SQ_Loop,...
    'Units','normalized','Position',[0.05 0.3 0.8 0.2],'String','Outer Loop'  );
allrb.SQ_Loop_Inner = uicontrol('Style','radiobutton','Parent',allbg.SQ_Loop,...
    'Units','normalized','Position',[0.05 0.8 0.8 0.2],'String','Inner Loop'  );

allrb.SQ_Unit_nm = uicontrol('Style','radiobutton','Parent',allbg.SQ_Unit,...
    'Units','normalized','Position',[0.05 0.3 0.8 0.2],'String','nm'  );
allrb.SQ_Unit_base = uicontrol('Style','radiobutton','Parent',allbg.SQ_Unit,...
    'Units','normalized','Position',[0.05 0.8 0.8 0.2],'String','base'  );

allrb.SQ_Cyl_Draw = uicontrol('Style','radiobutton','Parent',allbg.SQ_Cylinder,...
    'Units','normalized','Position',[0.05 0.8 0.8 0.2],'String','Draw Outline'  );
allrb.SQ_Cyl_Select = uicontrol('Style','radiobutton','Parent',allbg.SQ_Cylinder,...
    'Units','normalized','Position',[0.05 0.5 0.8 0.2],'String','Select Cylinder'  );
allrb.SQ_Cyl_AddOrDel = uicontrol('Style','radiobutton','Parent',allbg.SQ_Cylinder,...
    'Units','normalized','Position',[0.05 0.2 0.8 0.2],'String','Add or Delete'  );



TypeList={'Axes','UIControl'}; 
for Ti=1:length(TypeList)
    child_handles_Axes = findall(tgroup,'Type',TypeList{Ti});
    for k=1:length(child_handles_Axes)
    child_handles_Axes(k).FontSize=16 ;
    end
end
% handles.PB_Show3D.Parent=handles.tab1;
% handles.PB_Show3D.Position=[0.37 0.26+dy 0.1  0.05];



%------------













allaxes.HC_CylModel.ButtonDownFcn=@(src,evn)SetAxRotate(src,evn) ;
allaxes.HC_Extrude.ButtonDownFcn=@(src,evn)SetAxRotate(src,evn) ;
allaxes.HC_Section.ButtonDownFcn=@(src,evn)SetAxRotate(src,evn) ;
allaxes.SQ_CylModel.ButtonDownFcn=@(src,evn)SetAxRotate(src,evn) ;
allaxes.SQ_Extrude.ButtonDownFcn=@(src,evn)SetAxRotate(src,evn) ;
allaxes.SQ_Section.ButtonDownFcn=@(src,evn)SetAxRotate(src,evn) ;


% allaxes.HC_CylModel.Title.ButtonDownFcn=@(src,evn)SetAxRotateoff(src,evn) ;
fH=gcf; 
% fH.WindowButtonMotionFcn =@(src,evn)SetAxRotateoff(src,evn) ;
fH.KeyPressFcn=@(src,evn)testKeyPress(src,evn) ;
function SetAxRotate(src,evn)
% sdfsf=3
if evn.Button==1
    h=rotate3d; h.Enable = 'on'; 

end

end

function testKeyPress(src,evn)
evn;

end

function ssOH_Button(src,~)
ax= findobj(gcf,'Tag','MechOH3D');
axes(ax);
fH=gcf ;

fH.UserData.MovieStep =6 ; % play the 6th movie in ChangeAxesLimit

fH.KeyPressFcn=@(fH,evn)ChangeAxesLimit(fH,evn) ; % temporary change to other keypress fcn
fH.UserData.saveKeyMove3 = @(evn)ChangeAxesLimit(fH,evn) ; % temporary change to other keypress fcn

end
function ssSTEP_Button(src,~)
ax= findobj(gcf,'Tag','SS','Type','Axes');
axes(ax);
fH=gcf ;

src;
if isequal(src, findobj(gcf,'Tag','ss_Scaffold') )
    fH.UserData.MovieStep =4 ; % play the 4th movie in ChangeAxesLimit
elseif isequal(src, findobj(gcf,'Tag','ss_Staple') )
    fH.UserData.MovieStep =5 ; % play the 5th movie in ChangeAxesLimit
elseif isequal(src, findobj(gcf,'Tag','sss_Confforce') )
    fH.UserData.MovieStep =7 ; % play the 7th movie in ChangeAxesLimit
elseif isequal(src, findobj(gcf,'Tag','sss_Traj') )
    fH.UserData.MovieStep =8 ; % play the 8th movie in ChangeAxesLimit
end

fH.KeyPressFcn=@(fH,evn)ChangeAxesLimit(fH,evn) ; % temporary change to other keypress fcn
fH.UserData.saveKeyMove3 = @(evn)ChangeAxesLimit(fH,evn) ; % temporary change to other keypress fcn

end
function ssJSON_Button(src,~)
ax= findobj(gcf,'Tag','json3D','Type','Axes');
axes(ax);
fH=gcf ;

fH.KeyPressFcn=@(fH,evn)JsonKeyPress(fH,evn) ; % temporary change to other keypress fcn
fH.UserData.saveKeyMove4 = @(evn)JsonKeyPress(fH,evn) ; % temporary change to other keypress fcn

end
% function ssJSON_Button(src,~)
% ax= findobj(gcf,'Tag','json3D','Type','Axes');
% axes(ax);
% fH=gcf ;
% 
% fH.KeyPressFcn=@(fH,evn)ChangeAxesLimit(fH,evn) ; % temporary change to other keypress fcn
% fH.UserData.saveKeyMove3 = @(evn)ChangeAxesLimit(fH,evn) ; % temporary change to other keypress fcn
% 
% end

