function f= main()

%     <MagicDNA (Multi-component Assembly in a Graphical Interface guided by Computation for DNA origami) is a software for designing multi-component DNA origami structures.>
%     Copyright (C) <2020>  <Chao-Min Huang, Hai-Jun Su, and Carlos E. Castro>
%     The Ohio State University 
%     
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% main script, Create GUI


clear; clc ; % close all; 
filepath=mfilename('fullpath') ;
Inds = strfind(filepath,filesep) ;
cd( filepath(1: Inds(end)-1));


addpath(strcat(pwd, filesep, 'class' )) ;
addpath(strcat(pwd, filesep, 'functions' )) ;
addpath(strcat(pwd, filesep, 'parts' )) ;
addpath(strcat(pwd, filesep, 'scripts' )) ;
addpath(strcat(pwd, filesep, 'functions',filesep,'jsonlab' )) ;
addpath(strcat(pwd, filesep, 'InstructionMovie' )) ;

%get version to assign Tooltip or TooltipString 
% vv=version('-release') ;
% if str2num(vv(1:4))==2017
    ToolTipName='TooltipString' ;
% else
%     ToolTipName='Tooltip' ; 
% end

%% initialize
InitialGUI ;
h=rotate3d(gcf) ;
h.Enable ='off' ;
%% hide unnecessary toolbar
a = findall(gcf) ;
SStings={'Save Figure','New Figure','Open File' ,'Print Figure','Link Plot','Insert Colorbar','Insert Legend',...
    'Hide Plot Tools' ,'Show Plot Tools and Dock Figure'} ;
for k=1:length(SStings)
b = findall(a,'ToolTipString',SStings{k}) ;
set(b,'Visible','off')
end

%% -------Mechanism tab
%% -------STEP tab
allaxes.STEP=axes(ss_STEP,'Position',[0.05 0.1 0.6 0.85 ],'Tag','SS');title('');
allbtn.StartSTEP= uicontrol('Style','pushbutton','Parent',ss_STEP,...
    'Units','normalized','Position',[0.1 0.1 0.15  0.05],'String','StartSTEP','Tag','SS' ,'Visible' ,'off');
allbtn.StartSTEP.(ToolTipName)='Use a STEP file from CAD software to provide line models';

allbtn.StartPoints= uicontrol('Style','pushbutton','Parent',ss_STEP,...
    'Units','normalized','Position',[0.25 0.1 0.15  0.05],'String','StartPoints','Tag','SS','Visible' ,'off');
allbtn.StartPoints.(ToolTipName)='Provide a XYZ-point file or use Cartesian/cylindrical points to sketch lines.';

allbtn.FFCurve= uicontrol('Style','pushbutton','Parent',ss_STEP,...
    'Units','normalized','Position',[0.7 0.02 0.1 0.13],'String','FFCurve','Tag','SS');
% allbtn.FFCurve= uicontrol('Style','pushbutton','Parent',ss_STEP,...
%     'Units','normalized','Position',[0.4 0.1 0.15  0.05],'String','FFCurve','Tag','SS');
allbtn.FFCurve.(ToolTipName)='Spline sketch GUI or sketch on Cartesian/Cylindrical grids.';
% [0.85 0.02 0.1 0.13]

setGOfontsize( gctab , 18 , {'UIControl'} )  % set fontsize for uicontrol in this tab

sld_Font = uicontrol('Style', 'slider','Parent',ss_STEP,'Units','normalized',....
    'Min',4,'Max',18,'Value',8,'Position', [0.15 0.02 0.4  0.03] ,'Tag','SS');   
sld_Font.Callback=@(src,evn)set_CadDOM_FontSize(src,evn) ;
sld_Font.(ToolTipName) = 'Use this slider to change all the fontsizes of UI components in this figure.';
textFont = uicontrol('Style','text','Parent',ss_STEP,'Units','normalized','Position', [0.05 0.02 0.1  0.03],'String','Font Size' ,'Tag','SS');   
AssignIcon( allbtn.FFCurve,'SplineGUI.jpg' ) ;  

%% -------Assembly tab
allaxes.Assemblymain= axes(ss_Assembly,'Position',[0.05 0.07 0.6 0.88 ],'Tag','AssemblyMain' ); hold on;
allaxes.AssemblymainHiden = axes(ss_Assembly,'Position',[0.05 0.07 0.6 0.88 ],'Visible','off','hittest','off','Tag','HidenAssemblyMain'); % Invisible axes
hlink =linkprop([allaxes.Assemblymain allaxes.AssemblymainHiden],{ 'XLim' 'YLim' 'ZLim' ,'PlotBoxAspectRatio','View' }); % The axes should stay aligned
allaxes.Assemblymain.UserData.hlink= hlink ;
allaxes.AssemblymGraph= axes(ss_Assembly,'Position',[0.7,0.3,0.23,0.23],'Tag','AssemblyGraph');  hold on;
setAllowAxesRotate(h,allaxes.AssemblymGraph,0) ;
%% ------Parameter tab
% % Parameter = uibuttongroup(ss_RoutingParameter,'Position',[0.6 0.7 0.15 0.2],'Title','Prevent stacking','FontSize',12,'BorderWidth',2) ;
% 
% Parameter.scafXoverTol = uibuttongroup(ss_RoutingParameter,'Position',[0.6 0.55 0.35 0.15],'Title','Ignore scaf Xovers from ends by ','FontSize',12,'BorderWidth',2) ;
% Parameter.stapXoverTol = uibuttongroup(ss_RoutingParameter,'Position',[0.6 0.40 0.35 0.15],'Title','Ignore staple Xovers from ends by ','FontSize',12,'BorderWidth',2) ;
% Parameter.stapCutLimit = uibuttongroup(ss_RoutingParameter,'Position',[0.6 0.25 0.35 0.15],'Title','staple break limit','FontSize',12,'BorderWidth',2) ;
% 
% Parameter.choice.scafXoverTol.text=uicontrol('Style','text','Parent',Parameter.scafXoverTol,'Units','normalized','Position',[0.02 0.65 0.8 0.3],'String','Bundle index, separate by space or '',''',...
%     'HorizontalAlignment','left');
% Parameter.choice.scafXoverTol.edit=uicontrol('Style','edit','Parent',Parameter.scafXoverTol,'Units','normalized','Position',[0.05 0.2 0.8 0.2]  );


%% -------Scaf tab


allaxes.MechScaffold3D= axes(ss_Scaffold,'Position',[0.05 0.05 0.5 0.85],'Tag','MechScaffold3D');  hold on; title('Scaffold 3D'  );
%---------add legend as instructions, single row 
            ForLegend=surface(nan, nan,'Tag','test','EdgeColor','interp','FaceColor','none');
            % ForLegend=surface(nan, nan, 'Linestyle', 'none', 'Marker', 'none' );
%                         ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');

            hLg= legend(ForLegend,'Click me for instructions','Location','northwest' ) ; hLg.String={' Scaffold routing'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
%             ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_Scaf( src,evn,allaxes.MechScaffold3D );
            hLg.Title.String='Click me for instructions' ;
            hLg.Units='normalized'; %hLg.AutoUpdate ='off';
            hLg.Position=[0.0063 0.9728 0.1569 0.025];
%------------------



allaxes.MechScaffold2D= axes(ss_Scaffold,'Position',[0.6 0.35 0.35 0.55],'Tag','MechScaffold2D');  hold on; title('Scaffold 2D'  );
allaxes.MechScaffold2D.Visible='off';

ScafOption.prevStack = uibuttongroup(ss_Scaffold,'Position',[0.6 0.7 0.15 0.2],'Title','Prevent stacking','FontSize',12,'BorderWidth',2) ;
ScafOption.ScafXover = uibuttongroup(ss_Scaffold,'Position',[0.8 0.7 0.15 0.2],'Title','Scaffold Crossovers','FontSize',12,'BorderWidth',2) ;

ScafOption.BundleNoScafLoop = uibuttongroup(ss_Scaffold,'Position',[0.6 0.55 0.35 0.15],'Title','Bundles No Scaf Loop','FontSize',12,'BorderWidth',2) ;
ScafOption.PolyBundle = uibuttongroup(ss_Scaffold,'Position',[0.6 0.40 0.35 0.15],'Title','Poly Bundles','FontSize',12,'BorderWidth',2) ;
ScafOption.NoScafXover = uibuttongroup(ss_Scaffold,'Position',[0.6 0.25 0.35 0.15],'Title','No Scaf Xover','FontSize',12,'BorderWidth',2) ;
ScafOption.minDist_btwXover = uibuttongroup(ss_Scaffold,'Position',[0.6 0.15 0.35 0.1],'Title','Min distances between Xover','FontSize',12,'BorderWidth',2) ;

% allrb.HC_Loop_Outer = uicontrol('Style','radiobutton','Parent',allbg.HC_Loop,...
%     'Units','normalized','Position',[0.05 0.3 0.8 0.2],'String','Outer Loop'  );
ScafOption.choice.prevStack.scafloop=uicontrol('Style','radiobutton','Parent',ScafOption.prevStack,'Units','normalized','Position',[0.05 0.6 0.8 0.2],'String','scaf loop'  );
ScafOption.choice.prevStack.polyT=uicontrol('Style','radiobutton','Parent',ScafOption.prevStack,'Units','normalized','Position',[0.05 0.3 0.8 0.2],'String','poly T'  );
ScafOption.choice.prevStack.scafloop.(ToolTipName)='Leave single-stranded scaffold loops on the ends to prevent aggregation due to base-stacking.';
ScafOption.choice.prevStack.polyT.(ToolTipName)='Shrink the scaffold domain to cross-over locations on both end. Entend the staples with ssDNA value specified in ssDNA as single-stranded portion to prevent base-stacking. ';


ScafOption.choice.scafXover.Min=uicontrol('Style','radiobutton','Parent',ScafOption.ScafXover,'Units','normalized','Position',[0.05 0.6 0.8 0.2],'String','minimize'  );
ScafOption.choice.scafXover.Max=uicontrol('Style','radiobutton','Parent',ScafOption.ScafXover,'Units','normalized','Position',[0.05 0.3 0.8 0.2],'String','maximize'  );

ScafOption.choice.BundleNoScafLoop.text=uicontrol('Style','text','Parent',ScafOption.BundleNoScafLoop,'Units','normalized','Position',[0.02 0.65 0.8 0.3],'String','Bundle index, separate by space or '',''',...
    'HorizontalAlignment','left');
ScafOption.choice.BundleNoScafLoop.edit=uicontrol('Style','edit','Parent',ScafOption.BundleNoScafLoop,'Units','normalized','Position',[0.05 0.2 0.8 0.2]  );
ScafOption.choice.BundleNoScafLoop.text.(ToolTipName)='Specify bundles without scaffold loops.';
ScafOption.choice.BundleNoScafLoop.edit.(ToolTipName)='Specify bundles without scaffold loops.' ;


ScafOption.choice.PolyBundle.text=uicontrol('Style','text','Parent',ScafOption.PolyBundle,'Units','normalized','Position',[0.02 0.65 0.8 0.3],'String','Bundle index, separate by space or '',''',...
    'HorizontalAlignment','left');
ScafOption.choice.PolyBundle.edit=uicontrol('Style','edit','Parent',ScafOption.PolyBundle,'Units','normalized','Position',[0.05 0.2 0.8 0.2]  );
ScafOption.choice.PolyBundle.text.(ToolTipName)='Extend ssDNA scaffold loops to cross-over locations for polymerizing the structures with extra staples.';
ScafOption.choice.PolyBundle.edit.(ToolTipName)='Extend ssDNA scaffold loops to cross-over locations for polymerizing the structures with extra staples.';



ScafOption.choice.NoScafXover.text=uicontrol('Style','text','Parent',ScafOption.NoScafXover,'Units','normalized','Position',[0.02 0.65 0.8 0.3],'String','Bundle index, separate by space or '',''',...
    'HorizontalAlignment','left');
ScafOption.choice.NoScafXover.edit=uicontrol('Style','edit','Parent',ScafOption.NoScafXover,'Units','normalized','Position',[0.05 0.2 0.8 0.2]  );
ScafOption.choice.NoScafXover.text.(ToolTipName)='No internal scaffold cross-overs in these bundles for reconfiguration or linear actuations. May reduce the algorithm robustness in some case.';
ScafOption.choice.NoScafXover.edit.(ToolTipName)='No internal scaffold cross-overs in these bundles for reconfiguration or linear actuations. May reduce the algorithm robustness in some case.';


nums= 8:2:14;
ScafOption.choice.minDist_btwXover.popup=uicontrol('Style','popup','Parent',ScafOption.minDist_btwXover,'Units','normalized','Position',[0.02 0.65 0.4 0.3],'String',num2cell(nums),...
    'HorizontalAlignment','left');
ScafOption.choice.minDist_btwXover.popup.(ToolTipName)='The scaffold algorithm try to avoid the closeness between scaffold cross-overs. Both successful and failed results still get a routing but indicate in the command line. ';
         
allbtn.MechScafAgain= uicontrol(ss_Scaffold,'Style','pushbutton',...
    'Units','normalized','Position',[0.6 0.05 0.1 0.1],'String','GetScaffoldRouting','Tag','ScafAgain');
allbtn.MechScafDefault= uicontrol(ss_Scaffold,'Style','pushbutton',...
    'Units','normalized','Position',[0.85 0.05 0.1 0.1],'String','Default','Tag','DefaultScaf');
allbtn.MechScafDefault.Callback=@(src,evn)SetScafOptionDefault(src,evn,ScafOption)   ;
ss_Scaffold.UserData.ScafOption=ScafOption; %  move handle to userdata since this tab has assign Tag.
allbtn.PlotScaf= uicontrol(ss_Scaffold,'Style','pushbutton',...
    'Units','normalized','Position',[0.7 0.05 0.1 0.1],'String','PlotScaf','Tag','PlotScaf');
allbtn.PlotScaf.Callback=@(src,evn)PlotScaf_Fcn(src,evn)   ;

allbtn.MultiScaf= uicontrol(ss_Scaffold,'Style','pushbutton',...
    'Units','normalized','Position',[0.8 0.7 0.1 0.1],'String','MultiScaf','Tag','MultiScaf');
allbtn.MultiScaf.Callback=@(src,evn)SplitScafToMultiScaf(src,evn)   ;


allbtn.MechScafAgain.(ToolTipName)='Apply the scaffold algorithm to obtain a scaffold routing.' ;
allbtn.MechScafDefault.(ToolTipName)='Set the scaffold options to default values' ;
allbtn.PlotScaf.(ToolTipName)='If users have used ''caDNAno'' to update the routing, this button is to compare the routings.' ;

align([allbtn.MechScafAgain allbtn.PlotScaf allbtn.MechScafDefault],'distribute','bottom');


setGOfontsize( gctab , 16 , {'UIControl'} )  % set fontsize for uicontrol in this tab

ScafOption.ScafXover.Visible ='off'; % hide the new scaffold algorithm option

% setGOfontsize( gctab , 14 , {'ButtonGroup'} )  % set fontsize for uicontrol in this tab

%% -------Stap tab

allaxes.MechStaple3D= axes(ss_Staple,'Position',[0.05 0.05 0.5 0.85],'Tag','MechStaple3D');  hold on; title('Staple 3D'  );
%---------add legend as instructions, single row 
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            hLg= legend(ForLegend,'Click me for instructions','Location','northwest' ) ; hLg.String={'\bf{Click me for instructions}'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_Stap( src,evn,allaxes.MechStaple3D );
            hLg.Units='normalized'; hLg.AutoUpdate ='off'; hLg.Position=[0.0063 0.9728 0.1569 0.025];
%------------------
allaxes.MechStaple2D= axes(ss_Staple,'Position',[0.6 0.45 0.35 0.5],'Tag','MechStaple2D');  hold on; title('Staple 2D'  );
allbtn.MechStapleAgain= uicontrol(ss_Staple,'Style','pushbutton',...
    'Units','normalized','Position',[0.6 0.05 0.2 0.1],'String','GetStaple','Tag','StapleAgain');
allbtn.MechStapleJson= uicontrol(ss_Staple,'Style','pushbutton',...
    'Units','normalized','Position',[0.85 0.05 0.1 0.1],'String','Use cadnano','Tag','StapleJson');

allbtn.MechStapleAgain.(ToolTipName)='Apply the staple algorithm to obtain staple routings.' ;
allbtn.MechStapleJson.(ToolTipName)='Advanced use for manually modifying scaffold and staple routing in caDNAno software.' ;



StapOption.type = uibuttongroup(ss_Staple,'Position',[0.6 0.2 0.28 0.2],'Title','Type','FontSize',12,'BorderWidth',2) ;
StapOption.choice.type.straightuncut=uicontrol('Style','radiobutton','Parent',StapOption.type,'Units','normalized','Position',[0.05 0.7 0.4 0.2],'String','straight uncut'  );
StapOption.choice.type.zigzaguncut=uicontrol('Style','radiobutton','Parent',StapOption.type,'Units','normalized','Position',[0.05 0.2 0.4 0.2],'String','zigzag uncut'  );
StapOption.choice.type.zigzagcut=uicontrol('Style','radiobutton','Parent',StapOption.type,'Units','normalized','Position',[0.55 0.2 0.4 0.2],'String','zigzag cut'  );
StapOption.choice.type.straightcut=uicontrol('Style','radiobutton','Parent',StapOption.type,'Units','normalized','Position',[0.55 0.7 0.4 0.2],'String','straight cut'  );
StapOption.choice.type.zigzagcut.Value =1;

StapOption.choice.type.straightuncut.(ToolTipName)='Leave the staple stands without cross-overs and without breaking into short strands.' ;
StapOption.choice.type.zigzaguncut.(ToolTipName)='Apply all cross-overs to the staples without breaking into short strands.' ;
StapOption.choice.type.zigzagcut.(ToolTipName)='Apply all cross-overs to the staples and breaking into short strands.' ;
StapOption.choice.type.straightcut.(ToolTipName)='Leave the staple stands without cross-overs and  breaking into short strands.' ;



StapOption.halfXover = uibuttongroup(ss_Staple,'Position',[0.9 0.2 0.05 0.2],'Title','halfXover','FontSize',12,'BorderWidth',2) ;
StapOption.choice.halfXover.yes=uicontrol('Style','radiobutton','Parent',StapOption.halfXover,'Units','normalized','Position',[0.05 0.7 0.8 0.2],'String','Yes'  );
StapOption.choice.halfXover.no=uicontrol('Style','radiobutton','Parent',StapOption.halfXover,'Units','normalized','Position',[0.05 0.3 0.8 0.2],'String','No'  );
StapOption.choice.halfXover.yes.(ToolTipName)='Apply half cross-overs to both ends to reduce fraying if the orientation allows.' ;
StapOption.choice.halfXover.no.(ToolTipName)='No half cross-overs.' ;


ss_Staple.UserData.StapOption=StapOption; %  move handle to userdata since this tab has assign Tag.

setGOfontsize( gctab , 16 , {'UIControl'} )  % set fontsize for uicontrol in this tab
setAllowAxesRotate(h,allaxes.MechStaple2D,0) ;

%% -------overhang tab

allaxes.MechOH3D= axes(ss_OH,'Position',[0.05 0.05 0.63 0.85],'Tag','MechOH3D');  hold on; title('Over-hang design'  );
%---------add legend as instructions, single row 
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            hLg= legend(ForLegend,'Click me for instructions','Location','northwest' ) ; hLg.String={'\bf{Click me for instructions}'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_OH( src,evn,allaxes.MechOH3D);
            hLg.Units='normalized'; hLg.AutoUpdate ='off'; hLg.Position=[0.0063 0.9728 0.1569 0.025];
%------------------
allaxes.MechOH2D= axes(ss_OH,'Position',[0.7 0.78 0.28 0.2],'Tag','MechOH2D');  hold on; %title('OH 2D'  );
% imshow([pwd filesep 'OH.jpg']) ;axis off ;


allbtn.MechOH= uicontrol(ss_OH,'Style','pushbutton',...
    'Units','normalized','Position',[0.7 0.1 0.06 0.1],'String','OHfunction','Tag','OHfunction');
allbtn.MechOH2= uicontrol(ss_OH,'Style','pushbutton',...
    'Units','normalized','Position',[0.78 0.1 0.06 0.1],'String','Add','Tag','OHAdd','Enable','off');
allbtn.MechOH3= uicontrol(ss_OH,'Style','pushbutton',...
    'Units','normalized','Position',[0.86 0.1 0.06 0.1],'String','Apply 1st to All','Tag','ApplyAll','Enable','off');
allbtn.MechOH4= uicontrol(ss_OH,'Style','pushbutton',...
    'Units','normalized','Position',[0.92 0.1 0.06 0.1],'String','Add Single','Tag','Add Single','Enable','off');
align([allbtn.MechOH allbtn.MechOH2 allbtn.MechOH3 allbtn.MechOH4 ],'distribute','bottom');



allbtn.MechOH.(ToolTipName)='Initialize the tool or Clear the table data.' ;
allbtn.MechOH2.(ToolTipName)='Add a pair of the red and blue nodes to the table.' ;
allbtn.MechOH3.(ToolTipName)='For assigning the same setting to all overhang pairs. Sequences can still remain the same or assign differently later in caDNAno tab.' ;


OHtext1 = uicontrol(ss_OH,'Style','text','Units','normalized','Position',[0.7 0.05 0.1 0.02],'String','Text1','Tag','OHText1','HorizontalAlignment','left');
OHtext2 = uicontrol(ss_OH,'Style','text','Units','normalized','Position',[0.85 0.05 0.1 0.02],'String','Text2','Tag','OHText2','HorizontalAlignment','left');

str2 = cellfun(@num2str,num2cell([0, 6:12, 14 , 16 , 18, 20]),'uniformoutput',0);
str= {'3''','5'''}; str3= {'Connected','Free end' ,'Double crossover', 'Single crossover', 'Double overhang1', 'Double overhang2' };
t = uitable(ss_OH,'Units','normalized','Position',[ 0.7000    0.2500    0.2800    0.5800],'ColumnFormat',({[]  str [] [] str str2 str2 str3}),...
       'ColumnEditable', false,'ColumnName',{'A','B','C','D','E','F','G','H','I'},'Visible','off','Tag','OHTable');  
TStrO='<html><font size=5>Location</h1></html>' ; 
% t.ColumnEditable=[false false true false false] ;
 t.ColumnName{1}=replace(TStrO,"Location","P") ;
 t.ColumnName{2}=replace(TStrO,"Location","Ends") ;
 t.ColumnName{3}=replace(TStrO,"Location","Enable") ;
 t.ColumnName{4}=replace(TStrO,"Location","P") ;
 t.ColumnName{5}=replace(TStrO,"Location","Ends") ;
 t.ColumnName{6}=replace(TStrO,"Location","nA") ;
 t.ColumnName{7}=replace(TStrO,"Location","nB") ;
 t.ColumnName{8}=replace(TStrO,"Location","Ax on Closing S.") ;
 t.ColumnName{9}=replace(TStrO,"Location","Highlight") ;
t.ColumnEditable=[false true true false true true true true true] ; t.FontSize=12 ;

% t2 = uitable(ss_OH,'Units','normalized','Position',[ 0.7000    0.2200    0.2800    0.2],'ColumnFormat',({[]  str [] [] str str2 str2 str3}),...
%        'ColumnEditable', false,'ColumnName',{'A','B','C','D','E','F','G','H','I'},'Visible','off','Tag','OHTable2');  




TooltipStr = 'This table is used to specify the overhang options for each pair of overhangs.';
TooltipStr = [TooltipStr newline 'The overhangs will be extended from 5'' or 3'' site on PAx and PBx points.'] ;
TooltipStr = [TooltipStr newline 'Enable column allows to temporarily ignored in the algorith and easily recover if needed.'] ;
TooltipStr = [TooltipStr newline 'Depend on the selected primes, it had multiple combinations to design the corresponding closing strand(See the schematic).'] ;
TooltipStr = [TooltipStr newline 'See UI operations in the instructions.'] ;

% set(t,'Tooltip', TooltipStr) ;
% t.Tooltip=TooltipStr;

%%
%------------- cadnano tab
allaxes.cadnano3D= axes(ss_json,'Position',[0.05 0.35 0.4 0.6],'Tag','json3D');  %hold on; title('OH 3D'  );
%---------add legend as instructions, single row 
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            hLg= legend(ForLegend,'Click me for instructions','Location','northwest' ) ; hLg.String={'\bf{Click me for instructions}'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_cadnano( src,evn,allaxes.cadnano3D);
            hLg.Units='normalized'; hLg.AutoUpdate ='off'; hLg.Position=[0.0063 0.9728 0.1569 0.025];
%------------------
allaxes.cadnano2D= axes(ss_json,'Position',[0.5 0.35 0.4 0.6],'Tag','json2D');  %hold on; title('OH 3D'  );

t_json = uitable(ss_json,'Units','normalized','Position',[0.05 0.05 0.5 0.25],'ColumnFormat',({[]  str [] [] str [] }),...
       'ColumnEditable', false,'ColumnName',{'Start','End','Sequence','Length','Color','Note'},'Visible','on','Tag','t_json');  
allbtn.json1= uicontrol(ss_json,'Style','pushbutton',...
    'Units','normalized','Position',[0.58 0.05 0.08 0.1],'String','Preview','Tag','btn_json1');
allbtn.json2= uicontrol(ss_json,'Style','pushbutton',...
    'Units','normalized','Position',[0.78 0.05 0.08 0.1],'String','ExportJSON','Tag','btn_json2');
allbtn.json3= uicontrol(ss_json,'Style','pushbutton',...
    'Units','normalized','Position',[0.88 0.05 0.08 0.1],'String','ExportCSV','Tag','btn_json3');
allbtn.json4= uicontrol(ss_json,'Style','pushbutton',...
    'Units','normalized','Position',[0.68 0.05 0.08 0.1],'String','OverHangOption','Tag','btn_json4','Enable','off');
allbtn.json1.(ToolTipName)='Preview the scaffold and staple routings in both 2D and 3D panels. Assign the sequence of the scaffold to obtain complementary sequences on staples.' ;
allbtn.json4.(ToolTipName)='If users have used the Overhagn design tool, assign the sequencess on the closing strands to get complementary staple overhang sequences.' ;
allbtn.json2.(ToolTipName)='Export the routing as .json for caDNAno software. Two files will be generated but actually the same routing due to hybrid-lattice design.' ;
allbtn.json3.(ToolTipName)='Export the data in the table for staple sequeces list. The indexes of cylinders have been converted as caDNAno software.' ;

t_json.(ToolTipName)='The staple sequence list and its corresponding positions in the caDNAno software. Click on a staple to see the 2D and 3D locations.' ;


t_json.ColumnEditable=[false false false false false  true ] ;
t_json.FontSize=10;

jsonSlider1=uicontrol(ss_json,'Style','slider',...
    'Units','normalized','Position',[0.60 0.2 0.12 0.04],'Tag','slider_json1');
jsonSlider2=uicontrol(ss_json,'Style','slider',...
    'Units','normalized','Position',[0.75 0.2 0.12 0.04],'Tag','slider_json2');
jsonSlider1.(ToolTipName)='Change the transparency of scaffold strands in both 2D and 3D panels' ;
jsonSlider2.(ToolTipName)='Change the transparency of staple strands in both 2D and 3D panels' ;



txtjson1 = uicontrol(ss_json,'Style','text','Unit','normalized','FontSize',14,'Position', [0.60 0.26 0.12 0.05],....   
        'String','Scaffold Transparency','Enable','Inactive','HorizontalAlignment','left','Tag','JsonTxt_scaf');   
txtjson1.(ToolTipName)='Click me to change scaffold representation between helical and chicken-wire.' ;      
txtjson2 = uicontrol(ss_json,'Style','text','Unit','normalized','FontSize',14,'Position', [0.75 0.26 0.12 0.05],....   
        'String','Staple Transparency','Enable','Inactive' ,'HorizontalAlignment','left','Tag','JsonTxt_stap'); 
 txtjson2.(ToolTipName)='Click me to change staple representation between helical and chicken-wire.' ;      
   
    str ={'Click route on 2D panel','keyboard to move scaffold base ' ,'right click to assign step ',...
        ' ? -> overhang seq', ' * -> exceed scaffold length' } ;
txtjson3 = uicontrol(ss_json,'Style','text','Unit','normalized','FontSize',12,'Position', [0.87 0.26 0.12 0.06],....   
        'String',str,'Enable','Inactive','HorizontalAlignment','left','Tag','JSON_Text3' );   
txtjson3.ButtonDownFcn=@(src,evn)uisetfont(txtjson3);
txtjson3.(ToolTipName)='Use to see where the selected base locates on the caDNAno software.' ;

% mytxt.ButtonDownFcn = 'disp(''Text was clicked'')';

strs={'Show num', 'Show ends', 'Both ','None' } ;
jsonPop=uicontrol('Style','popupmenu','Parent',ss_json,...
    'Units','normalized','Position',[0.9 0.2  0.05  0.04],'Value',3,'String',strs);
jsonPop.(ToolTipName)='Change representations for staple strands.' ;

jsonPop2=uicontrol('Style','popupmenu','Parent',ss_json,...
    'Units','normalized','Position',[0.9 0.15  0.05  0.04],'Value',1,'String',{'select scaffold'});
jsonPop2.(ToolTipName)='Change representations for staple strands.' ;

setAllowAxesRotate(h,allaxes.cadnano2D,0) ;

%% oxDNA tab 
allaxes.oxDNA= axes(sss_Confforce,'Position',[0.05 0.25 0.5 0.7],'Tag','oxDNA');  %hold on; 
%----------------
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            hLg= legend(ForLegend,'Click me for instructions','Location','northwest' ) ; hLg.String={'\bf{Click me for instructions}'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_oxDNA1( src,evn,allaxes.oxDNA);
            hLg.Units='normalized'; hLg.AutoUpdate ='off'; hLg.Position=[0.0063 0.9728 0.1569 0.025];


allbtn.oxDNA1= uicontrol(sss_Confforce,'Style','pushbutton',...
    'Units','normalized','Position',[0.6 0.05 0.08 0.1],'String','Initial oxDNA','Tag','btn_oxDNA1');
allbtn.oxDNA2= uicontrol(sss_Confforce,'Style','pushbutton',...
    'Units','normalized','Position',[0.7 0.05 0.08 0.1],'String','Export conf+top','Tag','btn_oxDNA2');
allbtn.oxDNA3= uicontrol(sss_Confforce,'Style','pushbutton',...
    'Units','normalized','Position',[0.8 0.05 0.08 0.1],'String','Export BUILD','Tag','btn_oxDNA3');
allbtn.oxDNA1.(ToolTipName)='Initialize the configuration in the plot. Later allow users can adjust the configuration with the slider.' ;
allbtn.oxDNA2.(ToolTipName)='Export the topology and the configuration file. Automatically open another UI for fine tuning with rigid-body-transformation.' ;
allbtn.oxDNA3.(ToolTipName)='Export the current configuration to .bild file for Chimera for better quality.' ;



oxDNASlider=uicontrol(sss_Confforce,'Style','slider',...
    'Units','normalized','Position',[0.05 0.05 0.5 0.05],'Tag','slider_oxdna');
oxDNASlider.Max =4 ;


oxDNACheck=uicontrol(sss_Confforce,'Style','checkbox','String','Include Closing strand',...
    'Units','normalized','Position',[0.7 0.15 0.2 0.1],'Tag','check_oxdna');

% sss_pattern--------------
% allaxes.oxDNA2Pat= axes(sss_pattern,'Position',[0.05 0.25 0.5 0.7],'Tag','oxDNA2Pat');  %hold on; 
% allbtn.oxDNAPat1= uicontrol(sss_pattern,'Style','pushbutton',...
%     'Units','normalized','Position',[0.6 0.05 0.08 0.1],'String','Initial oxDNAPattern','Tag','btn_oxDNAPat1');
% allbtn.oxDNAPat2= uicontrol(sss_pattern,'Style','pushbutton',...
%     'Units','normalized','Position',[0.7 0.05 0.08 0.1],'String','Export conf+top','Tag','btn_oxDNAPat2');
% allbtn.oxDNAPat3= uicontrol(sss_pattern,'Style','pushbutton',...
%     'Units','normalized','Position',[0.8 0.05 0.08 0.1],'String','Export BUILD','Tag','btn_oxDNAPat3');
% str2 = cellfun(@num2str,num2cell(1:3),'uniformoutput',0);
% popOxdnaPattern   = uicontrol(sss_pattern,'Style', 'listbox',...
%     'String', str2,'Unit','normalized','Position', [0.6 0.2 0.08 0.1],'Tag','popOxdnaPattern'); 
% checkH_oxDNAPat = uicontrol(sss_pattern,'Style', 'checkbox','String', 'Trans/Rotate','Unit','normalized','Position', [0.7 0.2 0.08 0.1],'Tag','checkH_oxDNAPat'  ); 
% editH_oxDNAPat = uicontrol(sss_pattern,'Style', 'edit','String', '1','Unit','normalized','Position', [0.9 0.2 0.08 0.1],'Tag','editH_oxDNAPat'   ); 
% %-----------------------

% sss_Traj
allaxes.oxDNATrajAxe1= axes(sss_Traj,'Position',[0.05 0.25 0.5 0.7],'Tag','oxDNATrajAxe1');  %hold on; 
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            hLg= legend(ForLegend,'Click me for instructions','Location','northwest' ) ; hLg.String={'\bf{Click me for instructions}'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_oxDNA2( src,evn,allaxes.oxDNATrajAxe1);
            hLg.Units='normalized'; hLg.AutoUpdate ='off'; hLg.Position=[0.0063 0.9728 0.1569 0.025];



allaxes.oxDNATrajAxe2= axes(sss_Traj,'Position',[0.6 0.45 0.3 0.5],'Tag','oxDNATrajAxe2');  %hold on; 

allbtn.oxDNATraj1= uicontrol(sss_Traj,'Style','pushbutton',...
    'Units','normalized','Position',[0.6 0.05 0.08 0.1],'String','Load Initail+Traj','Tag','btn_oxDNATraj1');
allbtn.oxDNATraj2= uicontrol(sss_Traj,'Style','pushbutton',...
    'Units','normalized','Position',[0.7 0.05 0.08 0.1],'String','Compute deform','Tag','btn_oxDNATraj2','Visible','off');
allbtn.oxDNATraj3= uicontrol(sss_Traj,'Style','pushbutton',...
    'Units','normalized','Position',[0.8 0.05 0.08 0.1],'String','Export BUILD','Tag','btn_oxDNATraj3');
allbtn.oxDNATraj4= uicontrol(sss_Traj,'Style','pushbutton',...
    'Units','normalized','Position',[0.8 0.2 0.08 0.1],'String','RMSD RMSF','Tag','btn_oxDNATraj4');
allbtn.oxDNATraj5= uicontrol(sss_Traj,'Style','pushbutton',...
    'Units','normalized','Position',[0.9 0.2 0.08 0.1],'String','Export RMSF','Tag','btn_oxDNATraj5');

allbtn.oxDNATraj1.(ToolTipName)='Load the topology, initial configuration, and trajectory file to visualize the result.' ;
allbtn.oxDNATraj3.(ToolTipName)='Export the selected configuration as the popup menu to .bild file for Chimera.' ;
allbtn.oxDNATraj4.(ToolTipName)='Analyze the trajectory with RMSD and RMSF calculation. Also compute the average configuration and color-coded with RMSF values.' ;
allbtn.oxDNATraj5.(ToolTipName)='Export the average configuration with colors as RMSF value to .bild file for Chimera.' ;



popOxdnaTraj   = uicontrol(sss_Traj,'Style', 'popup',...
    'String', 'Frame','Unit','normalized','Position', [0.6 0.2 0.08 0.1],'Tag','popOxdnaTraj'); 
popOxdnaTraj.(ToolTipName)='Select the desired configuration to visualize. ''0'' means the initial configuration.  ' ;
listOxdnaTraj   = uicontrol(sss_Traj,'Style', 'listbox',...
    'String', 'Frame','Unit','normalized','Position', [0.7 0.2 0.08 0.1],'Tag','listOxdnaTraj'); 

listOxdnaTraj.(ToolTipName)='Select the desired configurations to visualize. Multiple selection is possible by ''Ctrl'' and ''Shift''. ''0'' means the initial configuration.  ' ;



allGraphics.allaxes=allaxes;
allGraphics.allbtn=allbtn;
allGraphics.allrb=allrb;
allGraphics.alledit=alledit;
allGraphics.allpop=allpop;

% UpdateConn


%% Setup Functions

SQobject= partSQ(allGraphics) ;
HCobject= partHC(allGraphics) ;

 mc= ?hyperbundle ;
%-------
% readSTEPFile(allGraphics,ss_STEP)
allbtn.StartSTEP.Callback=@(src,evn)readSTEPFile(src,evn,allGraphics,ss_STEP ,allbtn) ;
allbtn.StartPoints.Callback=@(src,evn)readSTEPFile(src,evn,allGraphics,ss_STEP ,allbtn) ;
allbtn.FFCurve.Callback=@(src,evn)FFCurve(src,evn,allGraphics,ss_STEP ,allbtn) ;

% allbtn.MechScafAgain.Callback= @(src,evn)SearchScaf_every13(src,evn) ;

allbtn.MechScafAgain.Callback= @(src,evn)SearchScaf(src,evn) ;

allbtn.MechStapleAgain.Callback= @(src,evn)SearchStap(src,evn) ;
allbtn.MechStapleJson.Callback= @(src,evn)UseCadnano(src,evn) ;
allbtn.MechOH.Callback= @(src,evn)OverhangInitial(src,evn) ;

allbtn.json1.Callback= @(src,evn)cadnanoInitial(src,evn,jsonSlider1,jsonSlider2,jsonPop,jsonPop2) ;
allbtn.oxDNA1.Callback= @(src,evn)oxDNAInitial(src,evn) ;
allbtn.oxDNAPat1.Callback= @(src,evn)oxDNAPatInitial(src,evn) ;
allbtn.oxDNATraj1.Callback= @(src,evn)LoadoxDNATraj(src,evn) ;



XY = GiveHCinPlaneP( 1 ) ;    nCyl= 6 ; 
AddBundle=  BundleCylinderHC(1,[],50*ones(1,6),200*ones(1,6),XY) ;

DefaultHB=hyperbundle(1,AddBundle);   % never use latter. Only for accessing the function to initialize

DefaultHBPath=[ pwd filesep 'Previous Mechanism' filesep 'Hinge.mat'];   % default hyper bundle
% DefaultHBPath=[ pwd '\Previous Mechanism\Hinge.mat'];   % default hyper bundle



S=open(DefaultHBPath) ; DefHB=S.S ;
DefaultHB=DefaultHB.loadobj(DefHB) ;            
btn_IntLoad = uicontrol(ss_Assembly,'Style', 'pushbutton', 'String', 'Initialize','Unit','normalized', 'Position', [0.87 0.55 0.05 0.05] ,...
            'Callback', {@(src,evn)loadIniMech( src,evn,DefaultHB )} );  
btn_IntLoad.UserData.keepme = 1 ;
        
f.SizeChangedFcn=@(src,evn)main_ChangeSize(src,evn) ;
f.UserData.MovieStep =0 ;

end

        


