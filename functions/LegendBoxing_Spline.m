function  LegendBoxing_Spline( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
       if evn.IntersectionPoint(2)>0
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';%CreateStruct.Units = 'normalized'; %CreateStruct.Position = 'normalized';
           STR={'\fontsize{11}\bf{Spline Sketch GUI}\rm:'; 
               'This GUI was newly added to enhance the design capability of freeform structure in MagicDNA. ';
               'Spline(s) could be loaded as XYZ file with one header line which specifying the initial nodes.';
               'These nodes are called break points(nodes) for constructing spline objects.  ';
               'In addition, moving these nodes allow users to control the freeform shape of the spline by using the following GUI interactions.';
               'Once finished, export the spline models (open or closed loop ) to the main sketch tab for further extrude or sweep operations.'
               '';
               'If there are multiple splines, select the spline of interest to edit using the popupmenu in the right.'
               'Multiple operations/modes on the spline are available in the other popupmenu.'
               '';
               'Mouse :' ;
               'Left Click -> Select one node.';
               'Middle Click -> Select a series of nodes between the previous to the current node.';
               'Right Click -> Add the one as selected.';
               '';
               'Align / Projection: '
               'For the selected spline, project all nodes to XY/XZ/YZ plane if select one or no node.'
               'If selected multiple nodes, align these nodes along the directions, i.e. assigning the same Z-value for the XY button.'
               '';
               'UI functions: ';
               '[m] : Move the center of all nodes to (0 0 0).'
               '[t] : print XYZ coordinate of the break nodes on the command line.'
               '[b] : In \bf{Break}\rm mode, break the spline into two lines on the selected node.'
               '[return]: In \bf{Insert}\rm mode, select a point on spline (under spline representation, not straight) and click ''return'' to add a break node. '
               '[backspace]: In \bf{Delete}\rm mode, delete the selected nodes. '
               '[l] : Connect two splines into one.'
%                'Viewing box:'; 'Use MATLAB default icons to rotate or zoom in/out globally. ' ;
%                'Whenever use the keyboard  to interact, remember to cancel rotate/zoom mode!! ';
%                'x y z limits can be changed individually by keys,[Q][W][E][A][S][D](case sensitive).'; 
%                '[Q][A] = +- X direction ' ; '[W][S] = +- Y direction '; '[E][D] = +- Z direction \bf{(Not for 2D diagram)}\rm'  ;'Lower cases: Shifting(- or +) the limit in the corresponding direction.' ; 'Upper cases: Expand(+)/shrink(-) the limit in the corresponding direction. '    
%                '';

              
               
               } ;
                      [icondata,iconcmap] = imread([pwd filesep 'images' filesep 'SplineGUI.jpg']);     
                      f = msgbox(STR ,'Instructions', 'custom',icondata,iconcmap,CreateStruct);    
%                        f = msgbox(STR ,'Instructions', 'help' ,CreateStruct);    
%                       f.Position
% f.Position(3) =0.5 ; f.Children(3).Units='normalized'; f.Children(3)
% f.Children(3).Position(3) = 1- f.Children(3).Position(1) -0.05;
% axis auto;
% movegui(f,'center');
       end

% \color[rgb]{0,0.5,0}
end

