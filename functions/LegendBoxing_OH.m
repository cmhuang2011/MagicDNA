function  LegendBoxing_OH( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
       if evn.IntersectionPoint(2)>0
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';
           STR={'\fontsize{11}This tab is for designing overhangs on staples directly in 3D model.';
               'Initialize the tool by clicking the ''OHfunction'' button(Later becomes ''Clear'' button).'
               'Where helixes point outward are possible sites for extending overhangs, shown in \bf{black dots}\rm.';
               'Use mouse left(\color{red}{red}\color{black}) and right(\color{blue}{blue}\color{black}) click to select two nodes.';
               'Hit "Add" botton to connect two nodes as a pair. Repeat until done.';
               'Modify the table to specify the parameters for the staple algorithm.';
               'Go to the "Staple" tab to find the staple routing with overhangs.';
               '';
               'UI functions';
               'Use \bf{Highlight column}\rm to see the data in the table mapping with 3D position in the right.';
               'If the pair of overhangs in \bf{Enable column}\rm is not checked, the algorithm will temporarily ignore it.';
               'Mouse: middle or right click to switch bolded(highlight) and thin lines. '
%                'Check the setting of a pair of overhang by clicking on thin lines.'
               'Click on thin lines to highlight the data in the table.';
               'Delete the pair of overhang by clicking on bolded lines.'
               
               '';
               'Viewing box:'; 'Use MATLAB default icons to rotate or zoom in/out globally.' ;
               'Whenever use the keyboard  to interact, remember to cancel rotate/zoom mode!! ';
               'x y z limits can be changed individually  by keys,[Q][W][E][A][S][D](case sensitive).'; 
               '[Q][A] = +- X direction ' ; '[W][S] = +- Y direction '; '[E][D] = +- Z direction '  ;'Lower cases: Shifting(- or +) the limit in the corresponding direction.' ; 'Upper cases: Expand(+)/shrink(-) the limit in the corresponding direction. '    

                '';
               '[x]: axis equal ';
               '[X]: axis auto ';
               '[p]: print current view under directory. ';
               '[h]: watch tutorial movie.';
               } ;
                                 
               f = msgbox(STR ,'Instructions', 'help' ,CreateStruct);    

       end


end

