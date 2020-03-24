function  LegendBoxing_cadnano( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
       if evn.IntersectionPoint(2)>0
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';%CreateStruct.Units = 'normalized'; %CreateStruct.Position = 'normalized';
           STR={'\fontsize{11}caDNAno tab: Preview the routings(both scaffold and staple) and specify scaffold sequences by ''Preview'' button.';
               'The 3D design model is juxtaposed with 2D design diagram, which mimics caDNAno design diagram.';
               'The purpose of this tab is for visualizing the design, including staple locations in the 3D structure, exporting staple sequences, and helping users to modify the design in caDNAno software.';
               '';
               'UI functions: ';
               '(1) Both scaffold and staples are able to change transparency(use sliders) and representation(switch chickenwire or helical representation by clicking the texts).';
               '(2) Click on the staple list to see the location of that staple in both 2D and 3D routings.  ' ;
               '(3) Click on the scaffold strand in 2D or 3D panel. Two red dots on 2D and 3D routing indicate the mapping location. Keyboard leftarrow and rightarrow can move along the scaffold strand. Right click on the 2D scaffold to change the step size.' ;
               '(4) Exporting Json files will generate two .json files with identical routings due to hybrid lattice.';
               '(5) Users can change staple colors by selecting staples(multiple) in the list, then clicking on somewhere in the tab not the table, and press \bf[t]\rm to change colors.';
               '(6) If users have used overhang tool, the ''OverHangOption'' button will be active. The ''Optimized Seq'' option generates sequences that have minimal complementarity to the scaffold, especially single stranded portions.';
               '(7) Staple colors can be assigned according to the 5'' ends in the bundles or the staple graph(considering nicks to form a long chain).'
               '(8) For multi-scaffold, the numbers of complementary bases between I scaffolds and J staples are shown in the GUI. Move the mouse to see the values and highlight the staple in the main window. Notice a analysis report in command line. '
               '';
               'Viewing box:'; 'Use MATLAB default icons to rotate or zoom in/out globally. ' ;
               'Whenever use the keyboard  to interact, remember to cancel rotate/zoom mode!! ';
               'x y z limits can be changed individually by keys,[Q][W][E][A][S][D](case sensitive).'; 
               '[Q][A] = +- X direction ' ; '[W][S] = +- Y direction '; '[E][D] = +- Z direction \bf{(Not for 2D diagram)}\rm'  ;'Lower cases: Shifting(- or +) the limit in the corresponding direction.' ; 'Upper cases: Expand(+)/shrink(-) the limit in the corresponding direction. '    
               '';
               '[o]: axis equal \bf{(recommended for 3D panel)}\rm';
               '[O]: axis auto ';
               '[n][N]: axis normal \bf{(recommended for 2D panel)}\rm';
               '[p]: print current view under directory.';      
               '[h]: watch tutorial movie.';
              
               
               } ;
                                 
                      f = msgbox(STR ,'Instructions', 'help' ,CreateStruct);    
%                       f.Position
% f.Position(3) =0.5 ; f.Children(3).Units='normalized'; f.Children(3)
% f.Children(3).Position(3) = 1- f.Children(3).Position(1) -0.05;
% axis auto;
% movegui(f,'center');
       end

% \color[rgb]{0,0.5,0}
end

