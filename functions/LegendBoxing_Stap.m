function  LegendBoxing_Stap( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
       if evn.IntersectionPoint(2)>0
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';
           STR={'\fontsize{11}Click \bf{GetStaple}\rm button to obtain staple routings.';
               'The default and recommend type is "zigzag cut", which applies all staple cross-overs to the staple strands, and breaks the long staples into short pieces.';
               'The halfXover option is to connect the staples ends if the 5'' end on a staple is close enough to the 3'' end on another staple \bfwithin the same bundle\rm for reducing fraying.';
               'If ssDNA on scaffold is set to 0-nt for end-to-end connections, the algorithm automatically connects them on staples for wire-frame and surface-modeling purposes.' ;
               'This is the only situation where staples may across two bundles, except using caDNAno to modify manually.';
               '';
               'It is optional to find staple routing with overhangs if users have completed the ''Overhang design'' tab.';
               '';
               '''Use cadnano'' button is for manually change the routing in caDNAno with the json file \bfwhich was exported from MagicDNA.\rm';
               'Please load the json file while the corresponding MagicDNA design is loaded.';
               '\bfLimitation: \rm';
               '1. Scaffold and staple strands can be cut, added Xovers, or modified in caDNAno but must have 3'' and 5'' ends (breaks).';
               '2. Do not change the helices indexes(Rnum) in caDNAno.' ;
               '3. If with suitable operations after modification, some preliminary analysis will show up in the command line, which means both scaffold and staples have been updated. '
                '';
               'Viewing box:'; 'Use MATLAB default icons to rotate or zoom in/out globally.' ;
               'Whenever use the keyboard to interact, \bfremember to cancel rotate/zoom mode!!\rm ';
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

% \color[rgb]{0,0.5,0}


end

