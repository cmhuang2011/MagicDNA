function  LegendBoxing_oxDNA2( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
       if evn.IntersectionPoint(2)>0
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';
           STR={'\fontsize{11}oxDNA trajectory:';
               'This tab is for visualizing the trajectory from oxDNA simulation. The loaded structure does not need to be the same as the one users are working on.';
               '1. Load the \bftopology, initial configuration, and the trajectory\rm(can be one frame) by clicking the Load button. The order in which files are selected is important.';
               '2. The first time to load the trajectory may take a longer time. Once it has been loaded, MagicDNA saves the trajectory into MATLAB format \bfunder the same folder\rm for faster loading in future. If users use this tab to visualize unfinished simulation earlier, remember to delete Alltraj.mat file when loading the finished trajectory.  '
               '3. Individual configuration in the trajectory can be visualized by selecting the popup menu, multiple selection by the listbox. User can export one configuration(assigned by the popup) to a .bild file for Chimera for better quality.';
               '4. Root-mean-square deviation(RMSD) and root-mean-square fluctuation(RMSF) analyses are included. The average configuration is computed by principle component analysis(PCA, required Statistics and Machine Learning Toolbox) and color-coded with RMSF level shown in the right. Similar to individual frame, this representation can be exported to a .bild file for Chimera as well. ';
               '';
               '[h]: watch tutorial movie.';
               '[g]: Call angle calculation GUI.'
               } ;
                                 
                      f = msgbox(STR ,'Instructions', 'help' ,CreateStruct);    

       end

% \color[rgb]{0,0.5,0}


end

