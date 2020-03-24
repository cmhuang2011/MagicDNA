function  LegendBoxing_SketchStep2( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
%        evn
       Rulers= linspace(0,1,7) ;
       Locats= sum(evn.IntersectionPoint(1)>=Rulers);
       
%        [icondata,iconcmap] = imread('trees.tif'); 
       
       if evn.IntersectionPoint(2)>0
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';
           STR={'\fontsize{11}This step involves converting line model to cylinder model(bundles) by assigning lengths, cross-sections, gradiences, and initial orientations.The 3D lines are related to the table shown in right.'  ;
               '\color{red}{Red(x)}, \color{green}{green(y)}, \color{black}and \color{blue}blue(z) lines\color{black} refer the local coordinate of each line.';
               '\color{blue}{Blue lines}\color{black} are pointing from Z1 sides to Z2 sides.'              
               'Each row in the table provides the required parameters of converting lines to bundles.';
               'Detailed descriptions are available in the table tooltip.' ;
               '';
               'UI functions:' ; '';
               'Local coordinate can be flipped in x and z directions by left clicks.';
               'Right clicks will rotate x and y directions along z by 90 degrees.';
               'The top slider can scale the lengths of each bundle.' ;
               'The bottom slider can explode all bundles for better visualization in Assembly tab.(suggested value between 1.25~1.35)'
               '';
               'Viewing box:'; 'Use MATLAB default icons to rotate or zoom in/out globally.' ;
               'Whenever use the keybroad to interact, remember to cancel rotate/zoom mode!! ';
               'x y z limits can be changed inividually by keys,[Q][W][E][A][S][D](case sensitive).'; 
               '[Q][A] = +- X direction ' ; '[W][S] = +- Y direction '; '[E][D] = +- Z direction '  ;'Lower cases: Shifting(- or +) the limit in the corresponding direction.' ; 'Upper cases: Expand(+)/shrink(-) the limit in the corresponding direction. '    

               '';
               '[x]: axis equal ';
               '[X]: axis auto ';
               '[p]: print current view under directory. ';
               '[h]: watch tutorial movie.'
               } ;
                                 
           f = msgbox(STR ,'Instructions', 'help' ,CreateStruct);    

       end
%        
%        if evn.Button==1
%            switch Locats
%                case 1
%                    ax.XLim= ax.XLim - interval ;      
%                case 2
%                      ax.XLim= ax.XLim + interval ;     
%                case 3
%                      ax.YLim= ax.YLim - interval ;
%                case 4
%                      ax.YLim= ax.YLim + interval ;
%                case 5
%                      ax.ZLim= ax.ZLim - interval ;
%                case 6
%                      ax.ZLim= ax.ZLim + interval ;
%            end
%        elseif evn.Button==2
%            switch Locats
%                case 1
%                    ax.XLim= (ax.XLim -mean(ax.XLim))/R + mean(ax.XLim)   ;   
%                case 2
%                    ax.XLim= R*(ax.XLim -mean(ax.XLim))+ mean(ax.XLim)   ;   
%                case 3
%                    ax.YLim= (ax.YLim -mean(ax.YLim))/R + mean(ax.YLim)   ;    
%                case 4
%                    ax.YLim= R*(ax.YLim -mean(ax.YLim))+ mean(ax.YLim)   ;     
%                case 5
%                      ax.ZLim= (ax.ZLim -mean(ax.ZLim))/R + mean(ax.ZLim)   ;     
%                case 6
%                      ax.ZLim= R*(ax.ZLim -mean(ax.ZLim))+ mean(ax.ZLim)   ;     
%            end 
%        else
%            ax.XLim=ax.XLim+0.01 ;
%        end

end

