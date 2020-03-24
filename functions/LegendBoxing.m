function  LegendBoxing( src,evn,axMain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% axMain.UserData
       axes(axMain);
       ax=gca; 
       interval =5; R=1.1;
%        evn
       Rulers= linspace(0,1,7) ;
       Locats= sum(evn.IntersectionPoint(1)>=Rulers);
       Locats(Locats==7)=6 ;
       if evn.IntersectionPoint(2)>0.5
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';
           STR={'\fontsize{11}Use the keyboard and the mouse to interact.'  ;'';
               '\fontsize{10}Bundle manipulation :';'Either use the listbox to select bundles or directly click the 3D bundles.'; 'Multiple selection is possible with ''ctrl'' or middle click. ' 
               'Use the manipulation panel to specify translation or rotation,' ;'step size, and global/local coordinate.'; 
               '[Q][A] = +- X direction ' ; '[W][S] = +- Y direction '; '[E][D] = +- Z direction '  ;'' ;
               'Connectivity: ';'Yellow togglebuttons are for showing the connectivity between bundles.';'To change the connectivity, either:(1)click on lines in graph. (2)edit the table. (3)Use the query panel. '            ;'';
               'Viewing: ' ; 'Use MATLAB default icons to rotate or zoom in/out. ' ; 'x y z limits can change by clicking(left or middle) on the bottom half of the legend.'; '';
               'Helpful tools:' ;'[U]: \bfU\rmndo the previous bundle manipulation'; '[L]: estimate scaffold \bfL\rmength, show in command line'; '[K]: estimate a portion of scaffold length, show in command line '
               '[C]: print the \bfC\rmylinder model as .bild to Chimera'; '[I]: generate an assembly with re-ordered bundle index.'; '[h]: watch tutorial movie.' ;
%                '[M]: save current bundle transformation \bfM\rmatrices and connectivity as Temp.mat under directory. When converting the line models from STEP to Assembly, if the number of edges matches, the program will load this data to expedite local modification.'
               } ;
                                 
           f = msgbox(STR ,'Instructions','help' ,CreateStruct);                      
           return
       end
       
       if evn.Button==1
           switch Locats
               case 1
                   ax.XLim= ax.XLim - interval ;      
               case 2
                     ax.XLim= ax.XLim + interval ;     
               case 3
                     ax.YLim= ax.YLim - interval ;
               case 4
                     ax.YLim= ax.YLim + interval ;
               case 5
                     ax.ZLim= ax.ZLim - interval ;
               case 6
                     ax.ZLim= ax.ZLim + interval ;
           end
       elseif evn.Button==2
           switch Locats
               case 1
                   ax.XLim= (ax.XLim -mean(ax.XLim))/R + mean(ax.XLim)   ;   
               case 2
                   ax.XLim= R*(ax.XLim -mean(ax.XLim))+ mean(ax.XLim)   ;   
               case 3
                   ax.YLim= (ax.YLim -mean(ax.YLim))/R + mean(ax.YLim)   ;    
               case 4
                   ax.YLim= R*(ax.YLim -mean(ax.YLim))+ mean(ax.YLim)   ;     
               case 5
                     ax.ZLim= (ax.ZLim -mean(ax.ZLim))/R + mean(ax.ZLim)   ;     
               case 6
                     ax.ZLim= R*(ax.ZLim -mean(ax.ZLim))+ mean(ax.ZLim)   ;     
           end 
       else
%            ax.XLim=ax.XLim+0.01 ;
       end

uistack(findobj(gcf,'Tag', 'HidenAssemblyMain'),'top') ;
end

