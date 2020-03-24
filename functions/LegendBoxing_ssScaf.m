function  LegendBoxing_ssScaf( src,evn,axMain )
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
       
%        [icondata,iconcmap] = imread('trees.tif'); 
       
       if evn.IntersectionPoint(2)>0.5
           CreateStruct.Interpreter = 'tex';
           CreateStruct.WindowStyle = 'modal';
           STR={'\fontsize{11}(1) For cylinder ends without connections, assign the extended ssDNA lengths from both ends to form scaffold loops in the top-right popup menu. The total lengths of scaffold loops will be two times of the number of bases at the end of each bundle.';
               'Recommended this value as least 4 nts.';
               '(2) Use the keyboard to zoom in/out and the mouse to select connections to specify lengths of scaffold ssDNA.';
             
               '';
               'Click on the texts in the right to change font size if necessary.';
               '';
               'Viewing box:'; 
               'Use MATLAB default icons to rotate or zoom in/out globally. ' ;
               'Whenever use the keyboard to interact, remember to cancel rotate/zoom mode!! ';'';
               'x y z limits can be changed individually by keys,[Q][W][E][A][S][D](case sensitive) or by clicking on the bottom half of the legend(left to right, -+x, -+y, -+z).'; 
               '[Q][A] = +- X direction '; 
               '[W][S] = +- Y direction '; 
               '[E][D] = +- Z direction '  ;
               'Lower cases or left clicks on the bottom half: Shifting(- or +) the limit in the corresponding direction.' ; 'Upper cases or middle clicks on the bottom half: Expand(+)/shrink(-) the limit in the corresponding direction. '    
               '';


               '[x]: axis equal ';
               '[X]: axis auto ';
               '[p]: print current view under directory. ';
               '[h]: watch tutorial movie.'
               } ;
                                 
%            f = msgbox(STR ,'Instructions', 'custom',icondata,iconcmap ,CreateStruct);    
                      f = msgbox(STR ,'Instructions', 'help' ,CreateStruct);    

%            f.
%            f.Children(2).Position(3:4) = [48 48] ;
%                       f = msgbox(STR ,'Instructions','help' ,CreateStruct ,'custom',icondata,iconcmap);                      

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

end

