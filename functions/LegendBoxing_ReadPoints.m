function  LegendBoxing_ReadPoints( src,evn,axMain )
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
           STR={'\fontsize{11}Left click to make a red node.'  ;
               'Right click to make a green node.';
               'Add a line between the red and green nodes by clicking the button.';
               'When done, click ''Send to Sketch''.' ;
               '';
               '[h]: watch tutorial movie.';
              
               
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

