function  LegendBoxing_Scaf( src,evn,axMain )
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
           STR={'\fontsize{11}Click \bf{GetScaffoldRouting}\rm button to obtain a scaffold routing.'  ;
               'The algorithm gives a new routing according to how users specified geometry and assembly.';
               'Routing results are different each time.';
               'Before doing this step, the mechanism must be assembled with desired lengths of ssDNA scaffold in Assembly and ssScaf tabs.';
               '';
               '\fontsize{11}Click \bf{MultiScaf}\rm button to split \bfone\rm long scaffold into N cycles by applying N-1 Xovers.'  ;
               'This is a stochastic method to find the Xovers satisfying the length constraints.'
               'This means the process may not be successful and requires multiple trials. '
               'Default setting is to use the specified N and percentage for the lower and upper bounds, which can be changed in the file \bffunctions/SplitScafToMultiScaf.m.\rm'
               
               '';
               'Viewing box:'; 'Use MATLAB default icons to rotate or zoom in/out globally.' ;
               'Whenever use the keyboard to interact, \bfremember to cancel rotate/zoom mode!!\rm ';
               'x y z limits can be changed individually  by keys,[Q][W][E][A][S][D](case sensitive).'; 
               '[Q][A] = +- X direction ' ; '[W][S] = +- Y direction '; '[E][D] = +- Z direction '  ;'Lower cases: Shifting(- or +) the limit in the corresponding direction.' ; 'Upper cases: Expand(+)/shrink(-) the limit in the corresponding direction. '    

               '';
               '[x]: axis equal ';
               '[X]: axis auto ';
               '[p]: print current view under directory. ';
               '[h]: watch tutorial movie.' ;
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

