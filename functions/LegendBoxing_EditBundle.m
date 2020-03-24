function  LegendBoxing_EditBundle( src,evn,axMain )
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
    CreateStruct.WindowStyle = 'model';
    STR={'\fontsize{12}Left click on the nodes on 2D or 3D panel to select. Selected nodes will turn red.'  ;
        'Right click to cancel or use the buttons.';
        'Use keyboard left or right arrow to move the selected nodes(red).';
        '2D and 3D nodes are synchronized. Change step size by "Ctrl" or "Shift".' ;
        '';
        '2D cylinder view :';
        '(1) For specifying internal crossovers, use the popups next to the buttons to select one cylinder or a pair of cylinder. The available crossovers from the cylinder or between these two cylinder will be shown. Click the desired one and press [x] to turn it into \color{red}\bf bold red \color{black}\rm as saved.';
        '(2) For ignore crossovers along a portion of bundle, drag the cursor and click [z] to extend the range to two sides with inputs. These two above operations can be canceled by right clicks.' ;
        %                '  '
        '';
        'Cross-section view :';
        '(1) Measure the distance between two cylinders on cross-section by clicking.' ;
        '(2) Zoom in 3D model to see the scaffold orientation.' ;
        '(3) If users want to change the cylinder pairing(scaffold loops), select two pairing lines and press ''Enter''. ';
        '(4) For ignore all crossovers between two cylinders, use left/right click on the circle between cylinders. \color{red}Red \color{black}means no crossover.  ';
        '';
        };
    
    STR2={'\fontsize{12}Viewing box:';
        'Use MATLAB default icons to rotate or zoom in/out globally. ' ;
        'Whenever use the keyboard to interact, remember to cancel rotate/zoom mode!! ';'';
        'x y z limits can be changed individually by keys,[Q][W][E][A][S][D](case sensitive) or by clicking on the bottom half of the legend(left to right, -+x, -+y, -+z).';
        '[Q][A] = +- X direction ';
        '[W][S] = +- Y direction ';
        '[E][D] = +- Z direction \bf{(Not for 2D and cross-section view)}\rm'  ;
        'Lower cases or left clicks on the bottom half: Shifting(- or +) the limit in the corresponding direction.' ; 'Upper cases or middle clicks on the bottom half: Expand(+)/shrink(-) the limit in the corresponding direction. '
        '';
        
        'Keys:';
        '[x][X]: Add forced internal crossovers. (Multi-scaffold, top-right)'
        '[z][Z]: deactivate the crossovers in certain Z range (Multi-scaffold, top-right)'
        '[o]: axis equal \bf{(recommended for 3D panel and cross-section)}\rm';
        '[O]: axis auto ';
        '[n][N]: axis normal \bf{(recommended for 2D panel)}\rm';
        '[p]: print current view under directory.';   
        '[h]: watch tutorial movie.';
        '\color{red}\fontsize{16}\bf{Need to update the entire mechanism(Use ''Add/Remove'' button in Assembly tab without adding or removing bundles) after editing all bundles!!}\rm ';
        
        } ;
    
    %            f = msgbox(STR ,'Instructions', 'custom',icondata,iconcmap ,CreateStruct);
    f = msgbox(STR ,'Instructions(1/2)', 'help' ,CreateStruct);      % f.Children(2).Position(3:4) = [48 48] ;
    f2 = msgbox(STR2 ,'Instructions(2/2)', 'help' ,CreateStruct);      % f.Children(2).Position(3:4) = [48 48] ;
    
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

