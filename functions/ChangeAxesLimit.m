function ChangeAxesLimit(src,evn)

% %     set(gcf,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) )

ax=gca;
interval =5; R=1.1;
switch  evn.Character
    case 'q'
        ax.XLim= ax.XLim + interval ;
    case 'a'
        ax.XLim= ax.XLim - interval ;
    case 'Q'
        ax.XLim= R*(ax.XLim -mean(ax.XLim))+ mean(ax.XLim)   ;
    case 'A'
        ax.XLim= (ax.XLim -mean(ax.XLim))/R + mean(ax.XLim)   ;
        %-----------
    case 'w'
        ax.YLim= ax.YLim + interval ;
    case 's'
        ax.YLim= ax.YLim - interval ;
    case 'W'
        ax.YLim= R*(ax.YLim -mean(ax.YLim))+ mean(ax.YLim)   ;
    case 'S'
        ax.YLim= (ax.YLim -mean(ax.YLim))/R + mean(ax.YLim)   ;
        %------------
    case 'e'
        ax.ZLim= ax.ZLim + interval ;
    case 'd'
        ax.ZLim= ax.ZLim - interval ;
    case 'E'
        ax.ZLim= R*(ax.ZLim -mean(ax.ZLim))+ mean(ax.ZLim)   ;
    case 'D'
        ax.ZLim= (ax.ZLim -mean(ax.ZLim))/R + mean(ax.ZLim)   ;
        
        %---------------
    case 'x'
        axis equal;
        fprintf('set axis equal \n')
    case 'X'
        axis auto
        fprintf('set axis auto \n')
    case 'P'
        fprintf('Capturing..... \n') ;
        print('-r100',gcf,'SnapShot_r100','-dpng')  ;
        fprintf('Captured current figure. \n') ;
    case 'p'
        fprintf('Capturing..... \n') ;        
        print('-r100',gcf,'SnapShot_r100','-dpng')  ;
        fprintf('Captured current figure. \n') ;
    case 'n'
        axis normal
        fprintf('set axis normal \n')
    case 'N'
        axis normal
        fprintf('set axis normal \n')
end

end