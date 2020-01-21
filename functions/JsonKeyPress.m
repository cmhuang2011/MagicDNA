function JsonKeyPress(src,evn)

% %     set(gcf,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) )

ax=gca ;
interval =5; R=1.1;
evn;
switch evn.Key
    case 'rightarrow'  % for cadnano tab, mapping 2D base with 3D position
        sH3D=findobj(gcf,'Tag','sH3D')  ;
        sH=findobj(gcf,'Tag','sH')  ;
        if isempty(sH);return;end
        stp=sH.UserData.Highlight.UserData.Step ;
        NewInd = sH.UserData.ind +stp;
        if NewInd<=length(sH.XData)
            NewXY =[sH.XData(NewInd) , sH.YData(NewInd)  ] ;
            evnSH.IntersectionPoint=NewXY;
            HighLightBase(sH,evnSH,sH3D) ;
        end
        return
    case 'leftarrow'
        sH3D=findobj(gcf,'Tag','sH3D')  ;
        sH=findobj(gcf,'Tag','sH')  ;
        if isempty(sH);return;end
        stp=sH.UserData.Highlight.UserData.Step ;
        NewInd = sH.UserData.ind -stp;
        if NewInd>0
            NewXY =[sH.XData(NewInd) , sH.YData(NewInd)  ] ;
            evnSH.IntersectionPoint=NewXY;
            HighLightBase(sH,evnSH,sH3D) ;
        end
        return
end

switch  evn.Character
    case 'h'
        implay('CadnanoInterface.mp4');
    case 't'
        t_json=findobj(gcf,'Tag','t_json') ;
        t_json.UserData;
        plotH=t_json.UserData.plotH;
        pStapleH=t_json.UserData.pStapleH;
        
        c = uisetcolor ;
        MultiSelectedStaple = unique( t_json.UserData.SelectedCell(:,1) ) ;
        for k=1:length(MultiSelectedStaple)
            Sti=MultiSelectedStaple(k) ;
            pStapleH{Sti}.Color =c ;
            plotH{Sti}.Color =c ;
            t_json.BackgroundColor(Sti,:)=c ;
            Scale256 = c* 255 ;
            HexCode=strcat( dec2hex(round(Scale256(1)),2) , dec2hex(round(Scale256(2)),2) , dec2hex(round(Scale256(3)),2 )   ) ;
            
            t_json.Data{Sti,5} = HexCode ;
            
        end
        
        %         sdsdf=3
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
        if isfield(ax.UserData, 'NoZ');return ;end
        ax.ZLim= ax.ZLim + interval ;
    case 'd'
        if isfield(ax.UserData, 'NoZ');return ;end
        ax.ZLim= ax.ZLim - interval ;
    case 'E'
        if isfield(ax.UserData, 'NoZ');return ;end
        ax.ZLim= R*(ax.ZLim -mean(ax.ZLim))+ mean(ax.ZLim)   ;
    case 'D'
        if isfield(ax.UserData, 'NoZ');return ;end
        ax.ZLim= (ax.ZLim -mean(ax.ZLim))/R + mean(ax.ZLim)   ;
        
        %---------------
    case 'o'
        axis equal;
        fprintf('set axis equal \n')
    case 'O'
        axis auto
        fprintf('set axis auto \n')
    case 'P'
        fprintf('Capturing..... \n') ;
        print('-r100',gcf,'SnapShotFromJson_r100','-dpng')  ;
        fprintf('Captured current figure. \n') ;
    case 'p'
        fprintf('Capturing..... \n') ;
        print('-r100',gcf,'SnapShotFromJson_r100','-dpng')  ;
        fprintf('Captured current figure. \n') ;
    case 'n'
        axis normal
        fprintf('set axis normal \n')
    case 'N'
        axis normal
        fprintf('set axis normal \n')
    case 'z'
        ggtab =gctab  ; %
        ss_json=findall(gcf,'Tag','ss_json') ;
        if isequal(ggtab,ss_json)             % curent tab on cadnano
            sH=findobj(ss_json,'Tag','sH') ;
            if ~isempty(sH)   % have initial the routing
                if isfield(sH.UserData,'ScafBCB') % have create red dot in cadnano tab
                    if sH.UserData.SFC(sH.UserData.ind,1)==sH.UserData.SFC(sH.UserData.ind-1,1)  % select (kth,k+1) base, avoid routing on single base
                        XY2d= [ sH.XData( sH.UserData.ind) , sH.YData( sH.UserData.ind) ] ;
                        
                        sH.UserData.Highlight.XData(2)=XY2d(1) ;
                        sH.UserData.Highlight.YData(2)=XY2d(2) ;
                        sH.UserData.Highlight.CData(1:2,1:3) =[1,0,0 ;  0 , 1 , 0.3];
                        sH.UserData.ScafBCBz=  sH.UserData.ScafBCB ;
                        %                                 sH.UserData
                        %                                 fprintf(' key z \n')
                    else
                        f = msgbox('Not valid base.');
                    end
                end
            end
        end
    case 'x'
        ggtab =gctab  ; %
        ss_json=findall(gcf,'Tag','ss_json') ;
        if isequal(ggtab,ss_json)             % curent tab on cadnano
            sH=findobj(ss_json,'Tag','sH') ;
            if ~isempty(sH)   % have initial the routing
                if isfield(sH.UserData,'ScafBCB') % have create red dot in cadnano tab
                    if sH.UserData.SFC(sH.UserData.ind,1)==sH.UserData.SFC(sH.UserData.ind-1,1)  % select (kth,k+1) base, avoid routing on single base
                        XY2d= [ sH.XData( sH.UserData.ind) , sH.YData( sH.UserData.ind) ] ;
                        if length(sH.UserData.Highlight.XData) ~= 1 % make sure use z key before x key
                            sH.UserData.Highlight.XData(3)=XY2d(1) ;
                            sH.UserData.Highlight.YData(3)=XY2d(2) ;
                            sH.UserData.Highlight.CData =[1,0,0;  0 , 1 , 0.3 ;  0 , 1 , 0.8];
                            sH.UserData.ScafBCBx=  sH.UserData.ScafBCB ;
                        end
                    else                                    f = msgbox('Not valid base.');
                        
                    end
                    %                                 fprintf(' key x \n')
                end
            end
        end
    case 'v'
        ggtab =gctab  ; %
        ss_json=findall(gcf,'Tag','ss_json') ;
        
        if isequal(ggtab,ss_json)             % curent tab on cadnano
            sH=findobj(ss_json,'Tag','sH') ;
            if ~isempty(sH)   % have initial the routing
                if isfield(sH.UserData,'ScafBCB') % have create red dot in cadnano tab
                    XY2d= [ sH.XData( sH.UserData.ind) , sH.YData( sH.UserData.ind) ] ;
                    if length(sH.UserData.Highlight.XData) == 3 % make sure have used 'z' and 'x'
                        if ~isfield(sH.UserData,'ExtraForcedScafXover')
                            axJson2D=findobj(gcf,'Tag','json2D') ; axes(axJson2D) ;
                            
                            MM =[sH.UserData.ScafBCBz ;  sH.UserData.ScafBCBx  ] ;
                            MM= [ MM(2,4:6), MM(1,1:3) ; MM(1,4:6), MM(2,1:3) ];
                            %                                                     MM= [ MM(1,1:3), MM(2,4:6) ; MM(2,1:3), MM(1,4:6) ];
                            sH.UserData.ExtraForcedScafXover = MM ;
                        else
                            axJson2D=findobj(gcf,'Tag','json2D') ; axes(axJson2D) ;
                            MM =[sH.UserData.ScafBCBz ;  sH.UserData.ScafBCBx  ] ;
                            MM= [ MM(2,4:6), MM(1,1:3) ; MM(1,4:6), MM(2,1:3) ];
                            %                                                     MM= [ MM(1,1:3), MM(2,4:6) ; MM(2,1:3), MM(1,4:6) ];
                            sH.UserData.ExtraForcedScafXover = [sH.UserData.ExtraForcedScafXover ; MM] ;
                        end
                        plot(sH.UserData.Highlight.XData(2:3),sH.UserData.Highlight.YData(2:3),'-o','Color','k' ,'MarkerFaceColor',ones(1,3)) ;
                        sH.UserData.Highlight.XData=sH.UserData.Highlight.XData(1) ;
                        sH.UserData.Highlight.YData=sH.UserData.Highlight.YData(1) ;
                        sH.UserData.Highlight.CData=[1,0 ,0];
                        sH.UserData = rmfield(sH.UserData,'ScafBCBz') ;
                        sH.UserData = rmfield(sH.UserData,'ScafBCBx') ;
                        
                    end
                end
            end
        end
        
end

end