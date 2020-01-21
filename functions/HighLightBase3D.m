function HighLightBase3D(src,evn,sH2D)
evn.IntersectionPoint ;
XYZ=[src.XData; src.YData; src.ZData]' ;
ax= findobj(gcf,'Tag','json3D') ; 
ax2=findobj(gcf,'Tag','json2D') ;
Json_Text=findobj(gcf,'Tag','JSON_Text3') ;
GetHyperB= sH2D.UserData.GetHyperB ;


evn.IntersectionPoint(1:3);
d= XYZ- ones(size(XYZ,1),1 )*evn.IntersectionPoint(1:3) ;
dd= sqrt(d(:,1).^2 +d(:,2).^2 +d(:,3).^2 ) ;
ind = find(dd==min(dd)) ;
% [~,ind] = ismember( evn.IntersectionPoint(1:3) ,XYZ ,'rows' )  ;
if ~isfield( sH2D.UserData , 'Highlight')
    axes(ax2) ;
    sH2D.UserData.Highlight= scatter(sH2D.XData(ind)  ,sH2D.YData(ind) ,120,'or' ,'filled') ;
%      m = uimenu(src.UserData.Highlight,'Text','Show');
     
    c = uicontextmenu;
    % Assign the uicontextmenu to the figure
    sH2D.UIContextMenu = c;
    % Create child menu of the uicontextmenu
    topmenu = uimenu('Parent',c,'Label','Step Size');    
    % Create submenu items
    m1 = uimenu('Parent',topmenu,'Label','1','Callback',@(source,callbackdata)changeSS(source,callbackdata,sH2D.UserData.Highlight));
    m2 = uimenu('Parent',topmenu,'Label','10','Callback',@(source,callbackdata)changeSS(source,callbackdata,sH2D.UserData.Highlight));
    m3 = uimenu('Parent',topmenu,'Label','50','Callback',@(source,callbackdata)changeSS(source,callbackdata,sH2D.UserData.Highlight));
     
    axes(ax) ;
     sH2D.UserData.Highlight3D= scatter3(src.XData(ind) ,src.YData(ind),src.ZData(ind),120,'or' ,'filled') ;
     sH2D.UserData.ind=ind ;
     sH2D.UserData.Highlight.UserData.Step= 1;
     
else
     sH2D.UserData.Highlight.XData(1) =sH2D.XData(ind)  ;
     sH2D.UserData.Highlight.YData(1) =sH2D.YData(ind) ;
     sH2D.UserData.Highlight3D.XData(1) = src.XData(ind) ;
     sH2D.UserData.Highlight3D.YData(1) = src.YData(ind) ;
     sH2D.UserData.Highlight3D.ZData(1) = src.ZData(ind) ;
     sH2D.UserData.ind=ind ;
end

C5ind =  sH2D.UserData.SFC(ind ,:) ; C5indp1 =  sH2D.UserData.SFC(ind+1 ,:) ;

BCB_C5ind=[ GetHyperB.RelateTable( GetHyperB.RelateTable(:,5)==C5ind(1) ,1:2) ,C5ind(2) ]  ;
BCB_C5indp1=[ GetHyperB.RelateTable( GetHyperB.RelateTable(:,5)==C5indp1(1) ,1:2) ,C5indp1(2) ] ;

sH2D.UserData.ScafBCB = [ BCB_C5ind,BCB_C5indp1];

[ BCB_C5ind,BCB_C5indp1];
% ind
StrShow= sprintf('Cadnano label= %i[%i] ' , sH2D.UserData.SFC_C4notation(ind ,:) ) ;
StrShow = [StrShow newline 'Scaf Position = ' num2str(ind)] ;
% StrShow{2} =' test '
Json_Text.String= StrShow ;

% sdfsf=3
end

    function changeSS(source,callbackdata,HLPscatter)
        switch source.Label
            case '1'
%                 source.Color = [1.0 0.80 0.80];
              HLPscatter.UserData.Step =1 ;
            case '10'
%                 source.Color = [0.80 1.0 0.80];
             HLPscatter.UserData.Step =10 ;

            case '50'
              HLPscatter.UserData.Step =50 ;

        end

    end
