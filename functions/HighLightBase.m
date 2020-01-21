function HighLightBase(src,evn,sH3D)
evn.IntersectionPoint ;
XY=[src.XData; src.YData]' ;
ax= findobj(gcf,'Tag','json3D') ; 
ax2=findobj(gcf,'Tag','json2D') ;
Json_Text=findobj(gcf,'Tag','JSON_Text3') ;
GetHyperB= src.UserData.GetHyperB ;

[~,ind] = ismember( evn.IntersectionPoint(1:2) ,XY ,'rows' )  ;
if ~isfield( src.UserData , 'Highlight')
    axes(ax2) ;
    src.UserData.Highlight= scatter(evn.IntersectionPoint(1) ,evn.IntersectionPoint(2),120,'or' ,'filled') ;
%      m = uimenu(src.UserData.Highlight,'Text','Show');
     
    c = uicontextmenu;
    % Assign the uicontextmenu to the figure
    src.UIContextMenu = c;
    % Create child menu of the uicontextmenu
    topmenu = uimenu('Parent',c,'Label','Step Size');    
    % Create submenu items
    m1 = uimenu('Parent',topmenu,'Label','1','Callback',@(source,callbackdata)changeSS(source,callbackdata,src.UserData.Highlight));
    m2 = uimenu('Parent',topmenu,'Label','10','Callback',@(source,callbackdata)changeSS(source,callbackdata,src.UserData.Highlight));
    m3 = uimenu('Parent',topmenu,'Label','50','Callback',@(source,callbackdata)changeSS(source,callbackdata,src.UserData.Highlight));
     
     axes(ax) ;
     src.UserData.Highlight3D= scatter3(sH3D.XData(ind) ,sH3D.YData(ind),sH3D.ZData(ind),120,'or' ,'filled') ;
     src.UserData.ind=ind ;
     src.UserData.Highlight.UserData.Step= 1;
     
else
     src.UserData.Highlight.XData(1) =evn.IntersectionPoint(1) ;
     src.UserData.Highlight.YData(1) =evn.IntersectionPoint(2) ;
     src.UserData.Highlight3D.XData(1) = sH3D.XData(ind) ;
     src.UserData.Highlight3D.YData(1) = sH3D.YData(ind) ;
     src.UserData.Highlight3D.ZData(1) = sH3D.ZData(ind) ;
     src.UserData.ind=ind ;
     
    
end
%  ind
C5ind =  src.UserData.SFC(ind ,:) ;
if ind+1<=size(src.UserData.SFC,1)
C5indp1 =  src.UserData.SFC(ind+1 ,:) ;
else
C5indp1 =  src.UserData.SFC(ind ,:) ;    
end

BCB_C5ind=[ GetHyperB.RelateTable( GetHyperB.RelateTable(:,5)==C5ind(1) ,1:2) ,C5ind(2) ] ;
BCB_C5indp1=[ GetHyperB.RelateTable( GetHyperB.RelateTable(:,5)==C5indp1(1) ,1:2) ,C5indp1(2) ] ;

src.UserData.ScafBCB = [ BCB_C5ind,BCB_C5indp1];
[ BCB_C5ind,BCB_C5indp1] ;

StrShow= sprintf('Cadnano label= %i[%i] ' , src.UserData.SFC_C4notation(ind ,:) ) ;
% StrShow = [StrShow newline 'Scaf Position = ' num2str(ind)] ;


[aa,~] = cellfun(@size,GetHyperB.ScafAllBase) ; 

Caa = cumsum(aa) ;
WhichScaf = ind<=cumsum(aa) ; 
WhichScaf = find(WhichScaf) ;
if isempty(WhichScaf)
    WhichScaf=1 ;
else
    WhichScaf=WhichScaf(1)  ;
end

aa=[0;Caa];

StrShow = [StrShow newline 'Scaf(' num2str(WhichScaf)  ')Position = ' num2str(ind-aa(WhichScaf))] ;

% StrShow{2} =' test '
% sdfsf=3
% 
Json_Text.String= StrShow ;
axes(ax2);
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
%         sdfsdf=3
    end
