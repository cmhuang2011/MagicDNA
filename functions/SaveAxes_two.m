%SaveAxes

ax_old = gca;

f_new = figure; 




ax_new = copyobj(ax_old,f_new)  ;
set(ax_new,'Position','default');
ax_new.Tag='';
set(gca,'FontSize',20)

fprintf('Use ex: print(''-r100'',f_new,''Image3_r100'',''-djpeg'') \n' ) ;
fprintf('Use ex: set(gcf,''KeyPressFcn'',@(src,evn)ChangeAxesLimit(src,evn) ) \n' ) ;

if isfield(ax_old.UserData ,'hlink')
ax_old_hidden = ax_old.UserData.hlink.Targets(2) ;
ax_new2 =  axes('Position',ax_new.Position,'Visible','off','hittest','off'); % Invisible axes
% ax_new2 =  axes('Parent', f_new,'Position',ax_new.Position,'Visible','off','hittest','off'); % Invisible axes
% ax_new2 =axes('Position',ax_new.Position,'Visible','off','hittest','off') ;
h_text = findobj(ax_old_hidden,'Type','Text')  ;
h_text_new = copyobj(h_text,ax_new2)  ;

hlink =linkprop([ax_new ax_new2],{ 'XLim' 'YLim' 'ZLim' ,'PlotBoxAspectRatio','View' }); % The axes should stay aligned
ax_new.UserData.hlink= hlink ;
% axes(ax_new) ;
end

set(f_new,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) ) 