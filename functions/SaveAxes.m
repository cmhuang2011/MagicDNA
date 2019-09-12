%SaveAxes

ax_old = gca;
% cbar_axes = colorbar;

f_new = figure;
ax_new = copyobj(ax_old,f_new)  ;
% cbar_new= copyobj(cbar_axes,f_new)
% cbar_axes.Parent=f_new
set(ax_new,'Position','default');
ax_new.Tag='';
set(gca,'FontSize',20)
% print('-r100',f_new,'Image3_r100','-djpeg')

fprintf('Use ex: print(''-r100'',f_new,''Image3_r100'',''-djpeg'') \n' ) ;
fprintf('Use ex: set(gcf,''KeyPressFcn'',@(src,evn)ChangeAxesLimit(src,evn) ) \n' ) ;

set(f_new,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) ) 