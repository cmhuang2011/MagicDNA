function printSTEPtable(varargin)

% printSTEPtable

%use for save parameters of exporting mechanism. output is a image with
%parameters later becoming a reference

if ~isempty(varargin)
    str=varargin{1} ;
else
    str='STEPtable' ;    
end

tableSTEP= findobj(0,'Tag','STEPtable') ;
f_new2 = figure;
drawnow;warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFig = get(handle(f_new2), 'JavaFrame'); jFig.setMaximized(true); drawnow;
t = uitable(f_new2,'Units','normalized','Position',[0.1 0.1 0.8 0.7],....
    'ColumnWidth','auto', 'ColumnEditable', true,'Data',tableSTEP.Data,....
    'ColumnName',{'Edge'; 'BP' ;'Nx/HC  ';'Ny';'dZ1/dNx';'dZ1/dNy';'dZ2/dNx';'dZ2/dNy' ;'shift'   },'FontSize',14 );
% print('-r300',f_new2,'STEPtable','-djpeg')

txtjson3 = uicontrol(gcf,'Style','text','Unit','normalized','FontSize',12,'Position', [0.1000    0.800    0.8000    0.1000],....   
        'String',str,'Enable','Inactive','HorizontalAlignment','left','BackgroundColor','w','FontSize',18);   

print('-r300',f_new2,str,'-djpeg')

close(f_new2)

%  set(gca,'FontSize' ,20)

end