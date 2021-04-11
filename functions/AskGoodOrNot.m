function  Result = AskGoodOrNot(Total_split_ForEachSlice)

% Result.Decision='Snap';
Result.N=Total_split_ForEachSlice;



fh = figure('units','pixels',...
    'position',[300 300 300 200],...
    'menubar','none',...
    'name',' Adjust local',...
    'numbertitle','off',...
    'resize','off');
Result.fh=fh;
        Result.Info='temp' ;

movegui(fh,'center')


edit1 = uicontrol('style','edit','Parent',fh, 'units','pixels','String',Total_split_ForEachSlice,...
    'position',[100 75 100 30]);
txtH2 = uicontrol('Style','text','Parent',fh,'FontSize',12,'Position', [0 110 150 80],....
    'String','Enter min # of bundle. weighted by curvature.');
% { }
btn2 = uicontrol('style','push','Parent',fh, 'unit','pix',...
    'position',[150 20 100 30], 'string','OK',...
    'Callback',@(src,evn)export(src,evn,edit1,fh)    );

btn1 = uicontrol('style','push','Parent',fh, 'unit','pix',...
    'position',[50 20 100 30], 'string','Assign',...
    'Callback',@(src,evn)export1(src,evn,edit1,fh)    );


uiwait(fh);
    
    function export(~,~,edit1,fh)
%         Result.Decision='Rotate';
        Result.N= str2num(edit1.String );       
        Result.fh= fh;       
        Result.Info='Close' ;
 
        close(fh)
    end

    function export1(~,~,edit1,fh)
%         Result.Decision='Rotate';
        Result.N= str2num(edit1.String );       
        Result.fh= fh;       
        Result.Info='Assign' ;
        close(gcf)
    end

end
