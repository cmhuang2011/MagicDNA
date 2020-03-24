function Result = AskLocalSnapOrRotationInterval()

Result.Decision='Snap';
Result.RotationInterval=5;


% sz = [400 200]; % figure size
% screensize = get(0,'ScreenSize') ;
% xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the
% % screen horizontally
% ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the
% screen vertically

% 300 300 300 110


fh = dialog('units','pixels',...
    'position',[300 300 300 150],...
    'menubar','none',...
    'name',' Adjust local',...
    'numbertitle','off',...
    'resize','off');

movegui(fh,'center')

% pop = uicontrol('style','pop','Parent',fh,...
%     'units','pixels',...
%     'position',[20 10 260 40],...
%     'string',{'All Random','Y only','X only'});

str=num2cell(1:7);
edit1 = uicontrol('style','edit','Parent',fh, 'units','pixels','String',5,...
    'position',[170 75 100 30]);
txtH2 = uicontrol('Style','text','Parent',fh,'FontSize',12,'Position', [140 110 150 20],....
    'String','Rotate along z by:(deg)');

btn1 = uicontrol('style','push','Parent',fh, 'unit','pix',...
    'position',[50 35 100 60], 'string','Snap x ',...
    'Callback',@(src,evn)btn1Call(src,evn)    );
% { }
btn2 = uicontrol('style','push','Parent',fh, 'unit','pix',...
    'position',[170 20 100 30], 'string','OK',...
    'Callback',@(src,evn)export(src,evn,edit1)    );



uiwait(fh);
    function btn1Call(~,~)
        Result.RotationInterval=nan ;
        fprintf('Snape local x direction to Global +-XYZ by applying the closest rotation matrixes within +-45 degrees.\n')
        close(gcf)
    end

    function export(~,~,edit1)
        Result.Decision='Rotate';
        Result.RotationInterval= str2num(edit1.String );
        
        close(gcf)
    end

end
