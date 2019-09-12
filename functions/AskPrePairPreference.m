function decision = AskPrePairPreference()

decision.pair=2;
decision.NConn=1;


sz = [400 200]; % figure size
% screensize = get(0,'ScreenSize') ;
% xpos = ceil((screensize(3)-sz(1))/2); % center the figure on the
% % screen horizontally
% ypos = ceil((screensize(4)-sz(2))/2); % center the figure on the
% screen vertically

% 300 300 300 110


fh = dialog('units','pixels',...
    'position',[300 300 300 150],...
    'menubar','none',...
    'name','Preference',...
    'numbertitle','off',...
    'resize','off');

movegui(fh,'center')

pop = uicontrol('style','pop','Parent',fh,...
    'units','pixels',...
    'position',[20 10 260 40],...
    'string',{'All Random','Y only','X only'});

str=num2cell(1:7);

pop2 = uicontrol('style','pop','Parent',fh, 'units','pixels',...
    'position',[170 75 100 30], 'string',str);
txtH2 = uicontrol('Style','text','FontSize',12,'Position', [140 110 150 20],....
    'String','# of Connection');

btn = uicontrol('style','push','Parent',fh, 'unit','pix',...
    'position',[50 75 100 30], 'string','OK',...
    'Callback',@(src,evn)export(src,evn,pop,pop2)    );




uiwait(fh);

    function export(~,~,pop,pop2)
        decision.pair=pop.Value;
        decision.NConn=pop2.Value;
        
        close(gcf)
    end

end
