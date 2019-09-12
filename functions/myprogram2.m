function myprogram2
    f = figure('WindowStyle','normal');
    c = uicontextmenu(f);

    % Assign the uicontextmenu to the figure
    f.UIContextMenu = c;

    % Create child menu of the uicontextmenu
    topmenu = uimenu('Parent',c,'Label','Change Color');
    
    % Create submenu items
    m1 = uimenu('Parent',topmenu,'Label','Red','Callback',@changecolor);
    m2 = uimenu('Parent',topmenu,'Label','Green','Callback',@changecolor);

    function changecolor(source,callbackdata)
        switch source.Label
            case 'Red'
                f.Color = [1.0 0.80 0.80];
            case 'Green'
                f.Color = [0.80 1.0 0.80];
        end
    end
end