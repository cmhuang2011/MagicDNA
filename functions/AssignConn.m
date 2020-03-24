function choice = AssignConn(B1,hyperB,AssFigH,patchH,lineIndex)


hyperB.TakeSeveralV3{lineIndex}=[];  %reset forcedConnection in this edge
% d = figure('Position',[300 300 800 250],'Name','Specify Connect Options');
d = dialog('Position',[300 300 800 250],'Name','Specify Connect Options');

movegui(d,'center') ;
hp1 = uipanel('Parent',d,'Title','Auto Searching','FontSize',12,...
    'Position',[0.05  0.05   0.4 0.4]);
hp2 = uipanel('Parent',d,'Title','Manually','FontSize',12,...
    'Position',[0.55  0.05   0.4 0.4]);



txt = uicontrol('Parent',d,'Style','text','Position',[180 140 250 50],...
    'String','Enter Number of Connections','FontSize',14);
ed = uicontrol('Parent',d,'style','popup','position',[450 130 150 40],...
    'string',{'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'17';'29 ';'?' },'FontSize',12);
% here to change number of connections options, integer only



str1=strcat('Bundle-  ', num2str(B1(1)));
txt1 = uicontrol('Parent',hp1,...
    'Style','text','Units','normalized','FontSize',12,...
    'Position',[0.05 0.5 0.4 0.3],...
    'String',str1);

popup1 = uicontrol('Parent',hp1,...
    'Style','popup','Units','normalized',...
    'Position',[0.05 0.2 0.4 0.25],'FontSize',12,...
    'String',{'Both';'Side only';'End only'},...
    'Callback',@popup_callback1);
%            popup_callback1

str2=strcat('Bundle-  ', num2str(B1(2)));
txt2 = uicontrol('Parent',hp1,...
    'Style','text','Units','normalized','FontSize',12,...
    'Position',[0.55 0.5 0.4 0.3],...
    'String',str2);

popup2 = uicontrol('Parent',hp1,...
    'Style','popup','Units','normalized',...
    'Position',[0.55 0.2 0.4 0.25],'FontSize',12,...
    'String',{'Both';'Side only';'End only'},...
    'Callback',@popup_callback2);

btn = uicontrol('Parent',d,...
    'Units','normalized','Position',[0.45 0.2 0.1 0.1],...
    'String','OK','FontSize',12,...
    'Callback',{@(src,evn)Exportresult(src,evn,ed)});


btnCallManualInput = uicontrol('Parent',hp2,'Units','normalized',...
    'Position',[0.1 0.1 0.7 0.7],...
    'String','Manually Assign','FontSize',12,...
    'Callback',{@(src,evn)CallManual(src,evn,hyperB,AssFigH,patchH,B1,ed,lineIndex)});


choice.Bundle1 = 1;  hyperB.choice.Bundle1=1;
choice.Bundle2 = 1;  hyperB.choice.Bundle2=1;
choice.NumberOfConnect=1;  hyperB.choice.NumberOfConnect=1;
choice.UseMethod='Auto';    hyperB.choice.UseMethod='Auto';
%  {@(src,evn)SortDist(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ)}
%
% setGOfontsize( gcf , 10 , {'UIControl'} )  % set fontsize for uicontrol in this tab
% 
% drawnow;
% print('-r100',gcf,'Image3_r100','-dpng');
% Wait for d to close before running to completion

uiwait(d);  

    function popup_callback1(popup,event)
        idx = popup.Value;
        popup_items = popup.String;
        choice.Bundle1 =idx;
        hyperB.choice.Bundle1=idx;  % output in class's temp space
    end
    function popup_callback2(popup,event)
        idx = popup.Value;
        popup_items = popup.String;
        choice.Bundle2 = idx;
        hyperB.choice.Bundle2=idx;  % output in class's temp space
    end
    function Exportresult(src,evn,ed)
        %         choice.NumberOfConnect=num2str(ed.String);
        if ~strcmp(ed.String{ed.Value}, '?')
        choice.NumberOfConnect=str2double(ed.String{ed.Value});
        hyperB.choice.NumberOfConnect=str2double(ed.String{ed.Value});  % output in class's temp space
        else
            answer = inputdlg('Enter a number:',  'Sample', [1 20]) ;
%             str2double(answer)
         choice.NumberOfConnect=str2double(answer);
         hyperB.choice.NumberOfConnect=str2double(answer);  % output in class's temp space
         
%          str2double(answer)
        end
        
%         hyperB.choice.NumberOfConnect=str2double(ed.String{ed.Value});  % output in class's temp space
%         choice;
        delete(gcf)
    end

end