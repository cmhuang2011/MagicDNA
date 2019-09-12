function main_ChangeSize(src,evn)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% src=gcf;
% src.Visible='off';
altBtn = findobj(gcf,'Style', 'pushbutton') ;
% fliplr(fH.Position(1,3:4))
% drawnow
% src.Position(1,3:4);
% tic
for k=1:length(altBtn)
    if ~isempty(altBtn(k).CData)       
         altBtn(k).Visible='off';

        altBtn(k).Units='pixels';
%         fliplr(altBtn(k).Position(1,3:4) );
% drawnow
        A=imresize(altBtn(k).UserData.OriCData ,fliplr(altBtn(k).Position(1,3:4)))  ;
        altBtn(k).Units='normalized';
        altBtn(k).CData=A;
%         size(A)
         altBtn(k).Visible='on';
    end
end

% toc
% src.Visible='on';

end

