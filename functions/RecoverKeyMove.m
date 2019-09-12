function RecoverKeyMove(src,evn,ccase)
%   Detailed explanation goes here
% ------------this will be very usefully to assign different
% keypressfunction so user can contorl orientation in difference tabs

fH=gcf;
try
    switch ccase
        case 1
            fH.KeyPressFcn=fH.UserData.saveKeyMove ;
        case 2
            fH.KeyPressFcn=fH.UserData.saveKeyMove2 ;
%         case 3
      
    end
catch
end
end

