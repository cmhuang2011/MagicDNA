function SwitchRotate(src,evn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% sdfsf=3;
% evn
if evn.Button==2   % middle key
    
    a = findall(gcf) ;
    b = findall(a,'ToolTipString','Rotate 3D')  ;
    if strcmp(b(1).State , 'off')
        
%         b(1).State='on' ;
            h =rotate3d ;
            
%             uiMode = getuimode(gcf,'Exploration.Rotate3d');
            h.ButtonDownFilter = @myfunction;
             h.Enable='on';
           dsgsdg=3
    else
%         b(1).State='off' ;
             h =rotate3d(gcf) ;
             h.Enable='off';
      
    end
    
%     h =rotate3d(gcf) ;
   
%     h.ActionPostCallback =@myfunction ;
%     sfs=3
end

end


function  myfunction(obj,event_obj)
obj
event_obj
% 
% sdfsf=3
end