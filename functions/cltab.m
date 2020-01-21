function  cltab
% clean tab's axes graphics 

%   Detailed explanation goes here

ax=gca;
if strcmp(ax.Parent.Type,'uitab')
%     tabH=ax.Parent ;
    siblingAxes=findobj(ax.Parent,'Type','Axes') ;
    for k=1:length(siblingAxes)
       axes( siblingAxes(k));
       cla;
    end
    
%     tables =findobj(ax.Parent,'Type','Table') ;
    
else
   error('Parent of gca is not a uitab ') 
end



end

