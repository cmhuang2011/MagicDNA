function tabH = gctab
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ax=gca;
if strcmp(ax.Parent.Type,'uitab')
    tabH=ax.Parent ;
else
   error('Parent of gac is not a uitab ') 
end


end

