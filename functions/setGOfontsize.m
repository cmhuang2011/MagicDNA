function setGOfontsize( tabH , FS , TypeList )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% TypeList={'Axes','UIControl'}; 

for Ti=1:length(TypeList)
    child_handles_Axes = findall(tabH,'Type',TypeList{Ti});
    for k=1:length(child_handles_Axes)
    child_handles_Axes(k).FontSize=FS ;
    end
end



end

