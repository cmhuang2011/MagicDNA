function [ result ] = FindBelongRange( arr1,arr2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

N=length(arr1)/2;
L2=arr2(1);R2=arr2(2);
result=[];
for i=1:N
   Left1=arr1(2*i-1); Right1=arr1(2*i);  
  heights=sort([Left1 Right1 L2 R2 ]);
  if max(heights)-min(heights)<= abs(L2-R2)+abs(Left1-Right1) 
   result=union(result,i);
  end 
end
% result shows intersect with which segments

end

