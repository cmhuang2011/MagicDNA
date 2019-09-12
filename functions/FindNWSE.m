function [ VV ] = FindNWSE(PosCyl,NegCyl,coor)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

Vec=zeros(1,2);
Vec(1)=int8(coor(NegCyl,1)-coor(PosCyl,1));
Vec(2)=int8(coor(NegCyl,2)-coor(PosCyl,2));
Vec=Vec/norm(Vec);
VV=[];
if (Vec(2)==1) &&  (Vec(1)==0)         %North
   VV=1     ;
end
if (Vec(2)==0) &&  (Vec(1)==1)         %East
   VV=4     ;
end
if (Vec(2)==-1) &&  (Vec(1)==0)         %South
   VV=3     ;
end
if (Vec(2)==0) &&  (Vec(1)==-1)         %West
   VV=2     ;
end

end

