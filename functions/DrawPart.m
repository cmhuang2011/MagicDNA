function  DrawPart( Part,color,YText )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Part
% figure;hold on;
if isempty(color)
    color=[1 0 0];
end
N=size(Part.CylInplanePosition,1);
for i=1:N 
   c=plot3( [Part.CylInplanePosition(i,1) Part.CylInplanePosition(i,1)],[Part.CylInplanePosition(i,2) Part.CylInplanePosition(i,2)],[Part.Z1(i) Part.Z2(i)],'.-r','Linewidth',4,'LineStyle','-');  
   if ismember(i, Part.AGroup)
   c.Color=color;
   else
    c.Color=[0 1  0];   
   end
   if YText==1
   text(Part.CylInplanePosition(i,1),Part.CylInplanePosition(i,2),Part.Z1(i)-2,num2str(i),'FontSize',10);
   end
end
axis auto ;
end

