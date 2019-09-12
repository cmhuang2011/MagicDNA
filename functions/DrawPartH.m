function [plotH,PlotXYZ]= DrawPartH( Part,color,YText,aH ,bundlei)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Part
% figure;hold on;
if isempty(color)
    color=[1 0 0];
end
N=size(Part.CylInplanePosition,1);
% plotH=cell(1,N);
% PlotXYZ=zeros(N,6);
PlotXYZ2=Part.CylinderXYZGlobal;
PlotXYZ=[PlotXYZ2(:,1),PlotXYZ2(:,4),PlotXYZ2(:,2),PlotXYZ2(:,5),PlotXYZ2(:,3),PlotXYZ2(:,6)];
    for i=1:N 
%        PlotXYZ(i,:)= [[Part.CylInplanePosition(i,1) Part.CylInplanePosition(i,1)],[Part.CylInplanePosition(i,2) Part.CylInplanePosition(i,2)],[Part.Z1(i) Part.Z2(i)]];
%        plotH{i}=plot3( [Part.CylInplanePosition(i,1) Part.CylInplanePosition(i,1)],[Part.CylInplanePosition(i,2) Part.CylInplanePosition(i,2)],[Part.Z1(i) Part.Z2(i)],'.-r','Linewidth',4,'LineStyle','-');  
%        if ismember(i, Part.AGroup)
% %        plotH{i}.Color=color;
%        else
% %         plotH{i}.Color=[0 1  0];   
%        end
       if YText==1
       text(Part.CylInplanePosition(i,1),Part.CylInplanePosition(i,2),Part.Z1(i)-2,num2str(i),'FontSize',10 );
       end
    end
    
    
plotH=plot3(aH, PlotXYZ(:,1:2)',PlotXYZ(:,3:4)',PlotXYZ(:,5:6)' ,'.-','Linewidth',4.5,'LineStyle','-','color',[0.2 0.4 0.2]); 
% plotH.UserData.belongBundle=bundlei*ones(4,1);
for cc=1:length(plotH)
  plotH(cc).UserData= [ bundlei,cc];
end

end

