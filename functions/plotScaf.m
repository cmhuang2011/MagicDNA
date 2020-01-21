
% f41H=figure(41);clf ;
% subplot(1,2,1);
ax= findobj(0,'Tag','MechScaffold3D');
axes(ax); cltab ; axes(ax);

hold on; axis  equal;

AllGXYZ=cell(1,length(obj.containBundle));
for k=1:length(obj.containBundle)
AllGXYZ{k}=obj.containBundle{k}.CylinderXYZGlobal ;
C1Center=mean([AllGXYZ{k}(1,1:3);AllGXYZ{k}(1,1:3);AllGXYZ{k}(1,4:6)]);
f41H.UserData.textH{k}=text( C1Center(1),C1Center(2),C1Center(3)+5, num2str(k),'FontSize',36);
end 

 SaveGHelix=cell(1,length(obj.containBundle));
for Bundlei=1:length(obj.containBundle)
 QQWithPM10Bases =obj.containBundle{Bundlei}.HelixXYZG;     
SaveGHelix{Bundlei}= QQWithPM10Bases;
end

 SaveGHelixBVec=cell(1,length(obj.containBundle)); %in old, before-simul coordinate
for Bundlei=1:length(obj.containBundle)
 QQWithPM10Bases =obj.containBundle{Bundlei}.HelixXYZGBVec;     
SaveGHelixBVec{Bundlei}= QQWithPM10Bases;
 end
 

SacfR=obj.ScafRouting ;
SacfR;
plotXYZ=zeros(size(SacfR));
for k=1:size(SacfR,1)
   bundle=SacfR(k,1);  Cyl=SacfR(k,2);
   alpha=SacfR(k,3)- obj.containBundle{bundle}.Zbase1(Cyl);
   beta=obj.containBundle{bundle}.Zbase2(Cyl)-SacfR(k,3);
   P= AllGXYZ{bundle}(Cyl,1:3);
   Q=AllGXYZ{bundle}(Cyl,4:6);
   XYZ=(beta*P + alpha*Q )/(alpha+beta);
   plotXYZ(k,:)=XYZ;   
end

SSplotXYZ=size(plotXYZ) ;
% plot3(plotXYZ(:,1), plotXYZ(:,2), plotXYZ(:,3) )
x = plotXYZ(:,1)';
y = plotXYZ(:,2)';
z = plotXYZ(:,3)';
col = (1:size(plotXYZ,1))*1000;  % This is the color
surface([x;x],[y;y],[z;z],[col;col], 'facecol','no', 'edgecol','interp', 'linew',2);
%     DSh=plot3(XYZdata(:,1),-XYZdata(:,2),XYZdata(:,3),'Linewidth',2);   %draw route
scatter3(plotXYZ(1,1),plotXYZ(1,2),plotXYZ(1,3),'s','filled','SizeData',100);   %mark head
scatter3(plotXYZ(end,1),plotXYZ(end,2),plotXYZ(end,3),'d','filled','SizeData',100);  %mark tail
%  xlim auto ; ylim auto; zlim auto;
axis equal;grid on;  xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;


  %------------------------------------------------------------------  