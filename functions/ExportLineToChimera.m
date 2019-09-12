function ExportLineToChimera( h_plot )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



fprintf('printing surface models for Chimera.\n')
file_name='Routing' ;   Radius =0.15;
fileID = fopen([pwd filesep file_name '.bild'],'w');
fprintf(fileID ,'\n' );
fprintf(fileID ,'.comment (If need colors for bundle ) .color  r g b\n' );

for  k = 1:length(h_plot)
% XYZ= [ h_plot.XData(1,:);h_plot.YData(1,:);h_plot.ZData(1,:) ]' ; % for suface object

XYZ = [ h_plot(k).XData(1,:);h_plot(k).YData(1,:);h_plot(k).ZData(1,:) ]' ;
if strcmp(get(h_plot(k), 'Type') ,'line')
fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  h_plot(k).Color );
Radius =0.13;
else
Radius =0.1;
    
end

for Bi = 1:size(XYZ ,1)-1
      fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',XYZ(Bi,:)',XYZ(Bi+1,:)',Radius )    ;
   
end


end

fclose(fileID);





end

