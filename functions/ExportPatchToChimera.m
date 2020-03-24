function ExportPatchToChimera( h_plot,varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



fprintf('printing patch models for Chimera.\n')
if nargin==0
    file_name='AlphaShape_' ;
else
    file_name=varargin{1} ;
end

Radius =0.15;
fileID = fopen([pwd filesep file_name '.bild'],'w');
fprintf(fileID ,'\n' );
fprintf(fileID ,'.comment (If need colors for bundle ) .color  r g b\n' );
%  fprintf(fileID ,'.transparency 0.5\n' );
for  k = 1:length(h_plot)
    % XYZ= [ h_plot.XData(1,:);h_plot.YData(1,:);h_plot.ZData(1,:) ]' ; % for suface object
    
    if strcmp(get(h_plot(k), 'Type') ,'patch')
%         fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  h_plot(k).FaceColor );
%          fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n', [0.2,0.8,0.3]);
           XYZ =  h_plot(k).Vertices ;
           Faces =  h_plot(k).Faces ;
%         Radius =0.2;
    else
        continue;
%         Radius =0.1;
%          fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  rand(1,3) );
    end
    
    for tri_i = 1:size(Faces ,1)
        Inds= Faces(tri_i ,:) ;
        
        xyzAll = XYZ(Inds,:)  ;
         OneArr =reshape(xyzAll', numel(xyzAll),1) ;
        
        fprintf(fileID , '.polygon  %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',OneArr )    ;
        
%         fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f open\n',XYZ(Bi,:)',XYZ(Bi+1,:)',Radius )    ;
        
%         fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',XYZ(Bi,:)',Radius )    ;
    end
    
    
end

fclose(fileID);
fprintf('Finish printing.\n')




end

