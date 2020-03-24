function ExportTrisToChimera( Vertices,Tri,varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Vertices: Nx3 ,  Tri: Kx4, indexes in Vertices to form tetrahedron 

fprintf('printing triangulation for Chimera.\n')
if nargin==0
    file_name='Trianulations_' ;
else
    file_name=varargin{1} ;
    
end

fileID = fopen([pwd filesep file_name '.bild'],'w');
fprintf(fileID ,'\n' );
fprintf(fileID ,'.comment (If need colors for bundle ) .color  r g b\n' );
%  fprintf(fileID ,'.transparency 0.5\n' );
% for  k = 1:size(Vertices ,1 )
    % XYZ= [ h_plot.XData(1,:);h_plot.YData(1,:);h_plot.ZData(1,:) ]' ; % for suface object
    
%     if strcmp(get(h_plot(k), 'Type') ,'patch')
% %         fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  h_plot(k).FaceColor );
% %          fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n', [0.2,0.8,0.3]);
% %            XYZ =  h_plot(k).Vertices ;
% %            Faces =  h_plot(k).Faces ;
% %         Radius =0.2;
%     else
%         continue;
% %         Radius =0.1;
% %          fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  rand(1,3) );
%     end
    
    for tri_i = 1:size(Tri ,1)
        
        Inds_vertices= Tri(tri_i ,:) ;
        
        xyzAll = Vertices(Inds_vertices,:)  ; 
         
       OneArr=[xyzAll(1,:),xyzAll(2,:),xyzAll(3,:),xyzAll(4,:),xyzAll(2,:)] ;
        
%         OneArr =reshape(xyzAll', numel(xyzAll),1) ;
        
        fprintf(fileID , '.polygon  %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n',OneArr )    ;
        
    end
    
    
% end

fclose(fileID);
fprintf('Finish printing.\n')




end

