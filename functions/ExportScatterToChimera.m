function ExportScatterToChimera( h_plot,varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



fprintf('printing surface models for Chimera.\n')
if nargin==1
    file_name='Scatter_fromMATLAB' ;
else
    file_name=varargin{1} ;
end

% Radius =0.25;
fileID = fopen([pwd filesep file_name '.bild'],'w');
fprintf(fileID ,'\n' );
fprintf(fileID ,'.comment (If need colors for bundle ) .color  r g b\n' );
%  fprintf(fileID ,'.transparency 0.2\n' );
for  k = 1:length(h_plot)
    % XYZ= [ h_plot.XData(1,:);h_plot.YData(1,:);h_plot.ZData(1,:) ]' ; % for suface object
    
    XYZ = [ h_plot(k).XData(1,:);h_plot(k).YData(1,:);h_plot(k).ZData(1,:) ]' ;
%      fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  [0 0 0] );
%      fprintf(fileID , '.marker  %4.2f %4.2f %4.2f\n',XYZ(1,:)')    ;
    
    if strcmp(get(h_plot(k), 'Type') ,'scatter')
%         fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  h_plot(k).CData );
        %          fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n', [0.2,0.8,0.3]);
        
        switch h_plot(k).SizeData
            case 96
                Radius = 0.3 ;
            case 36
                Radius = 0.1 ;
            case 58
                Radius = 0.5 ;
                
            case 6
                Radius = 0.1 ;
                
            otherwise
                Radius = 0.05 ;
        end
    else
        Radius =0.1;
%          fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  rand(1,3) );
    end
    Radius=Radius*3 ;
    
    Radius=0.2;
    for Bi = 1:size(XYZ ,1)-1   %
%         fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f open\n',XYZ(Bi,:)',XYZ(Bi+1,:)',Radius )    ;
        
        fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',XYZ(Bi,:)',Radius )    ;
    end
     fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',XYZ(Bi+1,:)',Radius )    ;
    
end

fclose(fileID);





end

