function ExportSurfaceToChimera( h_plot, cmap,varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



fprintf('printing surface models for Chimera.\n')
if nargin==2
    file_name='Routing_scaffold' ;
elseif nargin==3
    file_name=varargin{1} ;IndPrint=[];
    
else
    file_name=varargin{1} ; IndPrint=varargin{2} ;
    
end

Radius =0.1; %0.12 default
fileID = fopen([pwd filesep file_name '.bild'],'w');
fprintf(fileID ,'\n' );
fprintf(fileID ,'.comment (If need colors for bundle ) .color  r g b\n' );
%  fprintf(fileID ,'.transparency 0.2\n' );
cc=0;
for  k = 1:length(h_plot)
    % XYZ= [ h_plot.XData(1,:);h_plot.YData(1,:);h_plot.ZData(1,:) ]' ; % for suface object
    
%     if  strcmp(h_plot(k).Visible,'off')  ||  h_plot(k).LineWidth ==0.1
%        continue; 
%     end
    
%     if     h_plot(k).LineWidth ==0.5
%         Radius = 0.2;
% %     elseif h_plot(k).LineWidth ==0.1
% %         continue;  
%     else
%         Radius =0.08;
%     end
    
    XYZ = [ h_plot(k).XData(1,:);h_plot(k).YData(1,:);h_plot(k).ZData(1,:) ]' ;
    Radius = 0.6 ;
    %     XYZ(end+1,:)=XYZ(1,:) ; %CM
    %    RRcolor = rand(1,3) ;
    %      fprintf(fileID ,'.color %8.6f %8.6f %8.6f\n',  RRcolor );
    %      fprintf(fileID , '.marker  %8.6f %8.6f %8.6f\n',XYZ(1,:)')    ;
%     fprintf(fileID ,'.color %8.6f %8.6f %8.6f\n',  h_plot(k).Color );
    
%     if strcmp(get(h_plot(k), 'Type') ,'surface')
%         %         fprintf(fileID ,'.color %8.6f %8.6f %8.6f\n',  h_plot(k).Color );
%         %          fprintf(fileID ,'.color %8.6f %8.6f %8.6f\n', [0.2,0.8,0.3]);
%         
%         %         Radius =0.6 ;
%     else
%         %         Radius =0.1;
%         %          fprintf(fileID ,'.color %8.6f %8.6f %8.6f\n',  rand(1,3) );
%     end
    
    
    for Bi = 1:size(XYZ ,1)-1   %
%         if Bi>200
%             Radius = 0.2;
% %         elseif Radius==0.08 
% %             continue 
%         end
        s= Bi/ size(XYZ ,1)  ;
%         RefX = ([0:size(cmap ,1)-1]')/(size(cmap ,1)+1) ;
        RefX=linspace(0 , 1 , size(cmap ,1)) ;
        
        Vcolor = interp1(RefX, cmap , s ) ;
        
        if norm(diff(XYZ(Bi:Bi+1,:)) ,2)>1e-5
            
            
             fprintf(fileID ,'.color %8.6f %8.6f %8.6f\n', Vcolor);
            fprintf(fileID , '.cylinder %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f open\n',XYZ(Bi,:)',XYZ(Bi+1,:)',Radius )    ;
        end
        fprintf(fileID , '.sphere %8.6f %8.6f %8.6f %8.6f \n',XYZ(Bi,:)',Radius )    ;
    end
    fprintf(fileID , '.sphere %8.6f %8.6f %8.6f %8.6f \n',XYZ(Bi+1,:)',Radius )    ;
    
end

fclose(fileID);





end

