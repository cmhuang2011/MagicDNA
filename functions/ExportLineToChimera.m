function ExportLineToChimera( h_plot,varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here



fprintf('printing surface models for Chimera.\n')
if nargin==1
    file_name='Routing_stap' ;
elseif nargin==2
    file_name=varargin{1} ;IndPrint=[];
   
else
     file_name=varargin{1} ; IndPrint=varargin{2} ;
   
end

Radius =0.2;
fileID = fopen([pwd filesep file_name '.bild'],'w');
fprintf(fileID ,'\n' );
fprintf(fileID ,'.comment (If need colors for bundle ) .color  r g b\n' );
%  fprintf(fileID ,'.transparency 0.2\n' ); 
cc=0;
for  k = 1:length(h_plot)
    % XYZ= [ h_plot.XData(1,:);h_plot.YData(1,:);h_plot.ZData(1,:) ]' ; % for suface object
    
    
    
    if sum( h_plot(k).Color ==[1 0 0])~=3 
%         dsf=2
    elseif ismember(h_plot(k).UserData.Ind, IndPrint )
        continue
    else
        continue
    end
    
    
%     if ismember(h_plot(k).UserData.Ind, IndPrint ) && sum( h_plot(k).Color ==[1 0 0])==3 % hard code
% %         fprintf('slip some lines \n')
% %     continue  ; %ignore some lines.
%  cc=cc+1
%     else
%         continue  ; %ignore some lines.
%        
%     end
%     
    XYZ = [ h_plot(k).XData(1,:);h_plot(k).YData(1,:);h_plot(k).ZData(1,:) ]' ;
    
%     XYZ(end+1,:)=XYZ(1,:) ; %CM
%      fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  [0 0 0] );
%      fprintf(fileID , '.marker  %4.2f %4.2f %4.2f\n',XYZ(1,:)')    ;
    
    if strcmp(get(h_plot(k), 'Type') ,'line')
        fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  h_plot(k).Color );
%          fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n', [0.2,0.8,0.3]);
       
%         Radius =0.6 ;
    else
%         Radius =0.1;
%          fprintf(fileID ,'.color %4.2f %4.2f %4.2f\n',  rand(1,3) );
    end
    
    
    for Bi = 1:size(XYZ ,1)-1   %
        fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f open\n',XYZ(Bi,:)',XYZ(Bi+1,:)',Radius )    ;
        
        fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',XYZ(Bi,:)',Radius )    ;
    end
     fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',XYZ(Bi+1,:)',Radius )    ;
    
end

fclose(fileID);





end

