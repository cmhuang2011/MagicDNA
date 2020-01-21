function ColorfulRMSF
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

ax=gca;
% ss=findobj(gca,'Type','scatter') ;    %%--------
ss=findobj(gca,'Tag','plot3k')  ; % search for plot3k objects, only work for RMSF plot
%
SphereR =0.28 ; %0.28
if ~isempty(ss)
    Spwd=pwd;
    oxDNAOBJ= helpCADDOM(5) ;
    cd(oxDNAOBJ.PathName);            [status, msg, msgID] = mkdir('BILDs') ;
    
    
    cd([oxDNAOBJ.PathName 'BILDs'] );
        file2_name=strcat('RMSF' , '.bild')      ;
%     file2_name=strcat('RMSF_MagicDNALogo4_star' , '.bild')      ;  %---------
%     file2_name=strcat('RMSF_Fig5Use_v7' , '.bild')      ;  %---------
    fileID = fopen(file2_name,'w');
    %      fprintf(fileID , '.transparency  %4.2f \n',0.4  )    ;
    for k=1: length(ss)
        %         if k==3
                ColorRGB= ss(k).MarkerFaceColor ;   %-------Notice
        %         else
%                 ColorRGB=ss(k).Color;
        %         end
                fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
        
        %         if k<=2
        %          SphereR=2;
        %         else
        %          SphereR =0.35 ;
        %         end
        
        for Bi = 1 :length(ss(k).XData)
%             ColorRGB= ss(k).CData(Bi,:) ;  %-------Notice
%             fprintf(fileID , '.color %4.2f %4.2f %4.2f \n',ColorRGB  )    ;
%             fprintf(fileID , '.transparency  %8.6f\n',ss(k).UserData(Bi)  )    ;
            fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',ss(k).XData(Bi) ,ss(k).YData(Bi) ,ss(k).ZData(Bi) ,SphereR )    ;
            
        end
        
    end
    fclose(fileID);
    
    cd(Spwd);
    fprintf('Finish printing. \n')
else
    fprintf('the current axes is no RMSF. Can not export to Chimera. \n')
end

end

