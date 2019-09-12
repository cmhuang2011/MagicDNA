function  LineScatterObjectToChimera
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


ax=gca;
ss=findobj(gca,'Type' ,'Line') ; 
ss2=findobj(gca,'Type' ,'Scatter') ; 

SphereR=0.2;
if ~isempty(ss)
    [status, msg, msgID] = mkdir('BILDs') ;

    
    file2_name=strcat('MATLAB_lines' , '.bild')      ;
    fileID = fopen(file2_name,'w');
%      fprintf(fileID , '.transparency  %4.2f \n',0.4  )    ;
    for k=1:length(ss)
        LineWidths=ss(k).LineWidth/1.5 ;
        ColorRGB= ss(k).Color ; 
        fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
        
        if isempty(ss(k).ZData)
        XYZ= [ss(k).XData ; ss(k).YData ;zeros(size(ss(k).YData))]' ;
        else
        XYZ= [ss(k).XData ; ss(k).YData ;ss(k).ZData]' ;
        end
        indISnan = isnan(XYZ(:,1)) ;
         XYZ=XYZ(~indISnan ,:)  ; XYZ=XYZ +0.05*rand(size(XYZ)) ;
        for Bi = 1:1 :size(XYZ,1)-1
%         fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',ss(k).XData(Bi) ,ss(k).YData(Bi) ,ss(k).ZData(Bi) ,SphereR )    ;
if norm(XYZ(Bi,1:3)-XYZ(Bi+1,1:3)) >0.2
         fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',XYZ(Bi,1:3),XYZ(Bi+1,1:3) ,LineWidths )    ;
end
%          Center=[scaf.pScaf_center{1}.XData(k),scaf.pScaf_center{1}.YData(k),scaf.pScaf_center{1}.ZData(k)] ;
%         Backbone=[scaf.pScaf2{1}.XData(k),scaf.pScaf2{1}.YData(k),scaf.pScaf2{1}.ZData(k)] ;
%         fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',Center,Backbone,0.1 )    ;             
      
        end
        
    end
    
%     for k2 = 1:length(ss2)
%         ColorRGB= ss2(k2).CData(1,:) ; 
%         fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
%      for Bi = 1 :length(ss2(k2).XData)
%              fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',ss2(k2).XData(Bi) ,ss2(k2).YData(Bi) ,ss2(k2).ZData(Bi),1.4 )    ;
%      
%      end
%         
%     end
        
    fclose(fileID);
    
%     cd(Spwd)
else
    fprintf('the current axes is no lines. Can not export to Chimera. \n')
end




end

