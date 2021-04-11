function ColorfulRMSF
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

ax=gca; fH =gcf;
% ss=findobj(gca,'Type','scatter') ;    %%--------
ss=findobj(gca,'Tag','plot3k')  ; % search for plot3k objects, only work for RMSF plot
%
btn_LoadTraj= findobj(gcf,'Tag','btn_oxDNATraj1') ;   % get oxDNA object handle for topology
oxDNA_ex =  btn_LoadTraj.UserData.oxDNA_ex ;

% %-----shaded background in Chimera
ExportShade=0 ;
if ExportShade==1
tic;
fprintf('Calculating boundary patches........\n')
AlltrajXYZ= ax.UserData.AlltrajXYZ ;
shp = alphaShape(AlltrajXYZ(:,1),AlltrajXYZ(:,2),AlltrajXYZ(:,3)) ;
figure(25235);clf ;
pshp = plot(shp) ;
TargetFaces= 5000;
R=TargetFaces/size(pshp.Faces ,1);
reducepatch(pshp,R) ;
ExportPatchToChimera( pshp,'ShadedPatch' );
fprintf('Finish Calculating boundary nodes........\n')
T_calBound = toc
end
%---------------

figure(fH) ;

SphereR =0.28 ; %0.28
if ~isempty(ss)
    Spwd=pwd;
    oxDNAOBJ= helpCADDOM(5) ;
    cd(oxDNAOBJ.PathName);            [status, msg, msgID] = mkdir('BILDs') ;
    
    
    cd([oxDNAOBJ.PathName 'BILDs'] );
    file2_name=strcat('RMSF' , '.bild')      ;
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
    
    fprintf('Finish printing. \n') ;
    %-------------
    CColor = get(gca,'colororder')  ;
    Radius =SphereR*0.8 ; %0.12 default
    
    
    file3_name=strcat('MeanConf' , '.bild')      ;
    fileID2 = fopen(file3_name,'w');
    
    MeanConf = oxDNA_ex.MeanConf ;
    for k = 1: max(oxDNA_ex.Strand)
        BaseInd = oxDNA_ex.Strand==k ;
%         if k==1 % scaffold strand
        if sum(BaseInd) >1000 % scaffold strand
%             
%                     IndRemove = union(1:981, 7462:7560 );
%                     Ind = setdiff(1:7560, IndRemove);
%                     BaseInd = find(oxDNA_ex.Strand==k ) ;
%                     BaseInd = BaseInd(Ind) ;   % hard
            
            fprintf(fileID2 ,'.color %8.6f %8.6f %8.6f\n',  CColor(1,:) );
        else
%             if sum(BaseInd) == 36 % extended
%                BundleColor=[248,144,31]/255 ;   
%             elseif  sum(BaseInd) == 26 % shrink
%               BundleColor=[0.3,0.3,0.3] ;    
%             else
% %                 continue ; % hard
%               BundleColor=[0.5,0,0] ;    
%             end
                
            BundleColor=[0.8,0,0] ;

            fprintf(fileID2 , '.color %8.6f %8.6f %8.6f \n', BundleColor+0.1*(0.5-rand(1,3))  )    ;
%             fprintf(fileID2 ,'.color %8.6f %8.6f %8.6f\n',  CColor(randi(5)+1,:) );
        end
        XYZ = MeanConf(BaseInd ,1:3) ;
        
        for Bi = 1:size(XYZ ,1)-1   %
            if norm(diff(XYZ(Bi:Bi+1,:)) ,2)>1e-5
                fprintf(fileID2 , '.cylinder %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f open\n',XYZ(Bi,:)',XYZ(Bi+1,:)',Radius )    ;
            end
            fprintf(fileID2 , '.sphere %8.6f %8.6f %8.6f %8.6f \n',XYZ(Bi,:)',Radius )    ;
        end
        fprintf(fileID2 , '.sphere %8.6f %8.6f %8.6f %8.6f \n',XYZ(Bi+1,:)',Radius )    ;
    end
    fclose(fileID2);
    fprintf('Finish printing ribbon on backbones. \n') ;
    
    
    %     sdfs=3 ;
    if ExportShade==1
    ExportPatchToChimera( pshp,'ShadedPatch' );
    end
    %------------------
    cd(Spwd);
else
    fprintf('the current axes is no RMSF. Can not export to Chimera. \n')
end

end

