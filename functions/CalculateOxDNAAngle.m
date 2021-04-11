function CalculateOxDNAAngle()
% oo=gco;  button for load trajectory
oo=findobj(gcf,'String','Load Initail+Traj') ;
oxDNA_ex= oo.UserData.oxDNA_ex  ;

if isempty(oxDNA_ex.MeanConf)
    btn=findobj(gcf,'String','RMSD RMSF') ;
    feval(get(btn,'Callback'),btn,[]);
end
% load(strcat(oxDNA_ex.PathName, 'BM.mat'))
BM=oxDNA_ex.BM ;
d = figure('Position',[300 300 800 350],'Name','Specify Bundle indexes for calculating angles in trajectory','numbertitle','off');
% d = dialog('Position',[300 300 800 350],'Name','Specify Bundle indexes for calculating angles in trajectory');

movegui(d,'center') ;

strSS=strcat('Detect the base vector in the folder. Max =', num2str(max(BM))) ;
txt = uicontrol('Parent',d,'Style','text','Units','normalized','Position',[0.1750 0.8 0.8 0.15],...
    'String',strSS,'FontSize',14);


Str={'--'};
for k=1:max(BM);  Str{k}=strcat('Bundle  ',num2str(k));    end

popup1 = uicontrol('Parent',d,...
    'Style','popup','Units','normalized',...
    'Position',[ 0.1750 0.5 0.3 0.25],'FontSize',12,...
    'String',Str );

popup2 = uicontrol('Parent',d,...
    'Style','popup','Units','normalized',...
    'Position',[0.525 0.5 0.3 0.25],'FontSize',12,...
    'String',Str );

txt2 = uicontrol('Parent',d,'Style','text','Units','normalized','Position',[0.175 0.45 0.55 0.1],...
    'String','ratio of bases arround centers to be considered','FontSize',14);

sld1 = uicontrol('Style', 'slider','Parent',d,'Units','normalized',....
    'Min',0.6,'Max',1,'Value',1,'Position', [0.175 0.25 0.65 0.1]);     % scale


txt3 = uicontrol('Parent',d,'Style','text','Units','normalized','Position',[0.83 0.25 0.1 0.1],...
    'String','1','FontSize',14);
sld1.Callback= @(src,evn)sldscale1(src,evn,txt3) ;

checkH = uicontrol(d,'Style', 'checkbox','String', 'Refresh Angular plot','Unit','normalized','Position', [0.75 0.45 0.2 0.1] ,'FontSize',12); % [0.7 0.85 0.07 0.05]


btn = uicontrol('Parent',d,...
    'Units','normalized','Position',[0.45 0.1 0.1 0.1],...
    'String','OK','FontSize',12,... ;%...
    'Callback',{@(src,evn)CalAngles(src,evn,d,popup1,popup2,sld1, oxDNA_ex,checkH  )});  %



    function sldscale1(src,evn,txt3)
        txt3.String=num2str(src.Value) ;
    end



    function CalAngles(src,evn,d,popup1,popup2,sld1, oxDNA_ex,checkH)
        %
        Bi = popup1.Value ;
        Bj = popup2.Value ;
        
        if Bi==Bj
            opts = struct('WindowStyle','modal',...
                'Interpreter','tex');
            f = warndlg('\fontsize{15}Please specify two different bundles.','Warning',opts);
            return
        end
        
        % BM=oxDNA_ex.BM ;
        RatioBase = sld1.Value ;
        MeanConf = oxDNA_ex.MeanConf ;
        
        %--slider to filter ends to Ind_Bi Ind_Bj
        
        Ind_Bi = find(oxDNA_ex.BM==Bi) ; Ind_Bj = find(oxDNA_ex.BM==Bj) ;
        Ind_BiO =Ind_Bi ;        Ind_BjO =Ind_Bj ;  % Save for in case
        Center_Bi = mean( MeanConf(Ind_Bi,1:3) ) ;
        Center_Bj = mean( MeanConf(Ind_Bj,1:3) ) ;
        
        dXYZ_Bi = MeanConf(Ind_Bi,1:3)-  ones(length(Ind_Bi),1)*Center_Bi ;
        dXYZ_Bj = MeanConf(Ind_Bj,1:3)-  ones(length(Ind_Bj),1)*Center_Bj ;
        
        %         d_i= vecnorm(dXYZ_Bi ,2, 2)  ;   % not available for
        %         2017a
        d_i= dXYZ_Bi(:,1).^2+dXYZ_Bi(:,2).^2+dXYZ_Bi(:,3).^2 ;
        %         d_i
        [~,bbi ] =sort(d_i) ; NIndGetI=  round(length(d_i)*RatioBase) ;
%         d_j= vecnorm(dXYZ_Bj ,2, 2)  ; % not available for
        %         2017a
        d_j= dXYZ_Bj(:,1).^2+dXYZ_Bj(:,2).^2+dXYZ_Bj(:,3).^2 ;
        
        [~,bbj ] =sort(d_j) ; NIndGetJ=  round(length(d_j)*RatioBase) ;
        
        Ind_Bi=Ind_Bi(bbi(1:NIndGetI)) ;
        Ind_Bj=Ind_Bj(bbj(1:NIndGetJ)) ;
        
        % return
        %------------
        
        f523 = figure(523) ; %clf ; 
        f523.Name='Schematics and angle distribution' ; set(f523,'NumberTitle','off');

        subplot(2,1,1) ;cla;  hold on ;
        scatter3(MeanConf(:,1), MeanConf(:,2), MeanConf(:,3)  ,6, '.','CData',[0.5 0.5 0.5 ],'MarkerFaceAlpha',0.5 ) ;
        scatter3(MeanConf(Ind_Bi,1),MeanConf(Ind_Bi,2),MeanConf(Ind_Bi,3) ,'ro') ;
        scatter3(MeanConf(Ind_Bj,1),MeanConf(Ind_Bj,2),MeanConf(Ind_Bj,3) ,'bo') ;
        xlabel('X (nm)') ;ylabel('Y (nm)') ;zlabel('Z (nm)') ;
        title(' The average configuration and bases in consideration ')
        axis equal ;
        
        TrajOri =oxDNA_ex.Traj ;
        
        % StackTraj =reshape( permute(TrajOri(:,1:3,:)-0.4*TrajOri(:,4:6,:),[1, 2,3]) , size(TrajOri,1)*3 ,size(TrajOri,3) );
        % %-------------considering drifting
        % [coeff,score,latent,~,explained,PCA_mean] = pca(StackTraj') ;
        
        
        Bi_in_Mean = MeanConf(Ind_Bi ,:)  ;
        [coeff,score,latent,tsq1,explained,PCA_mean] = pca(Bi_in_Mean) ;
        Ref_BiAxis =  coeff(:,1)' ;
        % RefTwoPoints = [Bi_in_Mean(find(tsq1==max(tsq1)),:) ;Bi_in_Mean(find(tsq1==min(tsq1)),:) ];
        % scatter3(RefTwoPoints(:,1),RefTwoPoints(:,2),RefTwoPoints(:,3) , 86 ,'go')
        InnerProd_Bi = dot( Bi_in_Mean-ones(size(Bi_in_Mean,1),1)*PCA_mean , ones(size(Bi_in_Mean,1),1)*Ref_BiAxis ,2) ;
        RefTwoPoints_Bi = [Bi_in_Mean(find(InnerProd_Bi==max(InnerProd_Bi)),:) ;Bi_in_Mean(find(InnerProd_Bi==min(InnerProd_Bi)),:) ];
        scatter3(RefTwoPoints_Bi(:,1),RefTwoPoints_Bi(:,2),RefTwoPoints_Bi(:,3) , 86 ,'go','filled')
        RefInds_Bi = [find(InnerProd_Bi==max(InnerProd_Bi)) ,find(InnerProd_Bi==min(InnerProd_Bi)) ] ;
        RefInds_Bi=Ind_Bi(RefInds_Bi) ;
                
        Bj_in_Mean = MeanConf(Ind_Bj ,:)  ;
        [coeff,score,latent,tsq2,explained,PCA_mean] = pca(Bj_in_Mean) ;
        Ref_BjAxis =  coeff(:,1)' ;
        % RefTwoPoints = [Bj_in_Mean(find(tsq2==max(tsq2)),:) ;Bj_in_Mean(find(tsq2==min(tsq2)),:) ];
        % scatter3(RefTwoPoints(:,1),RefTwoPoints(:,2),RefTwoPoints(:,3) , 86 ,'go')
        InnerProd_Bj = dot( Bj_in_Mean-ones(size(Bj_in_Mean,1),1)*PCA_mean , ones(size(Bj_in_Mean,1),1)*Ref_BjAxis ,2) ;
        RefTwoPoints_Bj = [Bj_in_Mean(find(InnerProd_Bj==max(InnerProd_Bj)),:) ;Bj_in_Mean(find(InnerProd_Bj==min(InnerProd_Bj)),:) ];
        scatter3(RefTwoPoints_Bj(:,1),RefTwoPoints_Bj(:,2),RefTwoPoints_Bj(:,3) , 86 ,'go','filled')
        RefInds_Bj = [find(InnerProd_Bj==max(InnerProd_Bj)) ,find(InnerProd_Bj==min(InnerProd_Bj)) ] ;
        RefInds_Bj=Ind_Bj(RefInds_Bj) ;
        
        P1234 = MeanConf([RefInds_Bi, RefInds_Bj],: ) ;  %
        ddd1 = min ([norm(P1234(1,:)-P1234(3,:)) , norm(P1234(1,:)-P1234(4,:)) ] ) ;
        ddd2 = min ([norm(P1234(2,:)-P1234(3,:)) , norm(P1234(2,:)-P1234(4,:)) ] ) ;
        ddd3 = min ([norm(P1234(3,:)-P1234(1,:)) , norm(P1234(3,:)-P1234(2,:)) ] ) ;
        ddd4 = min ([norm(P1234(4,:)-P1234(1,:)) , norm(P1234(4,:)-P1234(2,:)) ] ) ;
        if ddd1<ddd2 
            RefInds_Bi =flip(RefInds_Bi) ;
        end
        if ddd3>ddd4 
            RefInds_Bj =flip(RefInds_Bj) ;
        end
        %---------  for the arrows
            Conf =MeanConf;
            Ref_BiAxis= diff( Conf(RefInds_Bi,1:3) ) ;
            Ref_BjAxis= diff( Conf(RefInds_Bj,1:3) ) ;
         Ref_BiAxis=Ref_BiAxis/norm(Ref_BiAxis) ;
         Ref_BjAxis=Ref_BjAxis/norm(Ref_BjAxis) ;
        
          Proj_i =dot(MeanConf(Ind_Bi,:) ,ones(length(Ind_Bi),1)*Ref_BiAxis ,2) ;
          Lmaxi = abs(min(Proj_i)- max(Proj_i))      ;      
          Oi = mean(MeanConf(Ind_Bi,:)) ;  
          Arrow_i = [Oi ;Oi+ 0.8*Lmaxi*Ref_BiAxis] ;
          quiver3(Arrow_i(1,1),Arrow_i(1,2),Arrow_i(1,3),diff(Arrow_i(:,1)),diff(Arrow_i(:,2)),diff(Arrow_i(:,3)),0 ,'LineWidth' ,3 ,'Color','r') ;
          
          Proj_j =dot(MeanConf(Ind_Bj,:) ,ones(length(Ind_Bj),1)*Ref_BjAxis ,2) ;
          Lmaxj = abs(min(Proj_j)- max(Proj_j))      ;      
          Oi = mean(MeanConf(Ind_Bj,:)) ;  
          Arrow_j = [Oi ;Oi+ 0.8*Lmaxj*Ref_BjAxis] ;
          quiver3(Arrow_j(1,1),Arrow_j(1,2),Arrow_j(1,3),diff(Arrow_j(:,1)),diff(Arrow_j(:,2)),diff(Arrow_j(:,3)),0 ,'LineWidth' ,3 ,'Color','b') ;
          
        %------------
        
        
%         return
        
        Save_AvgNvec = zeros( size(TrajOri,3) ,6) ;
        Angles = zeros( size(TrajOri,3) ,1) ;
        
%         AnglesOnRefPlane = zeros( size(TrajOri,3) ,1) ;
        for frame = 1 :  size(TrajOri,3)
            Conf = TrajOri(:,:,frame ) ;
            
            % Bi_bases = Conf(Ind_Bi ,1:3)  ;
            % [coeff,score,latent,~,explained,PCA_mean] = pca(Bi_bases) ;
            % Ref_BiAxis =  coeff(:,1)' ;
            Ref_BiAxis= diff( Conf(RefInds_Bi,1:3) ) ;
            
            % Bj_bases = Conf(Ind_Bj ,1:3)  ;
            % [coeff,score,latent,~,explained,PCA_mean] = pca(Bj_bases) ;
            % Ref_BjAxis =  coeff(:,1)' ;
            Ref_BjAxis= diff( Conf(RefInds_Bj,1:3) ) ;
            
            
            FlippingNvec_Bi = Conf(Ind_Bi , 7:9  ) ;
            InnerProd_Bi = dot( FlippingNvec_Bi , ones(size(FlippingNvec_Bi,1),1)*Ref_BiAxis ,2) ;
            FlipInds_i = InnerProd_Bi<0 ;
            AvgNvec_Bi = mean( [  FlippingNvec_Bi(~FlipInds_i,:) ; -FlippingNvec_Bi(FlipInds_i,:)] ) ;
            AvgNvec_Bi=AvgNvec_Bi/norm(AvgNvec_Bi) ;
            %
            %     QQ= [  FlippingNvec_Bi(~FlipInds_i,:) ; -FlippingNvec_Bi(FlipInds_i,:)] ;
            %      InnerProd_Test = dot( QQ , ones(size(QQ,1),1)*Ref_BiAxis ,2) ;
            %     FlipInds_test = InnerProd_Test<0 ; sum(FlipInds_test)
            %
            %
            
            FlippingNvec_Bj = Conf(Ind_Bj , 7:9  ) ;
            InnerProd_Bj = dot( FlippingNvec_Bj , ones(size(FlippingNvec_Bj,1),1)*Ref_BjAxis ,2) ;
            FlipInds_j = InnerProd_Bj<0 ;
            AvgNvec_Bj = mean( [  FlippingNvec_Bj(~FlipInds_j,:) ; -FlippingNvec_Bj(FlipInds_j,:)] ) ;
            AvgNvec_Bj=AvgNvec_Bj/norm(AvgNvec_Bj) ;
            
            
            Save_AvgNvec(frame,:) = [AvgNvec_Bi,AvgNvec_Bj ] ;
            Angles(frame) = acosd( dot(AvgNvec_Bi,AvgNvec_Bj)/norm(AvgNvec_Bi)/norm(AvgNvec_Bj)) ;
            
           %------------------
          
%             Conf ; 
%             [coeff,score,latent] = pca(Conf(:,1:3)) ;
%             NormalDir  =  coeff(:,3) ;
%             ProjAvg_Bi  = AvgNvec_Bi-  dot(AvgNvec_Bi, NormalDir)/norm(NormalDir) ;
%             ProjAvg_Bj  = AvgNvec_Bj-  dot(AvgNvec_Bj, NormalDir)/norm(NormalDir) ;
%             AnglesOnRefPlane(frame)= acosd( dot(ProjAvg_Bi,ProjAvg_Bj)/norm(ProjAvg_Bi)/norm(ProjAvg_Bj)) ;
            
%               ssdf=3
        end
%         save('Save_AvgNvec.mat','Save_AvgNvec');
        
        subplot(2,1,2) ; hold on ;
        ax=gca  ;
        if checkH.Value ==1
            cla ;
            ax.UserData.StrL = [];
        end
        
        if ~isfield(ax.UserData,'StrL' )
            ax.UserData.StrL= [] ;
        end
        
        
        pH= plot(Angles) ; 
%         pH2 = plot(AnglesOnRefPlane ,'-r') ;
        title(strcat('Selected bundles = [',num2str([Bi,Bj]),'] ' ,'{ }', 'R =', num2str(RatioBase),'{ }' ,'Mean =' , num2str(mean(Angles))  )  )   ;
        %         ax.UserData.StrL
        str = strcat( '[', num2str([Bi,Bj]) ,']' ,  ' ' , 'R =', num2str(RatioBase,'%3.1f') ) ;
        ax.UserData.StrL{end+1} = str;
        legend( ax.UserData.StrL ) ;
        
        % delete(d)
        ylim([0 180 ]) ;
        
        xlabel('frame'); ylabel('angle (degree)')
        
%         sdsf=3
        %----CM, save result into MATLAB .mat file 
%         oxDNA_ex.PathName
%             prompt = {'Enter Variable name:'};  dlg_title = 'Input';
%             num_lines = 1;   defaultans = {'Var1'};
%             file3_name = inputdlg(prompt,dlg_title,num_lines,defaultans); 
%             
%             %-------
%             prompt = {'Enter the value of grad'};
%             dlgtitle = 'Grad Value';
%             definput = {'0'};
%             dims = [1 40];
%             opts.Interpreter = 'tex';
%             answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
%             
%             SaveVar.str =  file3_name{1} ;
%             SaveVar.Angles =Angles ;
%             SaveVar.AnglesMean = mean(Angles) ;
%             SaveVar.AnglesStd = std(Angles) ;
%             SaveVar.grad = str2num( answer{1}) ;
%         
%             save(strcat('SQ2x4_r3_', file3_name{1},'.mat' ),'SaveVar')  ;
            
    end


end





