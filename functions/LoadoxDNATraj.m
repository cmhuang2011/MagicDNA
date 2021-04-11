function LoadoxDNATraj(src,evn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

rotate3d off ;
% ss_Assembly= findobj(0,'Tag','ss_Assembly') ;
% GetHyperB= ss_Assembly.UserData.HyperBundle ;
fH=gcf;
ax= findobj(fH,'Tag','oxDNATrajAxe1') ;
ax2= findobj(fH,'Tag','oxDNATrajAxe2')  ;

popOxdnaTraj= findobj(fH,'Tag','popOxdnaTraj') ;  popOxdnaTraj.Value=1;
listOxdnaTraj= findobj(fH,'Tag','listOxdnaTraj') ;  listOxdnaTraj.Value=1;

btn_oxDNATraj2= findobj(fH,'Tag','btn_oxDNATraj2') ;   % Compute deform
btn_oxDNATraj3= findobj(fH,'Tag','btn_oxDNATraj3') ;

btn_oxDNATraj4= findobj(fH,'Tag','btn_oxDNATraj4') ;   %use center
btn_oxDNATraj5= findobj(fH,'Tag','btn_oxDNATraj5') ;    %use orientaion

axes(ax);cltab; hold on ;axis equal;
axes(ax);


%----------------------------------
%             [obj.Topology_filename,PathName1,FilterIndex] = uigetfile({'*.dat;*.top','Topogyformat '},'Select the topology file');
%             obj.Spwd=pwd;     obj.PathName=PathName1;     cd(PathName1);
%             [Traj_filename,PathName2,FilterIndex]= uigetfile({'*.dat;*.conf','Trajectory file' },'Select the Trajectory file');
%             obj.Traj_filename=Traj_filename;

% SimulationFolder = uigetdir() ;
% oxDNA_ex = oxDNATrajObject({SimulationFolder});   % Only assign a folder with data, see code to change initial conf and traj file names .

oxDNA_ex = oxDNATrajObject({1,2,3});   % ask three inputs: topology, ini.Conf, traj .

% oxDNA_ex = oxDNATrajObject_polymer;   % control which class to use .
% oxDNA_ex = oxDNATrajObject_polymer_CM;   % control which class to use , Chao-Min use.





% ax.UserData.oxDNA_obj=oxDNA_ex ;
% oxDNA_ex=ax.UserData.oxDNA_obj ;

str2 = cellfun(@num2str,num2cell(0:oxDNA_ex.N_frame-1),'uniformoutput',0);
popOxdnaTraj.String=str2 ;
listOxdnaTraj.String=str2 ;
listOxdnaTraj.Max=oxDNA_ex.N_frame ;
src.UserData.oxDNA_ex=oxDNA_ex;
%------------


popOxdnaTraj.Callback =@(src,evn)popCallBack(src,evn,oxDNA_ex,ax);
listOxdnaTraj.Callback =@(src,evn)lstCallBack(src,evn,oxDNA_ex,ax);
btn_oxDNATraj2.Callback=@(src,evn)Mapping(src,evn,oxDNA_ex,ax,ax2);

btn_oxDNATraj5.Callback=@(src,evn)ColorfulRMSF;

btn_oxDNATraj4.Callback=@(src,evn)RMSD_RMSF(src,evn,oxDNA_ex,ax,ax2);  btn_oxDNATraj4.String='RMSD RMSF' ;


btn_oxDNATraj3.Callback=@(src,evn)ExportBILD(src,evn,oxDNA_ex,ax,ax2,popOxdnaTraj);

fprintf('finish loading trajectory \n ')

end

function  RMSD_RMSF(src,evn,oxDNA_ex,ax,ax2)
tic
TrajOri =oxDNA_ex.Traj ;

NoDriftTraj=TrajOri ;
% NoDriftTraj=zeros(size(TrajOri,1),9,size(TrajOri,3) ) ;

% StackTraj =reshape(permute(TrajOri(:,1:3,:),[1, 3,2]) ,size(TrajOri,1)*size(TrajOri,3),3 );
% MeanConf=oxDNA_ex.MeanConf ;
% MeanConf = oxDNA_ex.Traj(:,1:3,101) ; % temp
% figure; hold on; scatter3(MeanConf(:,1),MeanConf(:,2),MeanConf(:,3),'.')
% 1080
ExcInd = [1:981, 7462:7560,  size(oxDNA_ex.Strand,1)-1079:size(oxDNA_ex.Strand,1)];
IndAccount=ones(size(oxDNA_ex.Strand,1),1)==1 ;
IndAccount(ExcInd) =false ;

B3byNaLL=  transpose(TrajOri(:,1:3,1 ) );   % target

for frI= 1:size(TrajOri,3)
    A3byNaLL=  transpose(TrajOri(:,1:3,frI )) ;  % which will be transformed.
    
%     [regParams,~,~]=absor(A3byNaLL(:,IndAccount),B3byNaLL(:,IndAccount));
    [regParams,~,~]=absor(A3byNaLL,B3byNaLL);
    
    A_prime = transpose(regParams.R*A3byNaLL + regParams.t  ); % Base centers
    
    A_BVec=  transpose(TrajOri(:,4:6,frI )) ;  % orientations does not consider translation
    BVec_prime = transpose(regParams.R*A_BVec );
    A_NVec=  transpose(TrajOri(:,7:9,frI )) ;
    NVec_prime = transpose(regParams.R*A_NVec );
    %     A_prime=transpose(A_prime) ;
    %  scatter3(A_prime(:,1),A_prime(:,2),A_prime(:,3),'.') ;
    NoDriftTraj(:,:,frI) =[A_prime,BVec_prime,NVec_prime] ;
end
toc

oxDNA_ex.Traj=NoDriftTraj ;

% oxDNA 2 model
POS_MM_BACK1 = -0.3400 ;
POS_MM_BACK2 = 0.3408 ;
a2 = cross(NoDriftTraj(:,7:9,:),NoDriftTraj(:,4:6,:) ) ;

% Coeff=0.6 ;
% BackBoneT = NoDriftTraj(:,1:3,:)  + Coeff*NoDriftTraj(:,4:6,:)  ; % hard code for HB center

BackBoneT = NoDriftTraj(:,1:3,:)  + POS_MM_BACK1*NoDriftTraj(:,4:6,:) + POS_MM_BACK2*a2 ; % backbone



StackTraj =reshape(permute(BackBoneT,[1, 2,3]) ,size(NoDriftTraj,1)*3,size(NoDriftTraj,3) );
%-------------considering drifting
[coeff,score,latent,~,explained,PCA_mean] = pca(StackTraj') ;
MeanConf = reshape(PCA_mean, size(TrajOri,1) ,3);
oxDNA_ex.MeanConf=MeanConf;



%-------------RMSD
Ref_RMSD = NoDriftTraj(:,1:3,1 ); %initial configuration as reference
% Ref_RMSD = MeanConf ; % deviation from mean, for Receptor
RMSD_By_Frame = zeros(size(TrajOri,3),1 ) ;
% RMSD_By_xComp = zeros(size(TrajOri,3),1 ) ;
% RMSD_By_yComp = zeros(size(TrajOri,3),1 ) ;
% RMSD_By_zComp = zeros(size(TrajOri,3),1 ) ;
%
%
for k=1:length(RMSD_By_Frame)
    dxyz = NoDriftTraj(:,1:3,k ) - Ref_RMSD ;
    %     dxyz=dxyz(ReceptorPara.RecTraceBase(2),:) ;  % 
    %     %  dxyz=dxyz(:,1) ;  % 
    RMSD_By_Frame(k) = rms(dxyz(:)) ;
    %     RMSD_By_xComp(k) = rms(dxyz(:,1)) ;
    %     RMSD_By_yComp(k) = rms(dxyz(:,2)) ;
    %     RMSD_By_zComp(k) = rms(dxyz(:,3)) ;
    
end



figure(243) ; clf ; hold on ;
subplot(4 ,3 ,[1:6   ] ) ;
plot(0:length(RMSD_By_Frame)-1, RMSD_By_Frame,'LineWidth',2 ); title(strcat( 'RMSD. Avg=', num2str(mean(RMSD_By_Frame))  )) ;
xlabel('frame') ; ylabel('RMSD (nm)') ;
% toc
%------------------
Ref_RMSF = MeanConf ;
RMSF_By_Base = zeros(size(TrajOri,1),1 ) ;
for k2= 1 :length(RMSF_By_Base)
    r_ofTime = permute( NoDriftTraj(k2,1:3,:) , [3,2,1 ]) ;
    dxyz2 = r_ofTime-  repmat( Ref_RMSF(k2,:), size(r_ofTime,1) ,1) ;
    RMSF_By_Base(k2) =rms(dxyz2(:)) ;
    
end
subplot(4 ,3 ,[7 8 10 11   ] ); ax_RMSF =gca ;
plot(RMSF_By_Base);str={strcat( 'Avg = ', num2str(mean(RMSF_By_Base))  ) };

title(str) ;
xlabel('Base') ; ylabel('RMSF (nm)') ;
%++++
subplot(4 ,3 ,[ 9 12] ) ;
h2= histogram(RMSF_By_Base(1:end),20 ) ;
view([90 -90]);
xlabel('RMSF (nm)') ;  ylabel('N') ;
% return
% toc
%------------------
axes(ax2);cla; hold on ;axis equal;axes(ax2);
% figure(234);clf; hold on ;axis equal;


CClorData=RMSF_By_Base ;

%---------- custom color for fig 5
% % colormap jet ;
% CC=colormap ;
% BM=oxDNA_ex.BM  ;
% AllBundle=unique(BM) ;
% %  WireFrameBundles = [1:3]' ;
% LatticeBundles= [1:10]' ;
% % SurfaceBundlesl = [1:10] ;
% % WireFrameBundles = setdiff(AllBundle, union(SurfaceBundlesl,LatticeBundles)) ;
% 
%   WireFrameBundles=[21] ;
%   SurfaceBundlesl = setdiff(AllBundle, union(WireFrameBundles,LatticeBundles)) ;
% 
% Inds= ismember(BM ,LatticeBundles ) ;
% if sum(Inds)>0
% sH=scatter3( Ref_RMSF(Inds,1) , Ref_RMSF(Inds,2) ,Ref_RMSF(Inds,3)  ) ;
% sH.CData =CC(8,:);
% MaxV = max(RMSF_By_Base(Inds) ) ;MinV = min(RMSF_By_Base(Inds) ) ;
% vq = interp1([MinV;MaxV],[0 0.5 1 0 ;0.2 1 0.6 0.5],RMSF_By_Base(Inds),'linear','extrap') ;
% % vq = interp1([MinV;MaxV-0.2],[0 0 1 0 ;0 1 1 0.5],RMSF_By_Base(Inds),'linear','extrap') ;
% % vq = interp1([MinV;MaxV ],[ 0.9 0.9 1  ;0 0 1],RMSF_By_Base(Inds)) ;
% vq(vq>=1)=1; vq(vq<=0)=0;
% sH.CData =vq(:,1:3) ;
% sH.UserData=vq(:,4) ;
% hold on;
% end
% 
% Inds= ismember(BM,WireFrameBundles  ) ;
% if sum(Inds)>0
% sH2=scatter3( Ref_RMSF(Inds,1) , Ref_RMSF(Inds,2) ,Ref_RMSF(Inds,3)  ) ;    sH2.CData =CC(56,:);   %CC(56,:)
% MaxV = max(RMSF_By_Base(Inds) ) ; MinV = min(RMSF_By_Base(Inds) ) ;
% vq = interp1([MinV;MaxV ], [1 0.5 0  0 ; 0.8 0.4 0 0.5],RMSF_By_Base(Inds),'linear','extrap' ) ;
% % vq = interp1([MinV;MaxV ], [1 0.8 0.8 0 ; 1 0 1 0.5],RMSF_By_Base(Inds),'linear','extrap' ) ;
% % vq = interp1([MinV;MaxV ],[1 0.9 0.9  ; 1 0 0],RMSF_By_Base(Inds)) ;
% % [1 0.8 0.8 0 ; 1 0 1 0.5]
% vq(vq>=1)=1; vq(vq<=0)=0;
% sH2.CData =vq(:,1:3) ;
% sH2.UserData=vq(:,4) ;
% end
% 
% 
% Inds= ismember(BM,SurfaceBundlesl  ) ;
% if sum(Inds)>0
% sH3=scatter3( Ref_RMSF(Inds,1) , Ref_RMSF(Inds,2) ,Ref_RMSF(Inds,3)  ) ;    sH3.CData =[0,0.7,0  ];   %
% MaxV = max(RMSF_By_Base(Inds) ) ;MinV = min(RMSF_By_Base(Inds) ) ;
% vq = interp1([MinV;MaxV-0.5],[0 0.6 0.6 0.5; 1 1 0 0],RMSF_By_Base(Inds) ,'linear','extrap' ) ;
% % vq = interp1([MinV;MaxV],[0 0.5 0.4 0.5;1 1 0 0],RMSF_By_Base(Inds) ,'linear','extrap' ) ;
% vq(vq>=1)=1; vq(vq<=0)=0;
% % vq = interp1([MinV;MaxV ],[0.7 1 0.7 0; 1 1 0 0.5],RMSF_By_Base(Inds)) ;
% % [0.8 1 0.6 0; 0 0.6 0 0.5]
% % vq = interp1([MinV;MaxV ],[0.9 1 0.9  ; 0,1,0 ],RMSF_By_Base(Inds)) ;
% sH3.CData =vq(:,1:3) ;
% sH3.UserData=vq(:,4) ;
% end
% % sH=scatter3( Ref_RMSF(:,1) , Ref_RMSF(:,2) ,Ref_RMSF(:,3)  ) ;
% % sH.CData =CC(8,:);   % lattice use 8 :
% 
% % sH.CData =CC(56,:);   % lattice use 8 :
% %  figure; hist(RMSF_By_Base(Inds))
% 
% return
%---------
% CClorData=ones(size(RMSF_By_Base)) ;

[f,g,h]=plot3k({Ref_RMSF(:,1) , Ref_RMSF(:,2) ,Ref_RMSF(:,3)}, 'ColorData',CClorData ,'ColorRange',[min(RMSF_By_Base), max(RMSF_By_Base)] ) ;
str={'RMSF' ; strcat( 'Avg = ', num2str(mean(RMSF_By_Base))  ) };
ax=gca; ax.UserData.Ref_RMSF=Ref_RMSF; 
title(str) ;

%-----for shaded output
QQ =reshape(  permute(NoDriftTraj(:,1:3,:) , [1 3 2]) , size(NoDriftTraj,1)*size(NoDriftTraj,3) , 3) ;
% cc= 0 ; 
% for k = 1: size(NoDriftTraj ,3)
%     boundNode = boundary(NoDriftTraj(:,1,k) , NoDriftTraj(:,2,k) ,NoDriftTraj(:,3,k) ) ;
%     cc=cc+ length( unique(boundNode)) ;
% end
% QQ = zeros(cc ,3) ;cc2=0;
% for k = 1: size(NoDriftTraj ,3)
%     boundNode = boundary(NoDriftTraj(:,1,k) , NoDriftTraj(:,2,k) ,NoDriftTraj(:,3,k) ) ;
%     n = length( unique(boundNode)) ;
%     QQ(cc2+1:cc2+n ,:) = NoDriftTraj(unique(boundNode),1:3 ,k) ;
%     cc2=cc2+ length( unique(boundNode)) ;
% end
ax.UserData.AlltrajXYZ= QQ;
%---------------


end


function  ExportBILD(src,evn,oxDNA_ex,ax,ax2,popOxdnaTraj)
prompt = {'Enter BILD File name:'};    dlg_title = 'Input';  num_lines = 1;
defaultans = {'oxDNA_unrelaxed'};    answer = inputdlg(prompt,dlg_title,num_lines,defaultans) ;
opts.Interpreter = 'tex'; opts.Default = 'coarsed-grained';

answer_CG = questdlg({'\fontsize{15} Export CG or ribbong model ?' }  , ...
    'chimera options', 'coarsed-grained','ribbon',opts);
switch answer_CG
    case 'coarsed-grained'
        file2_name=strcat(answer{1},'_CG_', '.bild')     ;
        Ribbon=0;
    case 'ribbon'
        file2_name=strcat(answer{1},'_ribbon_' , '.bild')      ;
        Ribbon=1;
end
frame =popOxdnaTraj.Value;
oxDNA_ex.ExportBildRibbonTraj(Ribbon,frame,answer{1}) ;
end

function MappingUseCenter(src,evn,oxDNA_ex,ax,ax2)
popOxdnaTraj= findobj(0,'Tag','popOxdnaTraj') ;
SelectFrame = popOxdnaTraj.Value  ;

ProvaConf = oxDNA_ex.Traj(:,:,1) ;
TargeConf = oxDNA_ex.Traj(:,:,SelectFrame) ;
BM =oxDNA_ex.BM ;
axes(ax2);cla; hold on ;axis equal;axes(ax2);

A3byNaLL=  transpose(ProvaConf(:,1:3)) ;
B3byNaLL=  transpose(TargeConf(:,1:3) );

[regParamsaLL,~,~]=absor(B3byNaLL,A3byNaLL);
b_prime = transpose(regParamsaLL.R*B3byNaLL + regParamsaLL.t  );

b_NvecPrime = transpose(regParamsaLL.R*transpose(TargeConf(:,7:9))   );
BendingAngle = zeros(size(b_NvecPrime,1) ,1  ) ;
for sti =1 :max(oxDNA_ex.Strand)
    Inds= find(oxDNA_ex.Strand ==sti);
    
    VecA1 =b_NvecPrime( Inds(1)  ,:) ;
    VecA2 =b_NvecPrime( Inds(1+1)  ,:) ;
    BendingAngle(Inds(1) )=acosd(dot(VecA1,VecA2)) ;
    
    for basej = 2: length( Inds)-1
        Vecjm1 =b_NvecPrime( Inds(basej-1)  ,:) ;
        Vecj1 =b_NvecPrime( Inds(basej)  ,:) ;
        Vecjp1 =b_NvecPrime( Inds(basej)+1  ,:) ;
        
        BendingAngle(Inds(basej)  )=mean([ acosd(dot(Vecjm1,Vecj1)) ,       acosd(dot(Vecjp1,Vecj1)) ]  )   ;
    end
    Vecnm1 =b_NvecPrime( Inds(end-1)  ,:) ;
    Vecn =b_NvecPrime( Inds(end)  ,:) ;
    BendingAngle(Inds(end) )=acosd(dot(Vecnm1,Vecn)) ;
    
    %        sdfsff=3
end

plot3k({b_prime(:,1) , b_prime(:,2) ,b_prime(:,3)}, 'ColorData',BendingAngle ,'ColorRange',[min(BendingAngle), max(BendingAngle)] ) ;
str={strcat('Use NVec mean d = ',num2str(mean(BendingAngle))) ;  strcat('min= ',num2str(min(BendingAngle),3)  ,' max= ',num2str(max(BendingAngle),3 )) } ;
title( str  );

oxDNA_ex.BaseColor =BendingAngle ;
frpintf('End of MappingUSingCenter \n')
end

function MappingUseNVec(src,evn,oxDNA_ex,ax,ax2)
popOxdnaTraj= findobj(0,'Tag','popOxdnaTraj') ;
SelectFrame = popOxdnaTraj.Value  ;

ProvaConf = oxDNA_ex.Traj(:,:,1) ;
TargeConf = oxDNA_ex.Traj(:,:,SelectFrame) ;
BM =oxDNA_ex.BM ;
axes(ax2);cla; hold on ;axis equal;axes(ax2);

A3byNaLL=  transpose(ProvaConf(:,1:3)) ;
B3byNaLL=  transpose(TargeConf(:,1:3) );

[regParamsaLL,~,~]=absor(B3byNaLL,A3byNaLL);
b_prime = transpose(regParamsaLL.R*B3byNaLL + regParamsaLL.t  );

b_NvecPrime = transpose(regParamsaLL.R*transpose(TargeConf(:,7:9))   );
BendingAngle = zeros(size(b_NvecPrime,1) ,1  ) ;
for sti =1 :max(oxDNA_ex.Strand)
    Inds= find(oxDNA_ex.Strand ==sti);
    
    VecA1 =b_NvecPrime( Inds(1)  ,:) ;
    VecA2 =b_NvecPrime( Inds(1+1)  ,:) ;
    BendingAngle(Inds(1) )=acosd(dot(VecA1,VecA2)) ;
    
    for basej = 2: length( Inds)-1
        Vecjm1 =b_NvecPrime( Inds(basej-1)  ,:) ;
        Vecj1 =b_NvecPrime( Inds(basej)  ,:) ;
        Vecjp1 =b_NvecPrime( Inds(basej)+1  ,:) ;
        
        BendingAngle(Inds(basej)  )=mean([ acosd(dot(Vecjm1,Vecj1)) ,       acosd(dot(Vecjp1,Vecj1)) ]  )   ;
    end
    Vecnm1 =b_NvecPrime( Inds(end-1)  ,:) ;
    Vecn =b_NvecPrime( Inds(end)  ,:) ;
    BendingAngle(Inds(end) )=acosd(dot(Vecnm1,Vecn)) ;
    
    %        sdfsff=3
end

plot3k({b_prime(:,1) , b_prime(:,2) ,b_prime(:,3)}, 'ColorData',BendingAngle ,'ColorRange',[min(BendingAngle), max(BendingAngle)] ) ;
str={strcat('Use NVec mean d = ',num2str(mean(BendingAngle))) ;  strcat('min= ',num2str(min(BendingAngle),3)  ,' max= ',num2str(max(BendingAngle),3 )) } ;
title( str  );

oxDNA_ex.BaseColor =BendingAngle ;
end




function Mapping(src,evn,oxDNA_ex,ax,ax2)
popOxdnaTraj= findobj(0,'Tag','popOxdnaTraj') ;
SelectFrame = popOxdnaTraj.Value  ;

ProvaConf = oxDNA_ex.Traj(:,:,1) ;
TargeConf = oxDNA_ex.Traj(:,:,SelectFrame) ;
BM =oxDNA_ex.BM ;

axes(ax2);cla; hold on ;axis equal;axes(ax2);
dSave=zeros(size(ProvaConf,1),1) ;
for buni =1:max(BM)
    Inds = BM==buni ;
    
    A3byN=  transpose(ProvaConf(Inds,1:3)) ;
    B3byN=  transpose(TargeConf(Inds,1:3) );
    
    [regParams,~,~]=absor(A3byN,B3byN);
    
    A_prime = transpose(regParams.R*A3byN + regParams.t  );
    
    dxyz=A_prime-TargeConf(Inds,1:3) ;
    d= sqrt( dxyz(:,1).^2 + dxyz(:,2).^2 +dxyz(:,3).^2 ) ;
    dSave(Inds) =d ;
    
end
A3byNaLL=  transpose(ProvaConf(:,1:3)) ;
B3byNaLL=  transpose(TargeConf(:,1:3) );

[regParamsaLL,~,~]=absor(B3byNaLL,A3byNaLL);
b_prime = transpose(regParamsaLL.R*B3byNaLL + regParamsaLL.t  );
plot3k({b_prime(:,1) , b_prime(:,2) ,b_prime(:,3)}, 'ColorData',dSave ,'ColorRange',[min(dSave), max(dSave)] ) ;
%      plot3k({TargeConf(:,1) , TargeConf(:,2) ,TargeConf(:,3)}, 'ColorData',dSave ,'ColorRange',[min(dSave), max(dSave)] ) ;
str={strcat('mean d = ',num2str(mean(dSave))) ;  strcat('min= ',num2str(min(dSave),3)  ,' max= ',num2str(max(dSave),3 )) } ;
title( str  );
fprintf('End of Compute Deform, Bundle By bundle\n')

end

function popCallBack(src,evn,oxDNA_ex,ax)
axes(ax);cla; hold on ;axis equal; axes(ax);   %after cltab, remember to set the axe
oxDNA_ex.visualTraj(src.Value,1) ;
title(strcat('Frame =' , num2str(src.Value-1)) ) ;
end

function lstCallBack(src,evn,oxDNA_ex,ax)
axes(ax);cla; hold on ;axis equal; axes(ax);
frames =  src.Value ;
oxDNA_ex.visualTraj(frames,1) ;
title(strcat('Frame =' , num2str(frames-1)) ) ;
end


