classdef oxDNATrajObject < oxDNA_model
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        Spwd=[];
        PathName=[];
        Topology_filename=[];     % input file 1
        provaConf=[];                   % input file2
        Traj_filename=[];           % input file 3
        
        BM=[];
        %         objname='abc' ;                  % {cut from json name, can change after creating object}
        Boxs=[];
        
        Traj=[];
        N_frame=[];
        Strand=[];
        TSeq;
        MeanConf=[];
        %         colorRGB=[];
        %         BaseColor =[];
        %         N=[] ;
        Primes35 =[];
        MutualTrapPair =[]; % from dsremain
        
        
    end
    properties (Hidden = true)
        objname='abc' ;                  % {cut from json name, can change after creating object}
        PcaTraj= [];                    %if not, try run oxdna_SLL_cos_paper
        colorRGB=[];
        BaseColor =[];
        JsonFile=[];
        
    end
    properties (Constant,Hidden = true)
        %         CylR=0.06;
        %       SphereR=0.2 ;
        CylR=0.1;   %0.1
        SphereR=0.3 ;
        
        %         POS_BACK = -0.4 ;
        %         POS_MM_BACK1 = -0.34 ;
        %         POS_MM_BACK2 = 0.3408 ;
        %         POS_STACK = 0.34 ;
        %         POS_BASE = 0.4 ;
        %         GAMMA = 0.74 ; %POS_STACK - POS_BACK
        %         %-----------
        %         FENE_EPS =2.0;
        %         FENE_R0_OXDNA =0.7525;
        %         FENE_R0_OXDNA2 =0.7564;
        %         FENE_DELTA =0.25;
        %         FENE_DELTA2 =0.0625;
        %         %----------- excluded volume
        %         %1-> back-back 2->base-base 3->base-back 4->back-base
        %         EXCL_EPS =2.0;
        %         EXCL_S1 =0.7;  %sigma
        %         EXCL_S2 =0.33;
        %         EXCL_S3 =0.515;
        %         EXCL_S4 =0.515;
        %         EXCL_R1 =0.675;  %r star
        %         EXCL_R2 =0.32;
        %         EXCL_R3 =0.5;
        %         EXCL_R4 =0.5;
        %         EXCL_B1 =892.016223343;
        %         EXCL_B2 =4119.70450017;
        %         EXCL_B3 =1707.30627298;
        %         EXCL_B4 =1707.30627298;
        %         EXCL_RC1 =0.711879214356;
        %         EXCL_RC2 =0.335388426126;
        %         EXCL_RC3 =0.52329943261;
        %         EXCL_RC4 =0.52329943261;
        %-----------------
        
        
        
    end
    
    
    methods
        
        function Example(obj)
            return
            wbb=oxDNAobject;
            %             wbb.getInterpolate
            wbb.ExportBildRibbonTraj
            wbb.printCompChimera
            
            oxDNA=helpCADDOM(5) ;
            load('ReceptorPara.mat','ReceptorPara') ;
            ReceptorPara ;
            ReceptorPara.RecTraceBase ;
            oxDNA.Trace(ReceptorPara.RecTraceBase) ;
            oxDNA.Trace(ReceptorPara.Inds) ;
            oxDNA.Trace(ReceptorPara.RecTraceBase) ;
            figure ;hold on; Frm=101 ; oxDNA.visualTraj(Frm,1) ;
            scatter3(oxDNA.Traj(ReceptorPara.RecTraceBase,1 ,Frm) ,oxDNA.Traj(ReceptorPara.RecTraceBase,2 ,Frm),oxDNA.Traj(ReceptorPara.RecTraceBase,3 ,Frm),'filled') ;
            text(oxDNA.Traj(ReceptorPara.RecTraceBase,1 ,Frm) ,oxDNA.Traj(ReceptorPara.RecTraceBase,2 ,Frm),oxDNA.Traj(ReceptorPara.RecTraceBase,3 ,Frm),num2str([1:12]'))
        end
        function BaseSite = GetBaseSite(obj,T )
            BaseSite= T(:,1:3) + obj.POS_BASE*T(:,4:6) ;
        end
        function StackSite= GetStackSite(obj,T )
            StackSite= T(:,1:3) + obj.POS_STACK*T(:,4:6) ;
        end
        function BackboneSite= GetBackboneSite(obj,T )
            a2 = cross(T(:,7:9),T(:,4:6) ) ;
            BackboneSite= T(:,1:3) + obj.POS_MM_BACK1*T(:,4:6) + obj.POS_MM_BACK2*a2 ;
        end
        function Energies = Cal_FENE(obj )
            TotalFrame = obj.N_frame ; TotalBase = length(obj.Strand) ;
            Energies= zeros(TotalFrame,1) ;
            for k =1:TotalFrame
                BackboneSite= GetBackboneSite(obj,obj.Traj(:,:,k) ) ;
                for j = 1:max(obj.Strand)
                    BaseInd = obj.Strand==j ;
                    SubBB = BackboneSite(BaseInd,:) ;
                    
                    dxyz= diff(SubBB ,1) ;
                    dr = sqrt(dxyz(:,1).^2+dxyz(:,2).^2+dxyz(:,3).^2 ) ;
                    Vs_FENE = -obj.FENE_EPS*0.5*log(1- (dr-obj.FENE_R0_OXDNA2).^2/obj.FENE_DELTA2) ;
                    
                    Energies(k)= Energies(k) + sum(Vs_FENE) ;
                end
                fprintf('FENE Energy at frame %i = %7.6f \n',k, Energies(k)/TotalBase) ;
            end
            Energies=Energies/TotalBase;
        end
        
        
        function Energies = Cal_exc_bonded(obj)
            TotalFrame = obj.N_frame ;             TotalBase = length(obj.Strand) ;
            Energies= zeros(TotalFrame,1) ;
            for k =1:TotalFrame
                BackboneSite= GetBackboneSite(obj,obj.Traj(:,:,k) ) ;
                BaseSite=GetBaseSite(obj,obj.Traj(:,:,k) ) ;
                for j = 1:max(obj.Strand)
                    BaseInd = obj.Strand==j ;
                    SubBackbone = BackboneSite(BaseInd,:) ;
                    SubBase =  BaseSite(BaseInd,:) ;
                    %----2.Base-Base
                    dxyz= diff(SubBase ,1) ;   dr = sqrt(dxyz(:,1).^2+dxyz(:,2).^2+dxyz(:,3).^2 ) ;
                    LargerThanRC = dr> obj.EXCL_RC2 ;
                    SmallerThanRstar = dr <  obj.EXCL_R2 ;
                    Middle= and(LargerThanRC==0 ,SmallerThanRstar==0) ;
                    
                    Eng1=  get_Vlj(obj , dr(SmallerThanRstar ), obj.EXCL_EPS , obj.EXCL_S2) ;
                    Eng2=  get_Vsmooth(obj, dr(Middle ), obj.EXCL_B2, obj.EXCL_RC2, obj.EXCL_EPS ) ;
                    Energies(k)= Energies(k) + Eng1 +Eng2  ;
                    %                   %----3.base-back
                    dxyz=SubBase(1:end-1 ,:) - SubBackbone(2:end,:) ; dr = sqrt(dxyz(:,1).^2+dxyz(:,2).^2+dxyz(:,3).^2 ) ;
                    LargerThanRC = dr> obj.EXCL_RC3 ;
                    SmallerThanRstar = dr <  obj.EXCL_R3 ;
                    Middle= and(LargerThanRC==0 ,SmallerThanRstar==0) ;
                    
                    Eng1=  get_Vlj(obj , dr(SmallerThanRstar ), obj.EXCL_EPS , obj.EXCL_S3) ;
                    Eng2=  get_Vsmooth(obj, dr(Middle ), obj.EXCL_B3, obj.EXCL_RC3, obj.EXCL_EPS ) ;
                    Energies(k)= Energies(k) + Eng1 +Eng2  ;
                    %----4.back-base
                    dxyz=SubBackbone(1:end-1 ,:) - SubBase(2:end,:) ; dr = sqrt(dxyz(:,1).^2+dxyz(:,2).^2+dxyz(:,3).^2 ) ;
                    LargerThanRC = dr> obj.EXCL_RC4 ;
                    SmallerThanRstar = dr <  obj.EXCL_R4 ;
                    Middle= and(LargerThanRC==0 ,SmallerThanRstar==0) ;
                    
                    Eng1=  get_Vlj(obj , dr(SmallerThanRstar ), obj.EXCL_EPS , obj.EXCL_S4) ;
                    Eng2=  get_Vsmooth(obj, dr(Middle ), obj.EXCL_B4, obj.EXCL_RC4, obj.EXCL_EPS ) ;
                    Energies(k)= Energies(k) + Eng1 +Eng2  ;
                    %-----
                end
                fprintf('Exc.V Energy at frame %i = %7.6f \n',k, Energies(k)/TotalBase) ;
            end
            Energies=Energies/TotalBase;
        end
        
        
        
        
        function Trace(obj, WhichBase)
            figure;clf; hold on;
            subplot(3,3,1) ;
            CC=get(gca,'colororder') ;
            SaveG1=[];   SaveG2=[];   SaveG3=[]; Ns= 5;
            Colors=[ repmat(CC(1,:), 4,1) ; repmat(CC(2,:), 4,1) ; repmat(CC(3,:), 4,1)] ;
            for k=1:length(WhichBase)
                XArr= reshape(obj.Traj(WhichBase(k),1,:) ,size(obj.Traj,3) ,1) ;  XArr=XArr-104.454 ;
                if k<5
                    subplot(3,3,1) ;   hold on;
                    SaveG1=[SaveG1;XArr];
                    pHx= plot(XArr(1:Ns:end),'Color',Colors(k,:),'LineWidth',1 ); title('X')
                elseif k>8
                    subplot(3,3,2) ;      hold on;
                    SaveG2=[SaveG2;XArr];
                    pHx= plot(XArr(1:Ns:end),'Color',Colors(k,:),'LineWidth',1 ); title('X')
                else
                    subplot(3,3,3) ;     hold on;
                    SaveG3=[SaveG3;XArr];
                    pHx= plot(XArr(1:Ns:end),'Color',Colors(k,:),'LineWidth',1 ); title('X')
                end
                
                
                YArr= reshape(obj.Traj(WhichBase(k),2,:) ,size(obj.Traj,3) ,1) ;  YArr=YArr-YArr(1) ;
                ZArr= reshape(obj.Traj(WhichBase(k),3,:) ,size(obj.Traj,3) ,1) ;  ZArr=ZArr-ZArr(1) ;
                %                 subplot(3,1,1); hold on;
                subplot(3,3,4); hold on;   pHy= plot(YArr,'Color',Colors(k,:),'LineWidth',0.1 ); title('Y')
                subplot(3,3,7); hold on;   pHz= plot(ZArr,'Color',Colors(k,:),'LineWidth',0.1 ); title('Z')
            end
            subplot(3,3,1) ;  title( strcat('RMS = ', num2str(rms(SaveG1)) ) ) ;
            subplot(3,3,2) ;  title( strcat('RMS = ', num2str(rms(SaveG2)) ) ) ;
            subplot(3,3,3) ;  title( strcat('RMS = ', num2str(rms(SaveG3)) ) ) ;
        end
        function Trace2(obj, WhichBase)
            figure;clf; hold on;
            subplot(3,3,1) ;
            CC=get(gca,'colororder') ;
            SaveG1=[];   SaveG2=[];   SaveG3=[]; Ns= 1;
            Colors=[ repmat(CC(1,:), 4,1) ; repmat(CC(2,:), 4,1) ; repmat(CC(3,:), 4,1)] ;
            for k=1:length(WhichBase)
                XArr= reshape(obj.Traj(WhichBase(k),1,:) ,size(obj.Traj,3) ,1) ;  % XArr=XArr-104.454 ;
                if k<5
                    subplot(3,3,1) ;   hold on;
                    SaveG1=[SaveG1;XArr];
                    pHx= plot(movmean(XArr, Ns),'Color',Colors(k,:),'LineWidth',1 ); title('X')
                elseif k>8
                    subplot(3,3,2) ;      hold on;
                    SaveG2=[SaveG2;XArr];
                    pHx= plot(movmean(XArr, Ns),'Color',Colors(k,:),'LineWidth',1 ); title('X')
                else
                    subplot(3,3,3) ;     hold on;
                    SaveG3=[SaveG3;XArr];
                    pHx= plot(movmean(XArr, Ns),'Color',Colors(k,:),'LineWidth',1 ); title('X')
                end
                
                
                YArr= reshape(obj.Traj(WhichBase(k),2,:) ,size(obj.Traj,3) ,1) ;  YArr=YArr-YArr(1) ;
                ZArr= reshape(obj.Traj(WhichBase(k),3,:) ,size(obj.Traj,3) ,1) ;  ZArr=ZArr-ZArr(1) ;
                %                 subplot(3,1,1); hold on;
                subplot(3,3,4); hold on;   pHy= plot(YArr,'Color',Colors(k,:),'LineWidth',0.1 ); title('Y')
                subplot(3,3,7); hold on;   pHz= plot(ZArr,'Color',Colors(k,:),'LineWidth',0.1 ); title('Z')
            end
            subplot(3,3,1) ;  title( strcat('RMS = ', num2str(rms(SaveG1)) ) ) ;
            subplot(3,3,2) ;  title( strcat('RMS = ', num2str(rms(SaveG2)) ) ) ;
            subplot(3,3,3) ;  title( strcat('RMS = ', num2str(rms(SaveG3)) ) ) ;
        end
        
        function  obj=oxDNATrajObject(varargin)
%             nargin 
%             varargin
%             length(varargin{1})
%             return 
            %             addpath('jsonlab') ;
            if length(varargin{1})==1
                
                %-----use for batch input
                PathName1=[varargin{1}{1} filesep]; PathName2=[varargin{1}{1} filesep];
                obj.Topology_filename = 'prova.top' ;
                obj.Spwd=pwd;  obj.PathName=PathName1;   cd(PathName1);
%                 provaConf = 'stage4.conf' ; 
                provaConf = 'prova33.conf' ; % M
%                 provaConf='relaxed_go3.conf';
%                  provaConf='lastConf_simulation.conf';
%                  provaConf = 'ForClamp.conf' ;

                obj.provaConf=provaConf;
%                  Traj_filename ='traj1e8_f100.dat' ;  % 
%                 Traj_filename ='traj1e7.dat' ;  % 
%                 Traj_filename ='traj2e8_f200_dt2.dat' ;  % 
                 Traj_filename ='traj_Clamp1e8_f200_dt1.dat' ;  % 
% Traj_filename='traj_rpp.dat';
%                                 Traj_filename ='traj1e8_f100.dat' ;

                obj.Traj_filename=Traj_filename;
                %                  obj.Topology_filenam=
                %                  sddsf=3
            else
                [obj.Topology_filename,PathName1,~] = uigetfile({'*.dat;*.top','Topogyformat '},'Select the topology file');
                obj.Spwd=pwd;     obj.PathName=PathName1;     cd(PathName1);
                
                [provaConf,PathName2,~]= uigetfile({'*.dat;*.conf','Initial Conf file' },'Select the initial Configuration file');
                obj.provaConf=provaConf;
                
                [Traj_filename,PathName2,~]= uigetfile({'*.dat;*.conf','Trajectory file' },'Select the Trajectory file');
                obj.Traj_filename=Traj_filename;
            end
            tic
            
            fileTransInd_name=strcat( 'BM.mat') ;
            try
                load(fileTransInd_name, 'BM') ;
                obj.BM=BM ;
            catch
                %                     load(fileTransInd_name, 'BelongTransM')
                %                     BM=BelongTransM ;
            end
            %             BM=BelongTransM ;
            %                 obj.BM=BM ;
            
            %             sdfsf=3
            %-------  read topology file
            delimiterIn = ' ';
            A= importdata(strcat(PathName1,obj.Topology_filename),delimiterIn,1) ;
            
            FirstRow=strsplit(A.textdata{1,1});
            NBase=str2double(FirstRow{1});
            % NStrand=str2double(FirstRow{2});
            obj.TSeq= A.textdata((2:NBase+1)',2);
            
            QQ=A.textdata((2:NBase+1)',1);
            TStrand=zeros(size(QQ));
            for k=1:length(QQ)
                TStrand(k)=str2double(QQ{k});
            end
            obj.Primes35 = A.data;
            %
            %             %-------- read json file to get color
            %                         [Json_filename,PathName3,FilterIndex3]= uigetfile({'*.json','json file' },'Select the json file');
            %                         obj.JsonFile=Json_filename;
            %
            %                         dat=loadjson(strcat(PathName3, obj.JsonFile));
            %                         ColorCode=[-1 ,-1 ,-1];
            %                         EffecCyl=[];
            %                         for k=1:length(dat.vstrands)
            %                         if ~isempty(dat.vstrands{k}.stap_colors)
            %                         ColorCode=[ColorCode;   [k*ones(size(dat.vstrands{k}.stap_colors,1),1), dat.vstrands{k}.stap_colors   ]];
            %                         EffecCyl=union(EffecCyl,k);
            %                         end
            %                         end
            %                         ColorCode=setdiff(ColorCode,[-1,-1,-1],'rows') ;
            %                         [~,obj.colorRGB]=hist(ColorCode(:,3),unique(ColorCode(:,3)));
            %                         fprintf('Color codes are from JSON files. \n')
            %                         nBundle=length(obj.colorRGB);
            
            %             obj.colorRGB=rand(20,3) ;
            
            %------------get box size
            fffid=fopen(strcat(PathName2,Traj_filename));
            Oneline = fgetl(fffid);
            Oneline = fgetl(fffid);
            obj.Boxs=   str2num(Oneline(4:end)) ;
            fclose(fffid);
            %-------
            
            
            try
                %                     cd(PathName1);
                %                 load('AllTraj.mat','AllTraj')
                load(strcat(PathName1,'AllTraj.mat'),'AllTraj')
                %                 sss=ggg(3);   % intentional error to trigger catch
                %                                 load('AllTrajxxxx.mat','AllTraj')
                %                 AllTraj=AllTraj(:,:,1:end-1) ;  % ignore pca frame
                N_frame=size(AllTraj,3) ;
                
                %----------
                %                     First500Traj =AllTraj ;
                %                     [AllTrajname,PathAllTraj,~]= uigetfile({'*.mat','Add Traj file' },'Select the added AllTraj file');
                %
                %                     load(strcat(PathAllTraj,AllTrajname),'AllTraj')
                %                     Second1000Traj =AllTraj ;
                %                     Second1000Traj=Second1000Traj(:,:,2:end) ;
                %                     AllTraj= cat(3,First500Traj,Second1000Traj);
                %                     fprintf('hard coding for reading two saved trajectories.\n') ;
                %                     N_frame=size(AllTraj,3) ;
                %                     AllTraj=single(AllTraj) ;  % signle-digit precision
                %                     save('AllTraj.mat','AllTraj', '-v7.3')
                
                %-----
                
                %                     cd(obj.Spwd)
                %                 AllTraj  = FcnCancelDrift( AllTraj,obj.Boxs(1) ) ;
                
                %                load('PcaTraj.mat','PcaTraj') ;
                %                obj.PcaTraj=PcaTraj ;
                fprintf(' Found trajectory file in mat.  \n'  )
                
            catch
                fprintf(' could not find trajectory file in mat. searching one line by one \n'  )
                format3='%f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64 %f64' ;
                fffid=fopen(strcat(PathName2,Traj_filename));
                totalRow=0;
                while 1   % get totalRow for allocate,
                    Oneline = fgetl(fffid) ;
                    if    strcmp(Oneline(1),'t')
                        Oneline = fgetl(fffid); % header 2
                        Oneline = fgetl(fffid);  % header 3
                        C3 = textscan(fffid,format3,NBase);
                        Oneline = fgetl(fffid) ;
                        totalRow=totalRow+NBase +3;
                    else
                        break
                    end
                end
                fclose(fffid);
                
                N_frame= totalRow/(NBase+3) ;
                AllTraj=zeros(NBase, 9 ,  N_frame) ;
                fffid2nd=fopen(strcat(PathName2,Traj_filename));
                %     Oneline = fgetl(fffid2nd);
                count=0;
                nloop=0;
                whereist=[];
                %     extraFrame=40 ;
                %                     frame=1;
                frame=1; nw=0;
                while 1
                    Oneline = fgetl(fffid2nd) ;
                    count=count+1;
                    if    strcmp(Oneline(1),'t')
                        frame=frame+1;
                        Oneline = fgetl(fffid2nd); count=count+1;  % header 2
                        Oneline = fgetl(fffid2nd); count=count+1; % header 3
                        C3 = textscan(fffid2nd,format3,NBase);
                        A3=cell2mat(C3) ;
                        AllTraj( : , : ,frame )=A3(:,1:9) ;
                        %                             fprintf('GetOneMoreConf \n')
                        Oneline = fgetl(fffid2nd) ;
                    else
                        break
                    end
                    if nw>20000
                        break
                    end
                    nw=nw+1;
                end
                fclose(fffid2nd);
                
                N_frame=frame;
                toc
                %------------------prova configuration
                
                fffid3rd=fopen(strcat(PathName2,provaConf));
                count=0;
                nloop=0;
                whereist=[];
                %     extraFrame=40 ;
                frame=0;
                %                     tic
                while 1
                    Oneline = fgetl(fffid3rd); count=count+1;
                    if strcmp(Oneline(1),'t')
                        frame=frame+1;
                        Oneline = fgetl(fffid3rd); count=count+1;  % header 2
                        Oneline = fgetl(fffid3rd); count=count+1; % header 3
                        locaR=1;
                        continue
                    end
                    
                    if  length(Oneline)==1
                        if Oneline==-1
                            fclose(fffid3rd);
                            break
                        end
                    end
                    OneRowdata=str2num(Oneline) ;
                    AllTraj( locaR , : ,frame )=OneRowdata(1:9) ;
                    locaR=locaR+1;
                end
                
                
                %
                %                 AllTraj  = FcnCancelDrift( AllTraj,obj.Boxs(1) ) ;
                
                %%
                
                %                 N_frame=N_frame+1;
%                                     AllTraj=single(AllTraj) ;  % signle-digit precision
%                                 save('AllTraj.mat','AllTraj', '-v7.3')
                toc
            end  % end of try
            
            %               AllTraj  = FcnCancelDrift( AllTraj,obj.Boxs(1) ) ;
            %                         [a,b]=hist(TStrand,unique(TStrand)) ;a0=a  ;b0=b;
            %                         scafInitNotFirst=0;
            %
            %                         if 1~= find(a==max(a))
            %                             sacfInd=find(a==max(a));
            %                             QQT=AllTraj;
            %                             QQTStrand=TStrand;
            %                             PrevStappleInd = sum(a(1:sacfInd-1))  ;
            %                             ScafIndAll=PrevStappleInd+1:a(sacfInd)+PrevStappleInd  ;
            %                             AllTraj(ScafIndAll,:,:)=[];
            %                             TStrand(ScafIndAll,:)=[];
            %                             AllTraj=[QQT(ScafIndAll,:,:) ;AllTraj];
            %                             TStrand=[QQTStrand(ScafIndAll,:) ;TStrand];
            %                             QQQ=TStrand;
            %                             [a1,~,~]=unique(TStrand,'stable') ;
            %                             %     IndOrder= unique(TStrand,'stable') ;
            %                             for subs=1:length(a1)
            %                             TStrand(QQQ==a1(subs))=subs;
            %                             end
            %                             scafInitNotFirst=1;
            %                         end
%             size(AllTraj)
% obj.Boxs=[600 600 600]
            obj.Strand=TStrand ; CutLine =0.5 ;
%             AllTraj  = FcnCancelDrift( AllTraj,obj.Boxs(1) ) ;
            for k = 1: size(AllTraj ,3)
            XYZOneConf = AllTraj(:,1:3,k ) ;
            
% % %             if k==1
%              XYZOneConf(:,1)  = rem(XYZOneConf(:,1) , obj.Boxs(1)) ;
%              XYZOneConf(:,2)  = rem(XYZOneConf(:,2) , obj.Boxs(2)) ;
%              XYZOneConf(:,3)  = rem(XYZOneConf(:,3) , obj.Boxs(3)) ;
% %             else
%              XYZOneConf(:,1)  = mod(XYZOneConf(:,1)-0.5*obj.Boxs(1)  , obj.Boxs(1)) ;
%              XYZOneConf(:,3)  = mod(XYZOneConf(:,3)-0.5*obj.Boxs(3) , obj.Boxs(3)) ;
  %    

        
           if max(XYZOneConf(:,1)) - min(XYZOneConf(:,1)) >obj.Boxs(1)
             XYZOneConf(:,1)  = mod(XYZOneConf(:,1)-CutLine*obj.Boxs(1) , obj.Boxs(1)) ;
           end
%           Inds = XYZOneConf(:,1)<50 ;  % hard
%          XYZOneConf(Inds,1)=XYZOneConf(Inds,1)+obj.Boxs(1) ;

           if max(XYZOneConf(:,2)) - min(XYZOneConf(:,2)) >obj.Boxs(2)
             XYZOneConf(:,2)  = mod(XYZOneConf(:,2)-CutLine*obj.Boxs(2) , obj.Boxs(2)) ;
           end

           if max(XYZOneConf(:,3)) - min(XYZOneConf(:,3)) >obj.Boxs(3)
             XYZOneConf(:,3)  = mod(XYZOneConf(:,3)-CutLine*obj.Boxs(3) , obj.Boxs(3)) ;
           end

             
%              XYZOneConf(:,1)  = rem(XYZOneConf(:,1)-min(XYZOneConf(:,1))  , obj.Boxs(1)) ;
%              XYZOneConf(:,2)  = rem(XYZOneConf(:,2)-min(XYZOneConf(:,2)), obj.Boxs(2)) ;
%              XYZOneConf(:,3)  = rem(XYZOneConf(:,3)-min(XYZOneConf(:,3)) , obj.Boxs(3)) ;

             
%                 
% %             end
             
%              XYZOneConf(:,1)  = XYZOneConf(:,1)  - obj.Boxs(1)*ceil(XYZOneConf(:,1) /obj.Boxs(1) ) ;
%             XYZOneConf(:,2)  = XYZOneConf(:,2)  - obj.Boxs(2)*ceil(XYZOneConf(:,2) /obj.Boxs(2) ) ;
%             XYZOneConf(:,3)  = XYZOneConf(:,3)  - obj.Boxs(3)*ceil(XYZOneConf(:,3) /obj.Boxs(3) ) ;
            AllTraj(:,1:3,k )= XYZOneConf ;            
            end
            fprintf('Assume periodic boundary condition. Minimun image convention \n')
            
            
            obj.N_frame=N_frame ;
            Maxbox=max(AllTraj(:,1:3),[],3) ;
            Minbox=min(AllTraj(:,1:3),[],3) ;
            Center= mean(0.5*(Maxbox+Minbox)) ;
            
           ScafInd = obj.Strand==1 ;
            for k=1:1  %N_frame
               
%                 AllTraj(:,1:3,k)=AllTraj(:,1:3,k)- Center ;
                
                
                 ScafCenter =   mean(AllTraj(ScafInd,1:3,k) ) ;
                 AllTraj(:,1:3,k)=AllTraj(:,1:3,k)- ScafCenter ; %TdT 36
                 
                 
                %                 AllTraj(:,1:3,k)=AllTraj(:,1:3,k)*0.85 ;
                %
                %                   for strandi = 1 :max(obj.Strand)
                %                       BaseInd = obj.Strand== strandi ;
                %                       AllTraj(BaseInd,1:9,k) = flip(  AllTraj(BaseInd,1:9,k) ,1) ;
                %
                %                      if k ==1
                %                          obj.Primes35(BaseInd ,:) =flip( obj.Primes35(BaseInd ,:)  ,1) ;
                %                          obj.TSeq(BaseInd ,:) =flip( obj.TSeq(BaseInd ,:)  ,1) ;
                %                      end
                %
                %                   end
                
            end
            
%                         T1= AllTraj(:,1:3,1) ;
%                          for k=2:N_frame
%                            [regParams,~,~]=absor(AllTraj(:,1:3,k)' ,T1');
%             %                regParams.R*AllTraj(:,1:3,k)'+ regParams.t -T1' ;
%                            AllTraj(:,1:3,k)= transpose(regParams.R*AllTraj(:,1:3,k)'+ regParams.t ) ;
%                            AllTraj(:,4:6,k)= transpose(regParams.R*AllTraj(:,4:6,k)' ) ;
%                          end
%             
%             
            
            
            obj.Traj=AllTraj ;
            cd(obj.Spwd);
            
            %             obj.getPcaTraj ;
            
        end
        
        function ReaddSRemain(obj)
            
            fffid=fopen(strcat(obj.PathName,'dSRemain.conf'));
            if fffid==-1
                [DsFile,PathName2,FilterIndex]= uigetfile({'*.dat;*.conf','Mutual trap file' },'Select the mutual trap file');
                fffid=fopen(strcat(PathName2,DsFile));
            else
                PathName2=obj.PathName ;
                DsFile='dSRemain.conf';
            end
            
            
            n =0; EnterZone = false;
            while 1
                Oneline = fgetl(fffid);
                if strcmp(Oneline,'{')
                    EnterZone=true ;
                end
                if strcmp(Oneline,'}') && EnterZone
                    EnterZone=false ;
                    n=n+1 ;
                end
                if Oneline==-1
                    fclose(fffid);
                    break
                end
            end
            
            PairList = zeros(n,2) ;
            fffid2=fopen(strcat(PathName2,DsFile));
            n =1; EnterZone = false;
            while 1
                Oneline = fgetl(fffid2);
                if Oneline==-1
                    fclose(fffid2);
                    break
                    
                end
                if strcmp(Oneline,'{')
                    EnterZone=true ;
                end
                
                if contains(Oneline,'particle') && ~contains(Oneline,'ref_particle') && EnterZone
                    pp = find(Oneline == '=') ;
                    PairList(n,1) = str2num(Oneline(pp+1:end)) ;
                end
                if contains(Oneline,'ref_particle') && EnterZone
                    pp = find(Oneline == '=') ;
                    PairList(n,2) = str2num(Oneline(pp+1:end)) ;
                end
                
                if strcmp(Oneline,'}') && EnterZone
                    EnterZone=false ;
                    n=n+1 ;
                end
            end
            PairList=PairList+1 ;
            SwtichInd = PairList(:,2) <PairList(:,1) ;
            PairList(SwtichInd,: ) =flip(PairList(SwtichInd,: ), 2) ;
            PairList= unique(PairList, 'rows') ;
            
            obj.MutualTrapPair  =PairList  ;
            
            %             fclose(fffid);
            
            
        end
        
        function getPcaTraj2(obj)
            % without specifying two nodes as X-dir, direct compare with
            % the first frame using absor
            
            
            
            
            IndsInBundle = obj.BM==4 ;
            
            %------------align along a bundle
            Conf1_BundleAlign =  obj.Traj(IndsInBundle,:,1) ;
            [coeff,~,~,~,~,mu] = pca(Conf1_BundleAlign(:,1:3))  ;
            Ref_X =  coeff(:,1)'     ;
            [coeff2,~,~,~,~,mu] =pca(obj.Traj(:,1:3,1) ) ;
            Ref_Z =  coeff2(:,3)'    ;
            Ref_Z= cross(  cross( Ref_X,Ref_Z),Ref_X ) ;  % make sure these two are perpendicular
            R = [Ref_X ; cross(Ref_Z,Ref_X) ; Ref_Z]';
            obj.Traj(:,1:3,1)= obj.Traj(:,1:3,1)*R ;
            obj.Traj(:,4:6,1)= obj.Traj(:,4:6,1)*R ;
            obj.Traj(:,7:9,1)= obj.Traj(:,7:9,1)*R ;
            %-------------
            
            IndsInBundle= obj.BM>0 ;
            PcaTraj=zeros(size(obj.Traj)) ; PcaTraj(:,:,1)= obj.Traj(:,:,1) ; T1=PcaTraj(:,1:3,1) ;
            Traj= obj.Traj ;
            for k =2 :1 :obj.N_frame
                %              for k=2:N_frame
                [regParams,~,~]=absor(Traj(IndsInBundle,1:3,k)' ,T1(IndsInBundle,:)');
                %                regParams.R*AllTraj(:,1:3,k)'+ regParams.t -T1' ;
                T_absor13= transpose(regParams.R*Traj(:,1:3,k)'+ regParams.t ) ;
                T_absor46= transpose(regParams.R*Traj(:,4:6,k)' ) ;
                T_absor79= transpose(regParams.R*Traj(:,7:9,k)' ) ;
                
                %              end
                PcaTraj(:,1:9,k) =[T_absor13 ,T_absor46, T_absor79 ] ;
                
            end
            obj.PcaTraj=PcaTraj ;
            obj.Traj=PcaTraj ;
        end
        
        function getPcaTraj(obj)
            f2=figure(2) ;
            if isempty(f2.UserData)
                fprintf( ' Please use function visualConfAndAssignRef and mouse to assign two bases as ref \n' ) ;    return
            end
            if isempty(f2.UserData.RightClickInd) || isempty(f2.UserData.LeftClickInd)
                fprintf( ' Please use function visualConfAndAssignRef and mouse to assign two bases as ref \n' );    return
            end
            %-------------
            PcaTraj=zeros(size(obj.Traj)) ;
            pcaThirdSave=[];
            for iF=1 : obj.N_frame
                ConfT =obj.Traj(:,1:3, iF) ;
                [coeff,~,~,~,~,mu] = pca(ConfT)  ;
                PcaThirdVec=  coeff(:,3)' ;
                RefVecByTwoBases= ConfT(f2.UserData.RightClickInd,1:3) - ConfT(f2.UserData.LeftClickInd,1:3) ;
                RefVecByTwoBases=RefVecByTwoBases/norm(RefVecByTwoBases) ;
                pcaThirdA= dot(RefVecByTwoBases,PcaThirdVec) ;
                
                if isempty(pcaThirdSave)
                    
                    AbsRefZ=PcaThirdVec;
                    pcaThirdSave= PcaThirdVec ;  %only update in the first time
                else
                    if dot(pcaThirdSave,PcaThirdVec)>0
                        AbsRefZ=PcaThirdVec ;
                    else
                        AbsRefZ=-PcaThirdVec ;
                    end
                    pcaThirdSave= PcaThirdVec ;
                end
                
                
                %                 if pcaThirdA>0
                %                 AbsRefZ=PcaThirdVec;
                %                 else
                %                 AbsRefZ=-PcaThirdVec ;
                %                 end
                
                AbsRefX= RefVecByTwoBases- dot(RefVecByTwoBases,AbsRefZ)*AbsRefZ ;
                
                AbsRefY= cross(AbsRefZ,AbsRefX) ;
                RotM = [AbsRefX;AbsRefY;AbsRefZ ]' ;
                
                PcaConfT=zeros(size(ConfT,1) ,9) ;
                PcaConfT(:,1:3) =  (ConfT-mu)*RotM ;
                PcaConfT(:,4:6) =  obj.Traj(:,4:6, iF)*RotM ;
                PcaConfT(:,7:9) =  obj.Traj(:,7:9, iF)*RotM ;
                
                %                 scatter3(PcaConfT(:,1),PcaConfT(:,2) ,PcaConfT(:,3));
                %                 for Ppj=1 : 3    % checked
                %                 XYZ=[mu;mu +  50*coeff(:,Ppj)'] ;
                % %                 plot3( XYZ(:,1),XYZ(:,2),XYZ(:,3),'Linewidth',5)
                %                 end
                %                 XYZ=[mu;mu +  50*AbsRefX] ; plot3( XYZ(:,1),XYZ(:,2),XYZ(:,3),'r','Linewidth',5)
                %                 XYZ=[mu;mu +  50*AbsRefZ] ; plot3( XYZ(:,1),XYZ(:,2),XYZ(:,3),'b','Linewidth',5)
                
                PcaTraj(:,:,iF) =PcaConfT ;
                
            end
            obj.PcaTraj=PcaTraj ;
            
            
        end
        
        %         function visualConf(obj,frame)
        %              if nargin==1;frame=1;end
        %             f31=figure(31);clf; hold on ;
        %             f31.WindowScrollWheelFcn= @(src,evn)autoaxis(obj,src,evn) ;
        %
        %             T=obj.Traj(:,:,frame) ;
        %             cd(obj.PathName);
        % %             load('BM.mat')  ;  % read BelongTransM from previous export
        %             Coeff=-0.4;
        %             BackBoneT=  T(:,1:3) +Coeff* T(:,4:6) ;
        %             [a0,b0]=hist(obj.Strand,unique(obj.Strand)) ;
        %             for strandi=1:max(obj.Strand)
        %                 BaseInd=  find(obj.Strand==strandi) ;
        %                  if   a0(strandi)==max(a0)   %scaffold ribbon
        %                       plot3(BackBoneT(BaseInd,1) , BackBoneT(BaseInd,2) ,BackBoneT(BaseInd,3)  ,'Color', [0.2,0.5,0.8]) ;
        %                  else   %staple
        % %                         BundleInd=  unique(BelongTransM(BaseInd)) ;
        % %                         rgbDec=obj.colorRGB( BundleInd) ;
        % %                         RGBHex=dec2hex(rgbDec,6 ) ;
        % %                         Rc=hex2dec(RGBHex(1:2));
        % %                         Gc=hex2dec(RGBHex(3:4)) ;
        % %                         Bc=hex2dec(RGBHex(5:6))  ;
        % %                         BundleColor=[Rc,Gc ,  Bc]/256 ;
        %
        %                         plot3(BackBoneT(BaseInd,1) , BackBoneT(BaseInd,2) ,BackBoneT(BaseInd,3) ) ;
        %                  end
        %             end
        %             xlabel('x') ;ylabel('y') ;zlabel('z') ; axis equal ;
        %             title(strcat('Frame = ' ,num2str(frame)   )  )
        %             cd(obj.Spwd) ;
        %         end
        function pH=visualTraj(obj,frame,mode)
            if nargin==1;frame=1;end
            if   mode~=1 && isempty(  obj.PcaTraj)
                fprintf(' Need to do PCA first !!  \n') ;
                commandwindow
                return
            end
            %             sdfsf=1
            
            %                          load('BM.mat')  ;  % read BelongTransM from previous export
            
            %                         f11=figure(11);clf; hold on ;
            %                         f11.WindowScrollWheelFcn= @(src,evn)autoaxis(obj,src,evn) ;
            %             cd(obj.PathName);
            %             load('BM.mat')  ;  % read BelongTransM from previous export
            Coeff=-0.4;   % default -0.4 as oxDNA website,
            
            scatterCell= cell(length(frame),1) ;
            cc_scaf=0;  cc_stap=0;
            for iF=1:length(frame)
                %                 if mode==1
                T=obj.Traj(:,:,frame(iF)) ;
%                 T= obj.NoDriftTraj(:,:,frame(iF)) ;
                
                a2 = cross(T(:,7:9),T(:,4:6) ) ;
                BackBoneT = T(:,1:3)  +obj.POS_MM_BACK1*T(:,4:6) + obj.POS_MM_BACK2*a2 ;
                
                %                 POS_MM_BACK1 = -0.3400 ;
                %                 POS_MM_BACK2 = 0.3408 ;
                
                
                %                 BackBoneT=  T(:,1:3) +Coeff* T(:,4:6) ;
                [a0,b0]=hist(obj.Strand,unique(obj.Strand)) ;
                CC=[iF,0, length(frame)-iF]/ length(frame) ;
                pH=cell(length(frame),1) ;
                for strandi=1:max(obj.Strand)
                    BaseInd=  find(obj.Strand==strandi) ;
                    if  a0(strandi)>2000    %scaffold ribbon   CColor =  get(0,'defaultAxesColorOrder') ;   a0(strandi)==max(a0)
                        pH{strandi}= plot3(BackBoneT(BaseInd,1)+iF*0  , BackBoneT(BaseInd,2)+iF*0 ,BackBoneT(BaseInd,3) , 'b') ;
                        pH{strandi}.UserData.RenderAs = 'Scaf' ;  cc_scaf=cc_scaf+1; pH{strandi}.UserData.Ind=cc_scaf;
                    else   %staple
                        pH{strandi}=   plot3(BackBoneT(BaseInd,1)+iF*0 , BackBoneT(BaseInd,2)+iF*0 ,BackBoneT(BaseInd,3),'r'  ) ;
                        pH{strandi}.UserData.RenderAs = 'Stap' ; cc_stap=cc_stap+1; pH{strandi}.UserData.Ind=cc_stap;
                    end
                    pH{strandi}.UserData.oxDNAStrandInd = strandi ;
                    pH{strandi}.UserData.IniXYZ =BackBoneT(BaseInd,1:3) ;
                    pH{strandi}.UserData.BaseID = BaseInd ;
                end
                
                for strandi=1:max(obj.Strand)
                    pH{strandi}.UserData.pH=pH;
                end
                
                %               scatterCell{iF} = scatter3(T(obj.BM==1,1 )+iF*50,T(obj.BM==1,2 ),T(obj.BM==1,3 ),24,'o','filled'   ) ;
                
            end
            
            %             Inds= obj.Strand==38; sum(Inds)            % used to find which strand disociate
            %             scatter3(T(Inds,1)+Coeff* T(Inds,4),T(Inds,2)+Coeff* T(Inds,5),T(Inds,3)+Coeff* T(Inds,6) ,'filled')
            
            xlabel('x') ;ylabel('y') ;zlabel('z') ; axis equal ;
            title(strcat('Frame = ' ,num2str(frame)   )  )   ;
            %             legend( )
            %             cd(obj.Spwd) ;
        end
        function visualConfAndAssignRef(obj,frame)
            if nargin==1;frame=1;end
            f2=figure(2);clf;
            hold on ;
            f2.UserData.RightClick=[];            f2.UserData.LeftClick=[];
            f2.UserData.RightClickInd=[];            f2.UserData.LeftClickInd=[];
            
            f2.WindowScrollWheelFcn= @(src,evn)autoaxis(obj,src,evn) ;
            T=obj.Traj(:,:,frame) ;
            cd(obj.PathName);
            load('BM.mat')  ;  % read BelongTransM from previous export
            Coeff=-0.4;
            BackBoneT=  T(:,1:3) +Coeff* T(:,4:6) ;
            [a0,b0]=hist(obj.Strand,unique(obj.Strand)) ;
            scaterH=cell(1,max(obj.Strand)   ) ;
            for strandi=1:max(obj.Strand)
                BaseInd=  find(obj.Strand==strandi) ;
                if   a0(strandi)==max(a0)   %scaffold ribbon
                    scaterH{strandi}=scatter3(BackBoneT(BaseInd,1) , BackBoneT(BaseInd,2) ,BackBoneT(BaseInd,3),24  , [0.2,0.5,0.8] ,'.' ) ;
                    
                else   %staple
                    BundleInd=  unique(BelongTransM(BaseInd)) ;
                    rgbDec=obj.colorRGB( BundleInd) ;
                    RGBHex=dec2hex(rgbDec,6 ) ;
                    Rc=hex2dec(RGBHex(1:2));
                    Gc=hex2dec(RGBHex(3:4)) ;
                    Bc=hex2dec(RGBHex(5:6))  ;
                    BundleColor=[Rc,Gc ,  Bc]/256 ;
                    
                    scaterH{strandi}= scatter3(BackBoneT(BaseInd,1) , BackBoneT(BaseInd,2) ,BackBoneT(BaseInd,3),24  ,BundleColor ,'.') ;
                end
                scaterH{strandi}.ButtonDownFcn=@(src,evn)Hitscatter(obj,src,evn,f2) ;
                scaterH{strandi}.UserData=BaseInd ;
            end
            xlabel('x') ;ylabel('y') ;zlabel('z') ; axis equal ;
            title(strcat('Frame = ' ,num2str(frame)   )  )
            cd(obj.Spwd) ;
        end
        
        function  autoaxis(obj,src,evn)
            if evn.VerticalScrollCount==1
                axis equal ;
            elseif evn.VerticalScrollCount==3
                axis auto;
            end
        end
        
        
        function Hitscatter(obj,src,evn,f2)
            PotentialHitsBaseInd= src.UserData ;
            XYZList= [src.XData ; src.YData  ; src.ZData ]' ;
            dXYZ= XYZList-evn.IntersectionPoint ;
            d=dXYZ(:,1).^2 + dXYZ(:,2).^2 +dXYZ(:,3).^2  ;
            IndSelect= find(d==min(d)) ;
            switch  evn.Button
                case 1
                    
                    if isempty(f2.UserData.LeftClick)
                        f2.UserData.LeftClick= scatter3(XYZList(IndSelect,1),XYZList(IndSelect,2),XYZList(IndSelect,3),82,'or','filled' ) ;
                    else
                        f2.UserData.LeftClick.XData=XYZList(IndSelect,1) ;    f2.UserData.LeftClick.YData=XYZList(IndSelect,2) ;    f2.UserData.LeftClick.ZData=XYZList(IndSelect,3) ;
                    end
                    f2.UserData.LeftClickInd=PotentialHitsBaseInd(IndSelect) ;
                    
                case 3
                    if isempty(f2.UserData.RightClick)
                        f2.UserData.RightClick= scatter3(XYZList(IndSelect,1),XYZList(IndSelect,2),XYZList(IndSelect,3),82,'ob','filled' ) ;
                    else
                        f2.UserData.RightClick.XData=XYZList(IndSelect,1) ;    f2.UserData.RightClick.YData=XYZList(IndSelect,2) ;    f2.UserData.RightClick.ZData=XYZList(IndSelect,3) ;
                    end
                    f2.UserData.RightClickInd=PotentialHitsBaseInd(IndSelect) ;
            end
        end

        function ExportLastConf(obj)
            
            LastConf = obj.Traj(:,1:9 ,end) ;
%             LastConf(:,1:3) =   LastConf(:,1:3) -mean(  LastConf(:,1:3)) + [12 12 12 ] ; % hard
%             boxsize =obj.Boxs ;
            
            LastConf(:,1:3) = LastConf(:,1:3) - mean(LastConf(:,1:3)) ; 
            MMbox = max( LastConf(:,1:3)) ;
            mmbox = min( LastConf(:,1:3)) ;
            boxsize = 1.5 * max(MMbox- mmbox) ;
                       
%             boxsize=25 ; % hard
            fprintf('Printing the last configuration at %s \n',obj.PathName );
            cd(obj.PathName);       
            file2_name='lastConf_1e8.conf';
            fileID = fopen(file2_name,'w');
            fprintf(fileID,'t = 10000000\n'); E0=0;
            fprintf(fileID,'b = %9.6f %9.6f %9.6f\n',max(boxsize),max(boxsize),max(boxsize));
            fprintf(fileID,'E = %8.6f %8.6f %8.6f\n',E0,E0,E0);            
            for k=1:  size(LastConf,1)
                fprintf(fileID,'%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f\n',LastConf(k,:),zeros(1,6) );  % changed NVec for primes of oxDNA
            end
            
            fclose(fileID) ;
            cd(obj.Spwd) ;
        end
        
        
        function ExportBildRibbonTraj(obj,Ribbon,iF,filename)
            
            if nargin==1;Ribbon=0;end
            tic
            Maxbox=max(max(obj.Traj(:,1:3),[],3)) ; Maxbox=Maxbox(1:3) ;Maxbox=[max(Maxbox),max(Maxbox),max(Maxbox)] ;
            Minbox=min(min(obj.Traj(:,1:3),[],3)) ; Minbox=Minbox(1:3) ;Minbox=[max(Minbox),max(Minbox),max(Minbox)] ;
            [a0,b0]=hist(obj.Strand,unique(obj.Strand)) ;
            ColorRGB=[0.2,0.5,0.8];            %scaffold backbone
            Coeff=-0.2;
            cd(obj.PathName);            [status, msg, msgID] = mkdir('BILDs') ;
            
            %             load('BM.mat')  ;  % read BelongTransM from previous export
            cd([obj.PathName 'BILDs'] );
            %------------
            %                 ColorCode = (obj.BaseColor- min(obj.BaseColor) )/(max(obj.BaseColor)-min(obj.BaseColor))* [1,0,0];
            %                 ColorCode(:,3)=1-ColorCode(:,1) ;
            %------------
            PrintTraj=obj.Traj;
            T=PrintTraj(:,:,iF) ;
            
            if  Ribbon==1  %only ribbon
                file2_name=strcat(filename,'_ribbon_',num2str(iF-1) , '.bild') ;
                fprintf(' printing bild (ribbon) file %i \n' ,iF)
            else
                file2_name=strcat(filename,'_CG_',num2str(iF-1) , '.bild') ;
                fprintf(' printing bild (CG) file %i \n' ,iF)
            end
            
            fileID = fopen([obj.PathName 'BILDs' filesep file2_name],'w');
            countbase =1;
            BelongTransM=obj.BM ;
            
            %             randRGB=rand(1,3)
            %             randRGB=[1,1,1];
            CColor =  get(0,'defaultAxesColorOrder') ;
            for strandi=1:  max(obj.Strand)
                BaseInd=  find(obj.Strand==strandi) ;
                %                     if strandi==2
                %                         ColorRGB = [1,0,1];
                %                     end
                %                     strandi
                %                    CC  get(0,'defaultAxesColorOrder')
                if      a0(strandi)>1000 %    a0(strandi)>1000  % longer than 6000 nts                a0(strandi)==max(a0)   %scaffold ribbon
                    %                         ColorRGB= CColor( mod(strandi,size(CColor,1))+1  ,:) ;
                    %                         scafShow=11
                    ColorRGB=[0.2,0.5,0.8]+0.15*rand(1,3);            %scaffold backbone
                    %                     ColorRGB = [0.5 0.5 0.5 ] ;
                    %                                              ColorRGB=[0.8,0,0] ;
                    %                     ColorRGB= randRGB ;
                    fprintf('Scaffold length = %i \n', max(a0)  ) ;
                    fprintf(fileID , '\n'  )    ;
                    fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                    for kbs= 1:length(BaseInd)
                        
                        
                        
                        if Ribbon==1  %only ribbon
                            %                             fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                            
                            if kbs==length(BaseInd)
                                fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',BackBoneB,1.1*obj.SphereR/2 )    ;
                                continue ;
                            end  %# ribbon =N-1
                            BackBoneA=  T(BaseInd(kbs),1:3) +Coeff* T(BaseInd(kbs),4:6) ;
                            BackBoneB=  T(BaseInd(kbs+1),1:3) +Coeff* T(BaseInd(kbs+1),4:6) ;
                            %                             fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,obj.SphereR/2 )    ;
                            fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,1.1*obj.SphereR/2 )    ;
                        else        %with sphere on backbone and cylinder between
                            %                             fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                            %                                     ColorRGB=ColorCode(countbase,:) ; countbase=countbase+1;
                            
                            BackBoneA=  T(BaseInd(kbs),1:3) +Coeff* T(BaseInd(kbs),4:6) ;
                            BaseA= T(BaseInd(kbs),1:3)-Coeff* T(BaseInd(kbs),4:6) ;
                            %                             fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BaseA,obj.CylR/2 )    ;
                            fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,obj.SphereR/2 )    ;
                            %                                     fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',[1,0,0]  )    ;
                            %                                     fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n \n',BaseA,0.5*obj.SphereR )    ;
                            
                            if kbs==length(BaseInd)
                                continue ;
                            end  %# ribbon =N-1
                            BackBoneB=  T(BaseInd(kbs+1),1:3) +Coeff* T(BaseInd(kbs+1),4:6) ;
                            %                             fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,obj.CylR/2 )    ;
                        end
                    end
                else  %staple
                    
                    %                                                 BundleInd=  unique(BelongTransM(BaseInd)) ;
                    %                                                 if length(BundleInd) >1
                    %                                                     N2=countCM(BelongTransM(BaseInd) ,BundleInd)  ;
                    %                                                    BundleInd = BundleInd(N2==max(N2)) ;
                    %                                                 end
                    
                    %                                                 BundleInd=BundleInd(1) ;
                    %                                                 rgbDec=obj.colorRGB( BundleInd) ;
                    %                                                 RGBHex=dec2hex(rgbDec,6 ) ;
                    %                                                 Rc=hex2dec(RGBHex(1:2));
                    %                                                 Gc=hex2dec(RGBHex(3:4)) ;
                    %                                                 Bc=hex2dec(RGBHex(5:6))  ;
                    %                                                 BundleColor=[Rc,Gc ,  Bc]/256 ;
%                     if length(BaseInd) == 36 % extended
%                         BundleColor=[248,144,31]/255 ;
%                     elseif  length(BaseInd) == 26 % shrink
%                         BundleColor=[0.3,0.3,0.3] ;
%                     else
%                         BundleColor=[0.5,0,0] ;
%                     end
                    
                    
                    BundleColor=[0.8,0,0] ;
                    %                     BundleColor=[0.5 0.5 0.5] ;
                    %                                                 BundleColor= rand(1,3) ;
                    fprintf(fileID , '\n'  )    ;
                    fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',BundleColor+0*(0.5-rand(1,3))  )    ;
                    for kbs= 1:length(BaseInd)
                        
                        if Ribbon==1  %only ribbon
                            if kbs==length(BaseInd)
                                fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n ',BackBoneB,1.1*obj.SphereR/2 )    ;
                                continue ;
                            end  %# ribbon =N-1
                            BackBoneA=  T(BaseInd(kbs),1:3) +Coeff* T(BaseInd(kbs),4:6) ;
                            BackBoneB=  T(BaseInd(kbs+1),1:3) +Coeff* T(BaseInd(kbs+1),4:6) ;
                            %                             BaseA= T(BaseInd(kbs),1:3)-Coeff* T(BaseInd(kbs),4:6) ;
                            %                           fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',BundleColor  )    ;
                            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,obj.SphereR/2 )    ;
                            fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,1.2*obj.SphereR/2 )    ;
                        else  %with sphere on backbone and cylinder between
                            
                            %                                     BundleColor=ColorCode(countbase,:) ; countbase=countbase+1;
                            
                            
                            
                            BackBoneA=  T(BaseInd(kbs),1:3) +Coeff* T(BaseInd(kbs),4:6) ;
                            BaseA= T(BaseInd(kbs),1:3)-Coeff* T(BaseInd(kbs),4:6) ;
                            %                             fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',BundleColor  )    ;
                            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BaseA,obj.CylR/2 )    ;
                            fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,obj.SphereR/2 )    ;
                            if kbs==length(BaseInd); continue ;end  %# ribbon =N-1
                            BackBoneB=  T(BaseInd(kbs+1),1:3) +Coeff* T(BaseInd(kbs+1),4:6) ;
                            %                             fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',BundleColor  )    ;
                            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,obj.CylR/2 )    ;
                        end
                    end
                end
            end   %strandi
            fclose(fileID);
            %             end
            fprintf(' finish printing bild ribbon file \n' )
            cd(obj.Spwd);
            toc
        end
        
        function ExportBildRibbonTraj_Ori(obj,Ribbon)
            
            if nargin==1;Ribbon=0;end
            Ribbon
            tic
            Maxbox=max(max(obj.Traj(:,1:3),[],3)) ; Maxbox=Maxbox(1:3) ;Maxbox=[max(Maxbox),max(Maxbox),max(Maxbox)] ;
            Minbox=min(min(obj.Traj(:,1:3),[],3)) ; Minbox=Minbox(1:3) ;Minbox=[max(Minbox),max(Minbox),max(Minbox)] ;
            [a0,b0]=hist(obj.Strand,unique(obj.Strand)) ;
            ColorRGB=[0.2,0.5,0.8];            %scaffold backbone
            Coeff=-0.4;
            cd(obj.PathName);            [status, msg, msgID] = mkdir('BILDs') ;
            
            load('BM.mat')  ;  % read BelongTransM from previous export
            cd([obj.PathName 'BILDs'] );
            if ~isempty(obj.PcaTraj)
                PrintTraj=obj.PcaTraj;
            else
                PrintTraj=obj.Traj;
            end
            
            for iF=1:1:obj.N_frame  %        1:1000:obj.N_frame-1
                T=PrintTraj(:,:,iF) ;
                %                                 T=obj.Traj(:,:,iF) ;
                
                %                 shiftCemter =  mean(T(:,1:3)) - 0.5*(Maxbox+Minbox) ;
                %                 T(:,1:3)=T(:,1:3) - shiftCemter;
                
                if  Ribbon==1  %only ribbon
                    file2_name=strcat(obj.objname,'_ribbon_',num2str(iF) , '.bild') ;
                    fprintf(' printing bild (ribbon) file %i \n' ,iF)
                else
                    file2_name=strcat(obj.objname,'_CG_',num2str(iF) , '.bild') ;
                    fprintf(' printing bild (CG) file %i \n' ,iF)
                end
                
                fileID = fopen([obj.PathName 'BILDs\' file2_name],'w');
                for strandi=1:  max(obj.Strand)
                    BaseInd=  find(obj.Strand==strandi) ;
                    if   a0(strandi)==max(a0)   %scaffold ribbon
                        
                        for kbs= 1:length(BaseInd)
                            if Ribbon==1  %only ribbon
                                if kbs==length(BaseInd)
                                    fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n \n',BackBoneB,obj.SphereR/2 )    ;
                                    continue ;
                                end  %# ribbon =N-1
                                BackBoneA=  T(BaseInd(kbs),1:3) +Coeff* T(BaseInd(kbs),4:6) ;
                                BackBoneB=  T(BaseInd(kbs+1),1:3) +Coeff* T(BaseInd(kbs+1),4:6) ;
                                fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                                fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,obj.SphereR/2 )    ;
                                fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n \n',BackBoneA,obj.SphereR/2 )    ;
                            else        %with sphere on backbone and cylinder between
                                BackBoneA=  T(BaseInd(kbs),1:3) +Coeff* T(BaseInd(kbs),4:6) ;
                                BaseA= T(BaseInd(kbs),1:3)-Coeff* T(BaseInd(kbs),4:6) ;
                                fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                                fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BaseA,obj.CylR )    ;
                                fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n \n',BackBoneA,obj.SphereR )    ;
                                
                                if kbs==length(BaseInd); continue ;end  %# ribbon =N-1
                                BackBoneB=  T(BaseInd(kbs+1),1:3) +Coeff* T(BaseInd(kbs+1),4:6) ;
                                fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorRGB  )    ;
                                fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,obj.CylR )    ;
                            end
                        end
                    else  %staple
                        
                        %                         BundleInd=  unique(BelongTransM(BaseInd)) ;
                        %                         rgbDec=obj.colorRGB( BundleInd) ;
                        %                         RGBHex=dec2hex(rgbDec,6 ) ;
                        %                         Rc=hex2dec(RGBHex(1:2));
                        %                         Gc=hex2dec(RGBHex(3:4)) ;
                        %                         Bc=hex2dec(RGBHex(5:6))  ;
                        %                         BundleColor=[Rc,Gc ,  Bc]/256 ;
                        BundleColor=[1,0,0] ;
                        
                        for kbs= 1:length(BaseInd)
                            if Ribbon==1  %only ribbon
                                if kbs==length(BaseInd)
                                    fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n \n',BackBoneB,obj.SphereR/2 ) ;
                                    continue ;
                                end  %# ribbon =N-1
                                BackBoneA=  T(BaseInd(kbs),1:3) +Coeff* T(BaseInd(kbs),4:6) ;
                                BackBoneB=  T(BaseInd(kbs+1),1:3) +Coeff* T(BaseInd(kbs+1),4:6) ;
                                %                             BaseA= T(BaseInd(kbs),1:3)-Coeff* T(BaseInd(kbs),4:6) ;
                                fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',BundleColor  )    ;
                                fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,obj.SphereR/2 )    ;
                                fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n \n',BackBoneA,obj.SphereR/2 )    ;
                            else  %with sphere on backbone and cylinder between
                                
                                BackBoneA=  T(BaseInd(kbs),1:3) +Coeff* T(BaseInd(kbs),4:6) ;
                                BaseA= T(BaseInd(kbs),1:3)-Coeff* T(BaseInd(kbs),4:6) ;
                                fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',BundleColor  )    ;
                                fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BaseA,obj.CylR )    ;
                                fprintf(fileID , '.sphere %4.2f %4.2f %4.2f %4.2f \n \n',BackBoneA,obj.SphereR )    ;
                                if kbs==length(BaseInd); continue ;end  %# ribbon =N-1
                                BackBoneB=  T(BaseInd(kbs+1),1:3) +Coeff* T(BaseInd(kbs+1),4:6) ;
                                fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',BundleColor  )    ;
                                fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,obj.CylR )    ;
                            end
                        end
                    end
                end   %strandi
                fclose(fileID);
            end
            fprintf(' finish printing bild ribbon file \n' )
            cd(obj.Spwd);
            toc
        end
        function printCompChimera(obj,mode )
            if nargin==1;mode=1;end
            if mode==1
                str='_ribbon_';
            else
                str='_CG_';
            end
            obj.PcaTraj=obj.Traj ;
            %get  box
            Maxbox=max(max(obj.PcaTraj(:,1:3),[],3)) ; Maxbox=Maxbox(1:3) ;Maxbox=[max(Maxbox),max(Maxbox),max(Maxbox)] ;
            Minbox=min(min(obj.PcaTraj(:,1:3),[],3)) ; Minbox=Minbox(1:3) ;Minbox=[min(Minbox),min(Minbox),min(Minbox)] ;
            
            %             Maxbox=max(max(obj.Traj(:,1:3),[],3)) ; Maxbox=Maxbox(1:3) ;Maxbox=[max(Maxbox),max(Maxbox),max(Maxbox)] ;
            %             Minbox=min(min(obj.Traj(:,1:3),[],3)) ; Minbox=Minbox(1:3) ;Minbox=[min(Minbox),min(Minbox),min(Minbox)] ;
            
            
            file2B_name=strcat('TrajBox','.bild') ;
            fileID1 = fopen([obj.PathName 'BILDs\' file2B_name],'w');
            
            %             fprintf(fileID1 , '.color %8.6f %8.6f %8.6f \n',[0,0,0]  )    ;
            fprintf(fileID1 , '.color %8.6f %8.6f %8.6f \n',[0,0,0.8]  )    ;
            
            fprintf(fileID1 , '.transparency  %4.2f \n',0  )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',Minbox,[Minbox(1),Minbox(2),Maxbox(3)] , 1 )    ;
            fprintf(fileID1 , '.color %8.6f %8.6f %8.6f \n',[0.8,0,0]  )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',Minbox,[Maxbox(1),Minbox(2),Minbox(3)] , 1 )    ;
            fprintf(fileID1 , '.color %8.6f %8.6f %8.6f \n',[0,0.8,0]  )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',Minbox,[Minbox(1),Maxbox(2),Minbox(3)] , 1 )    ;
            fprintf(fileID1 , '.color %8.6f %8.6f %8.6f \n',[0,0,0]  )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Maxbox(1),Minbox(2),Maxbox(3)],[Maxbox(1),Minbox(2),Minbox(3)] , 1 )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Maxbox(1),Minbox(2),Maxbox(3)],[Maxbox(1),Maxbox(2),Maxbox(3)] , 1 )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Maxbox(1),Minbox(2),Maxbox(3)],[Minbox(1),Minbox(2),Maxbox(3)] , 1 )    ;
            
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Maxbox(1),Maxbox(2),Minbox(3)],[Maxbox(1),Maxbox(2),Maxbox(3)] , 1 )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Maxbox(1),Maxbox(2),Minbox(3)],[Maxbox(1),Minbox(2),Minbox(3)] , 1 )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Maxbox(1),Maxbox(2),Minbox(3)],[Minbox(1),Maxbox(2),Minbox(3)] , 1 )    ;
            
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Minbox(1),Maxbox(2),Maxbox(3)],[Maxbox(1),Maxbox(2),Maxbox(3)] , 1 )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Minbox(1),Maxbox(2),Maxbox(3)],[Minbox(1),Minbox(2),Maxbox(3)] , 1 )    ;
            fprintf(fileID1 , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',[Minbox(1),Maxbox(2),Maxbox(3)],[Minbox(1),Maxbox(2),Minbox(3)] , 1 )    ;
            fclose(fileID1);
            
            %--------------
            file2_name=strcat('AskChimeraExportPng ',str(1:end-1) , '.com') ;
            fileID = fopen([obj.PathName file2_name],'w');
            
            fprintf(fileID , 'close session  \n' )    ;
            fprintf(fileID , 'cd %s  \n' , [obj.PathName 'BILDs\'])    ;
            
            fprintf(fileID , 'set projection orthographic   \n')    ;
            fprintf(fileID , 'set bgColor white   \n')    ;
            fprintf(fileID , 'open %s  \n', file2B_name)    ;
            
            %             fprintf(fileID , 'scale 1.6 \n')    ;
            fprintf(fileID , 'scale 1.6 \n')    ;
            
            fprintf(fileID , 'windowsize 1920 360 \n')    ;
            %                         fprintf(fileID , 'windowsize 960 960 \n')    ;
            
            fprintf(fileID , '~modeldisplay  #0 \n')    ;
            
            %             fprintf(fileID , '2dlab create struct5 text char(39)actuating to Triangle 1char(39) size 30 xpos .5 ypos 0.15 color .0,.05,.02 \n')    ;
            %              fprintf(fileID, '2dlab create struct2 text %sActuating to Triangle 1%s size 30 xpos .4 ypos 0.25 color .0,.05,.02 \n' ,char(39),char(39))    ;
            %              fprintf(fileID, '2dlab create struct1 text %sClosing strands bonded %s size 30 xpos .4 ypos 0.17 color .8,.8,.1 \n',char(39),char(39) )    ;
            
            %             fprintf(fileID, '2dlab create struct2 text %sAfter adding closing strands %s size 30 xpos .4 ypos 0.25 color .0,.05,.02 \n' ,char(39),char(39))    ;
            %              fprintf(fileID, '2dlab create struct1 text %sClosing strands unbonded, free rotatation %s size 30 xpos .4 ypos 0.17 color .8,.8,.1 \n',char(39),char(39) )    ;
            
            
            fprintf(fileID , ' \n')    ;
            for k=1:1:obj.N_frame
                %                 bild_name=strcat('Wbb_ribbon_',num2str(k) , '.bild') ; str
                %                 bild_name=strcat(obj.objname,'_ribbon_',num2str(k) , '.bild') ;
                bild_name=strcat(obj.objname,str,num2str(k) , '.bild') ;
                
                fprintf(fileID , 'open %s  \n', bild_name)    ;
                if k==1
                    fprintf(fileID , 'scale 2.8 \n')    ;
                end
                fprintf(fileID , 'turn y  models #1 90   \n')    ;
                %                fprintf(fileID , 'turn z  models #1 90   \n')    ;
                %                 fprintf(fileID , 'turn z  models #1 180   \n')    ;
                fprintf(fileID , 'turn x  models #1 180   \n')    ;
                png_name=  strcat('image_',num2str(k) , '_.png') ;
                %                 fprintf(fileID , 'wait   \n')    ;
                fprintf(fileID , 'copy file %s  \n',  [obj.PathName 'BILDs\' png_name])    ;
                %----side view (2nd view)
                %                 fprintf(fileID , 'turn y  models #1 90   \n')    ;
                %                 png_name2=  strcat('Sideimage_',num2str(k) , '_.png') ;
                %                 %                 fprintf(fileID , 'wait   \n')    ;
                %                 fprintf(fileID , 'copy file %s  \n',  [obj.PathName 'BILDs\' png_name2])    ;
                %
                %                 %----top view
                %                 fprintf(fileID , 'turn y  models #1 -90   \n')    ;
                fprintf(fileID , 'turn x  models #1 90   \n')    ;
                png_name3=  strcat('Topimage_',num2str(k) , '_.png') ;
                fprintf(fileID , 'copy file %s  \n',  [obj.PathName 'BILDs\' str png_name3])    ;
                %
                
                
                if k~=obj.N_frame
                    fprintf(fileID , 'close #1 \n')    ;    fprintf(fileID , ' \n')    ;
                end
            end
            
            %                fprintf(fileID , 'close session  \n' )    ;  %----------------optional
            %             fprintf(fileID , 'stop  \n' )    ;  %----------------optional
            fclose(fileID);
            
            fprintf('finish printCompChimera \n'  )
        end %end of printCompChimera
        
        function askChimeraRotate(obj,mode )
            if nargin==1;mode=1;end
            if mode==1
                str='_ribbon_';
            else
                str='_CG_';
            end
            %-------------------------
            
            
            file2_name=strcat('AskChimeraRotate ',str(1:end-1) , '.com') ;
            fileID = fopen([obj.PathName file2_name],'w');
            
            fprintf(fileID , 'close session  \n' )    ;
            fprintf(fileID , 'cd %s  \n' , [obj.PathName 'BILDs\'])    ;
            
            fprintf(fileID , 'set projection orthographic   \n')    ;
            fprintf(fileID , 'set bgColor white   \n')    ;
            
            %             fprintf(fileID , 'scale 1.6 \n')    ;
            fprintf(fileID , 'scale 1.0 \n')    ;
            
            fprintf(fileID , 'windowsize 1920 1080 \n')    ;
            
            fprintf(fileID , ' \n')    ;
            
            bild_name=strcat(obj.objname,str,num2str(1) , '.bild') ;
            fprintf(fileID , 'open %s  \n', bild_name)    ;
            
            
            for k=0:10:360
                
                
                
                %                fprintf(fileID , 'turn x  models #1 -90   \n')    ;
                %                fprintf(fileID , 'turn z  models #1 90   \n')    ;
                %                 fprintf(fileID , 'turn z  models #1 180   \n')    ;
                fprintf(fileID , 'turn x  models #1 %i   \n',k  )    ;
                png_name=  strcat('image_',num2str(k) , '_.png') ;
                %                 fprintf(fileID , 'wait   \n')    ;
                fprintf(fileID , 'copy file %s  \n',  [obj.PathName 'BILDs\' png_name])    ;
                %                 %----side view (2nd view)
                %                  fprintf(fileID , 'turn y  models #1 90   \n')    ;
                %                 png_name2=  strcat('Sideimage_',num2str(k) , '_.png') ;
                % %                 fprintf(fileID , 'wait   \n')    ;
                %                 fprintf(fileID , 'copy file %s  \n',  [obj.PathName 'BILDs\' png_name2])    ;
                %
                %                 %----top view
                %                    fprintf(fileID , 'turn y  models #1 -90   \n')    ;
                %                    fprintf(fileID , 'turn x  models #1 90   \n')    ;
                %                    png_name3=  strcat('Topimage_',num2str(k) , '_.png') ;
                %                   fprintf(fileID , 'copy file %s  \n',  [obj.PathName 'BILDs\' str png_name3])    ;
                %
                %                 if k~=obj.N_frame
                %                 fprintf(fileID , 'close #1 \n')    ;    fprintf(fileID , ' \n')    ;
                %                 end
            end
            
            %                fprintf(fileID , 'close session  \n' )    ;  %----------------optional
            %             fprintf(fileID , 'stop  \n' )    ;  %----------------optional
            fclose(fileID);
            
            fprintf('finish askChimeraRotate \n'  )
            
            
        end
        
    end
    
end

