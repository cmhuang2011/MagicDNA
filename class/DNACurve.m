classdef DNACurve <handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % use for post-algorithm, may use hyperbundle and freeformcurve class
    properties
        
        Ori_freeformcurve = [];  %
        Scaled_FFC = [];  %
        %         st_line_model = [];  %
        
        MagicDNA_fH =[] ;
        MagicDNA_design = [] ; % hyperbundle
        
        initial_oxDNA= [] ;
        oxDNA_varients = [];
        BaseC5inoxDNA=[];
        %--------
        MatrixForObjFunc =[] ; % use extrudePerfectModel to calculate,  %[C5 ,CylPos, ds_NVec, ds_ei ]
        BaseC5_inoxDNA_sBase = [];
        %-----grachics
        f_align =[];
        ph_Domain = [] ;
        Domain_info =[];
        ph_PerfectRoutingStap = [];
        ph_PerfectRoutingScaf = [];
        ph_BaseOrientation = [];
        ph_BaseOrient_domain= [];
        ph_SuggestStapXover= [];
        ph_SelectXover= []; 
        %-----Improve insertion and skip
        AddInsertBase=[];
        AddSkipBase =[] ;
        
    end
    
    methods
        function ExampleOnCommandline(obj)
            Curve1 = DNACurve(fH) ;
            Curve1.visulize_initial ;
            
            Curve1.extrudePerfectModel_V2 ;
            Curve1.DrawDomainOnPerfectModel ;
%             Curve1.DrawRoutingOnPerfectModel
            Curve1.DrawBaseOrientationOnPerfectModel
            
            Curve1.UpdateBundleSkip_Insert ;
            %-------
            Curve1.MagicDNA_design.getDomainGraph ;
            Domains =Curve1.MagicDNA_design.DuplexDomain ;
            
            Curve1.MagicDNA_design.DomainGraph.neighbors(85)
            highlight(h,85 ,'NodeColor',[1 0 0])
            highlight(h, Curve1.MagicDNA_design.DomainGraph.neighbors(85) ,'NodeColor',[0 1 0])
            
            
        end
        
        function UpdateBundleSkip_Insert(obj)
            obj.MagicDNA_design.CustomSkipAndInsertion = true ;
            
            DomainL = obj.Domain_info ;
            CurrentDomainR=obj.MagicDNA_design.DuplexDomain;
            
            DomainModL = DomainL(:,1)*0.34-DomainL(:,2) ;
            Threshold= [ 0.3, -0.3 ]  ; % upper and lower bounds, nm
            TooLong = find(DomainModL>Threshold(1) );
            TooShort = find(DomainModL<Threshold(2) );
            VeryFewBaseDomain = find(DomainL(:,1)<5 ) ;
            size(TooLong)
            TooLong=setdiff(TooLong,VeryFewBaseDomain) ;
             size(TooLong)
            
            AddSkip =[-1 -1];           %C5
            for dm = 1: size(TooLong,1)
                DoLengthi =  size(  CurrentDomainR{TooLong(dm)}  ,1    ) ;
                ToolongFor_nm = DomainL(TooLong(dm),1)*0.34-DomainL(TooLong(dm),2) ;
                N_Skip =  round(ToolongFor_nm/0.34) ;
                Inds=  round(linspace(2,DoLengthi-1 ,N_Skip)) ;
%                 if Inds==0
%                     sdf=1
%                 end
%                 
                AddSkip= [AddSkip ;  CurrentDomainR{TooLong(dm)}(Inds,:) ]  ;
            end
            AddSkip=AddSkip(2:end,:) 
            %--------
            AddInsert =[-1 -1];           %C5
            for dm = 1: size(TooShort,1)
                DoLengthi =  size(  CurrentDomainR{TooShort(dm)}  ,1    ) ;
                ToolongFor_nm = -DomainL(TooShort(dm),1)*0.34+DomainL(TooShort(dm),2) ;
                N_Skip =  round(ToolongFor_nm/0.34) ;
                Inds=  round(linspace(2,DoLengthi-1 ,N_Skip)) ;
                AddInsert= [AddInsert ;  CurrentDomainR{TooShort(dm)}(Inds,:) ]  ;
            end
            AddInsert=AddInsert(2:end,:) ;
            %------
            AddInsert =setdiff(AddInsert , AddSkip,'rows') ;
            %--------
            Rtable =obj.MagicDNA_design.RelateTable ;
            for  k=1:size(AddSkip)
                [~, Ind] = ismember( AddSkip(k,1) , Rtable(:,5) ) ;
                BC = Rtable(Ind,1:2) ;
                obj.MagicDNA_design.containBundle{BC(1)}.skipPosition(end+1,:) =[BC , AddSkip(k,1:2)] ;
            end
            for k=1:length( obj.MagicDNA_design.containBundle)
%                 [k, size(obj.MagicDNA_design.containBundle{k}.skipPosition)]
                obj.MagicDNA_design.containBundle{k}.skipPosition=unique(obj.MagicDNA_design.containBundle{k}.skipPosition,'rows') ;
%                [k, size(obj.MagicDNA_design.containBundle{k}.skipPosition)]
            end
            
            for  k=1:size(AddInsert)
                [~, Ind] = ismember( AddInsert(k,1) , Rtable(:,5) ) ;
                BC = Rtable(Ind,1:2) ;
                obj.MagicDNA_design.containBundle{BC(1)}.InsertPosition(end+1,:) =[BC , AddInsert(k,1:2)] ;
            end
            for k=1:length( obj.MagicDNA_design.containBundle)
%                  [k, size(obj.MagicDNA_design.containBundle{k}.skipPosition)]
                obj.MagicDNA_design.containBundle{k}.InsertPosition=unique(obj.MagicDNA_design.containBundle{k}.InsertPosition,'rows') ;
%                 [k, size(obj.MagicDNA_design.containBundle{k}.skipPosition)]
            end
        end
        
        function DrawRoutingOnPerfectModel(obj)
            GetHyperB=obj.MagicDNA_design ;
            Domains=GetHyperB.DuplexDomain;
            scafC5= GetHyperB.scafC5 ;
            stapC5 =GetHyperB.StapList3 ;   %stapC5=cell2mat(stapC5) ;
            StapAllBase=GetHyperB.StapAllBase ;
            
            obj.ph_PerfectRoutingScaf=gobjects(length(stapC5) ,1) ;
            for k=1:length(scafC5)
                C5_base = scafC5{k} ;
                [aa, indexOnTable] = ismember( C5_base , obj.MatrixForObjFunc(:,1:2), 'rows') ;
                indexOnTable(indexOnTable==0)=[];  % may have single-stranded scaffold
                XYZ = obj.MatrixForObjFunc(indexOnTable,3:5) ;
                %                 obj.ph_PerfectRoutingScaf(k)=plot3(XYZ(:,1) ,XYZ(:,2) , XYZ(:,3) ,'-k','LineWidth',1) ;
            end
            
            obj.ph_PerfectRoutingStap=gobjects(length(stapC5) ,1) ;
            for k=1:length(stapC5)
                C5_base = StapAllBase{k} ;  % change to all base to have smooth curve
                [aa, indexOnTable] = ismember( C5_base , obj.MatrixForObjFunc(:,1:2), 'rows') ;
                XYZ = obj.MatrixForObjFunc(indexOnTable,3:5) ;
                obj.ph_PerfectRoutingStap(k)=plot3(XYZ(:,1) ,XYZ(:,2) , XYZ(:,3) ,'--k','LineWidth',1,'HitTest','off') ;
            end
            uistack(obj.ph_Domain,'top') ;
        end
        
        function DrawBaseOrientationOnPerfectModel(obj)
            GetHyperB=obj.MagicDNA_design ;
            [Domains, ~]=GetHyperB.DuplexDomain_neglectNickFcn;
            %             Domains=GetHyperB.DuplexDomain;
            scafC5= GetHyperB.scafC5 ;    ScafpAllBase=GetHyperB.ScafAllBase ;
            stapC5 =GetHyperB.StapList3 ; StapAllBase=GetHyperB.StapAllBase ;
            
            Isstap=0 ;  TM=1 ;
            [ScafHelixXYZ,ScafBasesCenter  ]=plotScaf2_Helix_V2_noGraphics( GetHyperB,scafC5,0 ,[0,0,1] ,TM ) ;     % get stap strands coordinate
            [StapHelixXYZ,StapBasesCenter  ]=plotScaf2_Helix_V2_noGraphics( GetHyperB,stapC5,1 ,[0,0,1] ,TM ) ;     % get stap strands coordinate
            
            Rtable= obj.MagicDNA_design.RelateTable;
            
            scafC5=cell2mat(scafC5) ;     stapC5=cell2mat(stapC5) ;
            ScafHelixXYZ=cell2mat(ScafHelixXYZ') ;    ScafBasesCenter=cell2mat(ScafBasesCenter') ;
            StapHelixXYZ=cell2mat(StapHelixXYZ') ;    StapBasesCenter=cell2mat(StapBasesCenter') ;
            ScafpAllBase=cell2mat(ScafpAllBase);       StapAllBase=cell2mat(StapAllBase);
            
            % MatrixForObjFunc: [C5 ,CylPos, ds_NVec, ds_ei ]
            Cons =1 ; IdealTwist =  360/10.5 ;
            UVW = zeros(10000,3) ; c= 1;
            UVW_domain = zeros(100000,3) ; c2= 1;
            DomainL_threshold = 24 ;
            DomainInd_AboveThres = 1:length(Domains) ;
            for k=1: length(Domains)
                ThisDomain = Domains{k} ;
                [~, indexOnTable] = ismember( ThisDomain([1 end],:) , obj.MatrixForObjFunc(:,1:2), 'rows') ;
                XYZ = obj.MatrixForObjFunc(indexOnTable,3:5) ;
                [aa2,IndInStap] = ismember( ThisDomain([1 end],:) , StapAllBase, 'rows')  ;%IndInStap(IndInStap==0) =[];
                [aa3,IndInScaf] = ismember( ThisDomain([1 end],:) , ScafpAllBase, 'rows')  ;%IndInStap(IndInStap==0) =[];
                %                 [aa2,IndInStap] = ismember( ThisDomain([1 end],:) , StapAllBase, 'rows')  ;%IndInStap(IndInStap==0) =[];
                %                 [aa3,IndInScaf] = ismember( ThisDomain([1 end],:) , ScafpAllBase, 'rows')  ;%IndInStap(IndInStap==0) =[];
                
                [~,indd]=ismember(ThisDomain(:,1),obj.MagicDNA_design.RelateTable(:,5) ) ;
                BC_rep = obj.MagicDNA_design.RelateTable(indd,1:2) ;
                if length(unique(BC_rep(:,2)) )>1
                    ShowError= 12345
                else
                    Cyl= unique(BC_rep(:,2)) ;
                    obj.ph_Domain(k).UserData.Cyl=Cyl ;
                end
                
                if IndInStap(1)==1;  st1 =true;
                else;  st1= StapAllBase(IndInStap(1),1)~= StapAllBase(IndInStap(1)-1,1); end
                
                if   st1      %aa2(1)==1
                    HeadOrentation = StapHelixXYZ(IndInStap(1) , :)- StapBasesCenter(IndInStap(1) , :)  ;
                elseif ScafpAllBase(IndInScaf(1),1)~= ScafpAllBase(IndInScaf(1)-1,1)
                    HeadOrentation = -(ScafHelixXYZ(IndInScaf(1) , :)- ScafBasesCenter(IndInScaf(1) , :))  ;
                else
                    ShowError =11;
                end
                
                if IndInStap(2)==size(StapAllBase,1);  st2=true;
                else;  st2=StapAllBase(IndInStap(2),1)~= StapAllBase(IndInStap(2)+1,1);  end
                
                
                if st2   %aa2(2)==1
                    TailOrentation = StapHelixXYZ(IndInStap(2) , :)- StapBasesCenter(IndInStap(2) , :)  ;
                elseif aa3(2)==1
                    TailOrentation = ScafHelixXYZ(IndInScaf(2) , :)- ScafBasesCenter(IndInScaf(2) , :)  ;
                else
                    ShowError =11;
                end
                %                 if k==86
                %                     eds=2
                %                 end
                
                Orientation =vec3norm_CM([HeadOrentation ; TailOrentation]) ;
                
                AddArrs = [XYZ(1,:) ; XYZ(1,:)+Cons*Orientation(1,:) ; nan nan nan  ;...
                    XYZ(2,:) ; XYZ(2,:)+Cons*Orientation(2,:) ; nan nan nan    ] ;
                UVW(c:c+5,:) =  AddArrs ; c=c+6 ;
                %--------------interpolate base orientation
                % MatrixForObjFunc: [C5 ,CylPos, ds_NVec, ds_ei ]
                %                 t_domain =
                [~ , EntireDomain_ind ]=  ismember( ThisDomain, obj.MatrixForObjFunc(:,1:2), 'rows') ;
                PVec=  obj.MatrixForObjFunc(EntireDomain_ind,3:5) ;
                %----equalize distance
                ssline = cscvn( PVec') ;
                PVec_eq=fnval(ssline, linspace(ssline.breaks(1),ssline.breaks(end),size(PVec,1)  ))';
                
                TVec = vec3norm_CM(gradient(PVec')') ;
                
                AllSegOnThisCyl_C5s =  Rtable(Rtable(:,2)==Cyl,5) ;
                IndAllBase= ismember(  obj.MatrixForObjFunc(:,1),AllSegOnThisCyl_C5s);
                
                XYZonCly= obj.MatrixForObjFunc(IndAllBase,3:5 );
                [~,DoInd_OnCyl ]= ismember(PVec , XYZonCly,'rows')  ;
                splineOfThiCylinder= cscvn(  obj.MatrixForObjFunc(IndAllBase,3:5 )') ;
                DomainBases_t = splineOfThiCylinder.breaks(DoInd_OnCyl) ;
                t_Corr_onfine = DomainBases_t/max(splineOfThiCylinder.breaks)*max(obj.Ori_freeformcurve.t_arrs_fine) ;
                Inds= round(interp1(obj.Ori_freeformcurve.t_arrs_fine,1:length(obj.Ori_freeformcurve.t_arrs_fine),t_Corr_onfine) );
                Ref_xy = obj.Ori_freeformcurve.PTN_Vec_onCurve(Inds,7:9) ;
                
                if size(ThisDomain,1) < DomainL_threshold
                    DomainInd_AboveThres=setdiff(DomainInd_AboveThres,k) ;
                    continue;
                end
                
                %-------
                a2Vec=Orientation(1,:)- dot(Orientation(1,:),TVec(1,:))*TVec(1,:);  a2Vec=a2Vec/norm(a2Vec);
                Refa2Vec=Ref_xy(1,:)- dot(Ref_xy(1,:),TVec(1,:))*TVec(1,:);  Refa2Vec=Refa2Vec/norm(Refa2Vec);
                if dot (cross(a2Vec,Refa2Vec),TVec(1,:))>0
                    Angle_Head = acosd( dot(a2Vec,Refa2Vec))  ;
                else
                    Angle_Head = -acosd( dot(a2Vec,Refa2Vec))  ;
                end
                
                b2Vec=Orientation(2,:)- dot(Orientation(2,:),TVec(end,:))*TVec(end,:);  b2Vec=b2Vec/norm(b2Vec);
                Refb2Vec=Ref_xy(end,:)- dot(Ref_xy(end,:),TVec(end,:))*TVec(end,:);  Refb2Vec=Refb2Vec/norm(Refb2Vec);
                if dot (cross(b2Vec,Refa2Vec),TVec(end,:))>0
                    Angle_Tail = acosd( dot(b2Vec,Refb2Vec))  ;
                else
                    Angle_Tail = -acosd( dot(b2Vec,Refb2Vec))  ;
                end
                Angle_Tail=-Angle_Tail ;
                Angle_Head=-Angle_Head ;
                
                
                TargetTailAng= mod(Angle_Head+ IdealTwist*(size(ThisDomain,1)-1) , 360) ;
                if mod(Angle_Tail-TargetTailAng,360)> mod(TargetTailAng-Angle_Tail,360)
                    CompAng=  -mod(TargetTailAng-Angle_Tail,360) ;
                else
                    CompAng=  mod(Angle_Tail-TargetTailAng,360) ;
                end
                EndAngle= Angle_Head+ IdealTwist*(size(ThisDomain,1)-1)+CompAng ;
                AngArr= linspace(Angle_Head, EndAngle,size(ThisDomain,1)) ;
                %                 mod(AngArr,360)
                %                 diff(  mod(AngArr,360))
                
                BaseOrientation=[Orientation(1,:);zeros(length(ThisDomain)-2 ,3) ; Orientation(2,:)  ] ;
                for bbs =2:length(ThisDomain)-1
                    RAxis =TVec(bbs,:) ;
                    theta = mod(AngArr(bbs),360) ;
                    RMat = RotationAxisToRotaionMatrix( -RAxis,theta ) ;
                    BaseOrientation(bbs,:)=  Ref_xy(bbs,:)*RMat ;
                end
                
                
                XYZall=PVec_eq ;
                AddArrsAll =zeros(3*size(BaseOrientation,1),3   ) ;
                for qqi = 1: length(ThisDomain)
                    AddArrsAll(3*qqi-2:3*qqi ,:) =[XYZall(qqi,:) ; XYZall(qqi,:)+0.8*Cons*BaseOrientation(qqi,:) ; nan nan nan ] ;
                end
                
                AddArrsAll=AddArrsAll(2:3:end ,:) ;  %helical Rep
                AddArrsAll(end+1,:) =[nan nan nan] ;
                obj.ph_Domain(k).UserData.BaseBackBone= AddArrsAll(1:end-1,:) ;
                
                
                %                 AddArrs = [XYZall(1,:) ; XYZall(1,:)+Cons*Orientation(1,:) ; nan nan nan  ;...
                %                 XYZall(2,:) ; XYZall(2,:)+Cons*Orientation(2,:) ; nan nan nan    ] ;
                %                 UVW(c:c+5,:) =  AddArrs ; c=c+6 ;
                %              plot3(AddArrsAll(:,1) ,AddArrsAll(:,2) ,AddArrsAll(:,3) ,'g' );
                
                UVW_domain(c2:c2+size(AddArrsAll,1)-1 ,:) =  AddArrsAll ; c2=c2+size(AddArrsAll,1)  ;
                
                %                 sdfsdf=3
                %                                             RMat = RotationAxisToRotaionMatrix( LocalAxis,theta*180/pi )     ;
                
                
            end
            
            UVW=UVW(1:c-1 ,:) ;
            UVW_domain=UVW_domain(1:c2-1 ,:) ;
            if isempty(  obj.ph_BaseOrientation) || ~isvalid( obj.ph_BaseOrientation)
                obj.ph_BaseOrientation = plot3(UVW(:,1) ,UVW(:,2) ,UVW(:,3) ,'r' ,'HitTest','off');
            else
                obj.ph_BaseOrientation.XData=UVW(:,1)';
                obj.ph_BaseOrientation.YData=UVW(:,2)';
                obj.ph_BaseOrientation.ZData=UVW(:,3)';
            end
            
            if isempty(  obj.ph_BaseOrient_domain) || ~isvalid( obj.ph_BaseOrient_domain)
                obj.ph_BaseOrient_domain = plot3(UVW_domain(:,1) ,UVW_domain(:,2) ,UVW_domain(:,3) ,'g' ,'LineWidth',1.5,'HitTest','off');
            else
                obj.ph_BaseOrient_domain.XData=UVW_domain(:,1)';
                obj.ph_BaseOrient_domain.YData=UVW_domain(:,2)';
                obj.ph_BaseOrient_domain.ZData=UVW_domain(:,3)';
            end
            
            set( obj.ph_Domain ,'Color',[0.2 0.2 0.2]  );  set( obj.ph_Domain ,'LineWidth',0.5  ) ;
            set( obj.ph_Domain(DomainInd_AboveThres) ,'Color',[0 0 0.8]  );  set( obj.ph_Domain(DomainInd_AboveThres) ,'LineWidth',2  ) ;
            
            
            uistack(obj.ph_Domain,'top') ;
        end
        
        function DrawDomainOnPerfectModel(obj)
            
            GetHyperB=obj.MagicDNA_design ;
%             [Domains, ~]=GetHyperB.DuplexDomain_neglectNickFcn;
                          Domains=GetHyperB.DuplexDomain;
            %             MergeTableInOri
            figure(obj.f_align) ;
            obj.ph_Domain= gobjects(length(Domains) ,1) ;
            obj.Domain_info = zeros(length(Domains) ,2 ) ; % [NBase , length ]
            for k= 1:  length(Domains)
                C5_base = Domains{k} ;
                [aa, indexOnTable] = ismember( C5_base , obj.MatrixForObjFunc(:,1:2), 'rows') ;
                if ~all(aa)
                    sdfsf=2
                end
                XYZ = obj.MatrixForObjFunc(indexOnTable,3:5) ;
                
                obj.ph_Domain(k)=plot3(XYZ(:,1) ,XYZ(:,2) , XYZ(:,3) ,'LineWidth',1) ;
                d = diff(XYZ);
                total_length = sum(sqrt(sum(d.*d,2))) +0.34;
                obj.Domain_info(k,:) = [size(XYZ,1),total_length  ]  ;
                
            end
            for k=1:length(Domains)
                obj.ph_Domain(k).UserData.ind =k ;
                obj.ph_Domain(k).UserData.C5Base =Domains{k} ;
            end
            set(  obj.ph_Domain, 'ButtonDownFcn', @(src,evn)SelectDomain(obj,src, evn)) ;
            %               end
            
            IgonreShortDomain =find( obj.Domain_info(:,1)<5  )  ;SS=  obj.Domain_info(:,1)<5 ;
            EffectDomain = find(~SS) ;
            
            figure(153);clf ; subplot(2,1,1) ;
            p1= plot(EffectDomain,obj.Domain_info(~SS,1)*0.34,'-x') ; hold on ;
            p2 = plot(EffectDomain,obj.Domain_info(~SS,2),'-o') ;
            
            legend([ p1 , p2] ,{ 'NBase*0.34','L on PerfectModel ' }) ;
            subplot(2,1,2) ;
            p3= plot(EffectDomain,obj.Domain_info(~SS,1)*0.34-obj.Domain_info(~SS,2),'-x')  ; hold on ;
            legend(p3, 'Difference (nm)') ;
            if isempty(obj.f_align.KeyPressFcn)
                set(obj.f_align,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) );
            end
            %------------
            %             R =sum(obj.Domain_info(:,1)*0.34)./sum(obj.Domain_info(:,2) ) ;
            %             obj.Domain_info(:,2) =R*obj.Domain_info(:,2)  ; % scale perfect model such that mean value =0
            
            DomainModL = obj.Domain_info(:,1)*0.34-obj.Domain_info(:,2) ;
            Threshold= [ 0.3, -0.3 ]  ; % upper and lower bounds, nm
            TooLong = find(DomainModL>Threshold(1) ); TooLong=setdiff(TooLong, IgonreShortDomain) ;
            TooShort = find(DomainModL<Threshold(2) ); TooShort=setdiff(TooShort, IgonreShortDomain) ;
            Other = setdiff(setdiff([1:length(obj.ph_Domain)],TooLong  ) ,TooShort) ;
            
            set( obj.ph_Domain(Other) ,'Color',[0.2 0.2 0.2]  );  set( obj.ph_Domain(Other) ,'LineWidth',0.2  ) ;
            set( obj.ph_Domain(TooLong) ,'Color',[1 0 0]  );  set( obj.ph_Domain(TooLong) ,'LineWidth',2  ) ;
            set( obj.ph_Domain(TooShort) ,'Color',[0 0 1]  );  set( obj.ph_Domain(TooShort) ,'LineWidth',2  ) ;
            fprintf(' (red)Domains %s are too long for more than 1nm. \n',num2str(TooLong')) ;
            fprintf(' Total (red)Domains are %s. \n',num2str(length(TooLong))) ;
            
            fprintf(' (blue)Domains %s are too short for more than 1nm. \n',num2str(TooShort')) ;
            figure(obj.f_align) ;
            title(strcat('Number of domain =',num2str( length(Domains))  )  ) ;
        end
        function  SelectDomain(obj,src, evn)
            %             evn
            %     Domains=obj.MagicDNA_design.DuplexDomain;
            if evn.Button ==1
                set( obj.ph_Domain ,'LineStyle',':'  );
                set(src  ,'LineStyle','-' ) ;
                %                 set( obj.ph_Domain ,'Color',[0.2 0.2 0.2],'LineWidth',1  );
                %                 set(src ,'Color',[0.8 0.2 0.2]  ,'LineWidth',2) ;
                ind= src.UserData.ind  ;
                obj.Domain_info(ind,:) ;
                src.UserData.C5Base'
                fprintf('Domain ID: %i \n', src.UserData.ind ) ;
                fprintf('Domain Length : %4.2f bp and %4.2f nm \n',  obj.Domain_info(ind,:) ) ;
            elseif evn.Button==3
                CColor = get( obj.ph_Domain,'Color') ; CColor=cell2mat(CColor);
                IndBlue = and(and(CColor(:,1)==0 ,CColor(:,2)==0),CColor(:,3)==0.8)';
                LineStysss= get( obj.ph_Domain,'LineStyle') ;
                IndSolid =  (cell2mat(LineStysss) =='-')' ;
                PrecSelect = intersect(find( IndBlue),find(  IndSolid)) ;
                obj.ph_Domain(PrecSelect).UserData
                
                if ~isempty(PrecSelect) && sum(src.Color==[0 0 0.8])==3
                                    TrySomething =1
                    BB1 = obj.ph_Domain(PrecSelect).UserData.BaseBackBone ;
                    BB2 = src.UserData.BaseBackBone ;
                    mmDist1 = BB1-mean(BB2) ;  mmDist1=min(sqrt(sum(mmDist1.^2 ,2))) ;
                    mmDist2 = BB2-mean(BB1) ;  mmDist2=min(sqrt(sum(mmDist2.^2 ,2))) ;
                    
                    if mmDist1<3 && mmDist2<3
                        Node1 = 0.5*(BB1(1:end-1,:)+BB1(2:end,:) ) ;Node1(:,end+1)= (1:size(Node1,1))' ;
                        Node2 = 0.5*(BB2(1:end-1,:)+BB2(2:end,:) ) ;Node2(:,end+1)= (1:size(Node2,1))' ;
                        
                        Dist_all = zeros(size(Node1,1), size(Node2,1) ) ;
                        for j1 = 1:size(Node1,1)
                            for j2 =1:size(Node2,1)
                                Dist_all(j1,j2)= norm(Node1(j1,1:3)-Node2(j2,1:3)  );
                            end
                        end
                        AvailablComb =  Dist_all<2 ;
                        [u,v] = find(AvailablComb) ; GGind = find(AvailablComb) ;
                        UVd_list = [u,v,Dist_all(GGind)] ;
                        UVd_list = sortrows(UVd_list ,3) ; UVd_list=[UVd_list,zeros(size(UVd_list,1) ,1)] ;
                        UVd_list(1,4) =1 ;
                        for filt= 2:size(UVd_list,1)
                            PrevRow = UVd_list(1:filt-1 ,:) ;
                            ccr1 =PrevRow(:,4)==1 ;
                            ccr2 = abs(PrevRow(:,1)-ones(filt-1,1)*UVd_list(filt,1)) <3 ;
                            ccr3 = abs(PrevRow(:,2)-ones(filt-1,1)*UVd_list(filt,2)) <3 ;
                            Final = or(ccr2,ccr3) ;
                            Final(ccr1==0) =0;
                            if sum(Final )==0
                                UVd_list(filt,4)=1 ;
                            end
                        end
                        PossibleCase = UVd_list(UVd_list(:,4)==1 ,:) ;
                        XYZ=zeros(3*size(PossibleCase,1) ,3) ;
                        for k =1:size(PossibleCase,1)
                            XYZ(3*k-2:3*k,:) =[ Node1(PossibleCase(k,1),1:3) ; Node2(PossibleCase(k,2),1:3); nan nan nan ];
                            %
                        end
%                         plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3));
                        %                          TrySomething =1
                        
                        if isempty(  obj.ph_SuggestStapXover) || ~isvalid( obj.ph_SuggestStapXover)
                            obj.ph_SuggestStapXover = plot3(XYZ(:,1) ,XYZ(:,2) ,XYZ(:,3) ,'k' ,'LineWidth',3);
                            set(  obj.ph_SuggestStapXover, 'ButtonDownFcn', @(src,evn)SelectXover(obj,src, evn)) ;

                        else
                            obj.ph_SuggestStapXover.XData=XYZ(:,1)';
                            obj.ph_SuggestStapXover.YData=XYZ(:,2)';
                            obj.ph_SuggestStapXover.ZData=XYZ(:,3)';
                            set(  obj.ph_SuggestStapXover, 'ButtonDownFcn', @(src,evn)SelectXover(obj,src, evn)) ;
                        end
                        
                    end
                end
            else
                
                ind= src.UserData.ind  ;
                obj.Domain_info(ind,:) ;
                src.UserData.C5Base'
                fprintf('Domain ID: %i \n', src.UserData.ind ) ;
                fprintf('Domain Length : %4.2f bp and %4.2f nm \n',  obj.Domain_info(ind,:) ) ;
                
            end
            
            %             sdfsf=3
        end
        
        function SelectXover(obj,src, evn)
            XYZ = [src.XData; src.YData; src.ZData ]' ; 
            clickXYZ= evn.IntersectionPoint ;
            d= XYZ-clickXYZ ;
            dd = sum(d.^2 ,2) ;
            ind = find(dd==min(dd)) ;
            
           ind= ceil(ind/3); inds=[3*ind-2,3*ind-1];
            
                        if isempty(  obj.ph_SelectXover) || ~isvalid( obj.ph_SelectXover)
                         obj.ph_SelectXover=plot3(XYZ(inds,1),XYZ(inds,2),XYZ(inds,3), 'r','LineWidth',2  )  ;
                        else
                         obj.ph_SelectXover.XData=  XYZ(inds,1) ;
                         obj.ph_SelectXover.YData=  XYZ(inds,2) ;
                         obj.ph_SelectXover.ZData=  XYZ(inds,3) ;

                        end
        end
        
        
        function obj = DNACurve(MagicDNA_fH )
            obj.MagicDNA_fH =MagicDNA_fH ;
            ss_Assembly= findobj(MagicDNA_fH,'Tag','ss_Assembly') ;
            obj.MagicDNA_design= ss_Assembly.UserData.HyperBundle ;
            
            ss_Step= findobj(MagicDNA_fH,'Tag','ss_STEP') ;
            obj.Ori_freeformcurve= ss_Step.UserData.FFC ;
        end
        
        
        function visulize_initial(obj)
            %             axAssembly = findobj(obj.MagicDNA_fH, 'type','Axes','Tag','AssemblyMain') ;
            %             axes(axAssembly) ;
            %             %             SaveAxes_two ;
            %             CollectBundles_L = zeros(length(obj.MagicDNA_design.containBundle) ,1 ) ; %unit nm
            %             for k = 1: length(CollectBundles_L)
            %                 CollectBundles_L(k)=   abs(mean( obj.MagicDNA_design.containBundle{k}.Z1)- mean(obj.MagicDNA_design.containBundle{k}.Z2) ) ;
            %             end
            %
            %             Ori_Stline_L = zeros( length(CollectBundles_L)-1 ,1) ;
            %             for k=1: length(CollectBundles_L)
            %                 Ori_Stline_L(k) = norm(diff(  obj.Ori_freeformcurve.PointOnCurve{1}(k:k+1,:))) ;
            %             end
            %
            %             LineModelRatio = mean( CollectBundles_L./Ori_Stline_L ) + 0.3;
            %
            %             ScaledPosition = LineModelRatio * obj.Ori_freeformcurve.PointOnCurve{1} ;
            obj.Scaled_FFC = freeformcurve(2,obj.Ori_freeformcurve.PointOnScaleCurve_afterSlicing) ; %initail FFC object
            %----------
            %             obj.initial_oxDNA = oxDNA_FFCuse() ; % initialize object ;
            %             obj.freeformcurve.PointOnCurve{1}
            
            
        end
        
        function align_oxDNA_w_curve(obj)
            for notgood =1: 1
                obj.Scaled_FFC.VisualizeFFC_noGUI ;
                obj.f_align = obj.Scaled_FFC.fH2 ;
                
                figure(obj.f_align) ; hold on ;
                pH = obj.initial_oxDNA.visualTraj(1,1) ;
                % ----------
                %             MagicDNA_design.ScafAllBase
                BaseC5_inoxDNA = [-1,-1 ,-1,-1] ; % [C5 , base , strand , baseInStrand]
                strand =1 ;
                for k =1:length(  obj.MagicDNA_design.ScafAllBase )
                    Ez=obj.MagicDNA_design.ScafAllBase{k} ;
                    NewLines =  [Ez, strand*ones(size(Ez,1) ,1),[1:size(Ez,1)]' ] ;
                    BaseC5_inoxDNA= [ BaseC5_inoxDNA ;  NewLines];
                    strand=strand+1 ;
                end
                for k =1:length(  obj.MagicDNA_design.StapAllBase )
                    Ez=obj.MagicDNA_design.StapAllBase{k} ;
                    NewLines =  [Ez, strand*ones(size(Ez,1) ,1),[1:size(Ez,1)]' ] ;
                    BaseC5_inoxDNA= [ BaseC5_inoxDNA ;  NewLines];
                    strand=strand+1 ;
                end
                BaseC5_inoxDNA ;
                BaseC5_inoxDNA = BaseC5_inoxDNA(2:end ,:) ;
                Rtable =obj.MagicDNA_design.RelateTable ;
                %-----------
                for Bi = 1 : max( obj.initial_oxDNA.BM )
                    AllBase_ind =  obj.initial_oxDNA.BM  == Bi ;
                    
                    Target =  obj.Scaled_FFC.PointOnCurve{1}(Bi:Bi+1 ,:) ;
                    obj.MagicDNA_design
                    
                    %-----Z1 ;
                    Arrs = [-1, -1] ;
                    for Cyli = 1: length(obj.MagicDNA_design.containBundle{Bi}.Zbase1 )
                        C5=    obj.MagicDNA_design.RelateTable(and(Rtable(:,1)==Bi, Rtable(:,2)==Cyli)  ,5) ;
                        Arrs=[Arrs ; [C5 ,  obj.MagicDNA_design.containBundle{Bi}.Zbase1(Cyli)] ] ;
                    end
                    Arrs = Arrs(2:end ,:) ;
                    
                    [~,EndBasesZ1_indOX  ]= ismember( Arrs ,BaseC5_inoxDNA(:,1:2) ,'rows') ;
                    Z1side_XYZ = zeros(  size(Arrs,1) ,3);
                    for k=1:size(Z1side_XYZ ,1)
                        Z1side_XYZ(k,1)= pH{BaseC5_inoxDNA(EndBasesZ1_indOX(k),3 )}.XData(BaseC5_inoxDNA(EndBasesZ1_indOX(k),4 )) ;
                        Z1side_XYZ(k,2)= pH{BaseC5_inoxDNA(EndBasesZ1_indOX(k),3 )}.YData(BaseC5_inoxDNA(EndBasesZ1_indOX(k),4 )) ;
                        Z1side_XYZ(k,3)= pH{BaseC5_inoxDNA(EndBasesZ1_indOX(k),3 )}.ZData(BaseC5_inoxDNA(EndBasesZ1_indOX(k),4 )) ;
                    end
                    mm_Z1side= mean(Z1side_XYZ) ;
                    %----- Z2
                    Arrs = [-1, -1] ;
                    for Cyli = 1: length(obj.MagicDNA_design.containBundle{Bi}.Zbase2 )
                        C5=    obj.MagicDNA_design.RelateTable(and(Rtable(:,1)==Bi, Rtable(:,2)==Cyli)  ,5) ;
                        Arrs=[Arrs ; [C5 ,  obj.MagicDNA_design.containBundle{Bi}.Zbase2(Cyli)] ] ;
                    end
                    Arrs = Arrs(2:end ,:) ;
                    
                    [~,EndBasesZ2_indOX  ]= ismember( Arrs ,BaseC5_inoxDNA(:,1:2) ,'rows') ;
                    Z2side_XYZ = zeros(  size(Arrs,1) ,3);
                    for k=1:size(Z2side_XYZ ,1)
                        Z2side_XYZ(k,1)= pH{BaseC5_inoxDNA(EndBasesZ2_indOX(k),3 )}.XData(BaseC5_inoxDNA(EndBasesZ2_indOX(k),4 )) ;
                        Z2side_XYZ(k,2)= pH{BaseC5_inoxDNA(EndBasesZ2_indOX(k),3 )}.YData(BaseC5_inoxDNA(EndBasesZ2_indOX(k),4 )) ;
                        Z2side_XYZ(k,3)= pH{BaseC5_inoxDNA(EndBasesZ2_indOX(k),3 )}.ZData(BaseC5_inoxDNA(EndBasesZ2_indOX(k),4 )) ;
                    end
                    mm_Z2side= mean(Z2side_XYZ) ;
                    %------
                    Src_P = [mm_Z1side ; mm_Z2side]  ;
                    
                    [regParams,~,Acc2]=absor(Src_P',Target');
                    XYZ_Ori =    obj.collect_pH(BaseC5_inoxDNA(AllBase_ind,3:4) , pH) ;
                    
                    %                Cyls_C5 = obj.MagicDNA_design.RelateTable(obj.MagicDNA_design.RelateTable(:,1)==Bi ,5) ;
                    NewXYZ =  regParams.R* XYZ_Ori'+ regParams.t ;
                    NewXYZ=NewXYZ' ;
                    obj.Deploy_pH( BaseC5_inoxDNA(AllBase_ind,3:4) , pH,NewXYZ) ;
                    
                end
                obj.BaseC5inoxDNA=BaseC5_inoxDNA ;
                set(obj.f_align,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) )  ;
            end
        end
        
        function extrudePerfectModel_V2(obj)
            obj.Scaled_FFC.VisualizeFFC_noGUI ;
            obj.f_align = obj.Scaled_FFC.fH2 ;
            Rtable =obj.MagicDNA_design.RelateTable ;
            figure(obj.f_align) ; hold on ; title(' This is vector field.')
            %---------  MagicDNA_design.ScafAllBase
            BaseC5_inoxDNA = [-1,-1 ,-1,-1] ; % [C5 , base , strand , baseInStrand]
            strand =1 ;
            for k =1:length(  obj.MagicDNA_design.ScafAllBase )
                Ez=obj.MagicDNA_design.ScafAllBase{k} ;
                NewLines =  [Ez, strand*ones(size(Ez,1) ,1),[1:size(Ez,1)]' ] ;
                BaseC5_inoxDNA= [ BaseC5_inoxDNA ;  NewLines];
                strand=strand+1 ;
            end
            ScafInd = 1: size(BaseC5_inoxDNA,1)-1 ;
            for k =1:length(  obj.MagicDNA_design.StapAllBase )
                Ez=obj.MagicDNA_design.StapAllBase{k} ;
                NewLines =  [Ez, strand*ones(size(Ez,1) ,1),[1:size(Ez,1)]' ] ;
                BaseC5_inoxDNA= [ BaseC5_inoxDNA ;  NewLines];
                strand=strand+1 ;
            end
            BaseC5_inoxDNA = BaseC5_inoxDNA(2:end ,:) ;
            obj.BaseC5_inoxDNA_sBase =BaseC5_inoxDNA ; % single base C5 and [strand base] in ini. oxDNA

            %-----------------
            SpineCurve = obj.Ori_freeformcurve.Scaled_Curve ;
            
            %----allocate larger array,
            All_C5_Cyl_traj = zeros(round(1.5* size(BaseC5_inoxDNA,1)) , 11 ) ;   %[C5 ,CylPos, ds_TVec, ds_NVec ]
            %             All_ek_OnTraj= zeros(round(1.5* size(BaseC5_inoxDNA,1)) , 7 ) ;
            ct_b = 1 ;
            All_Base_onCyl_Pos= obj.Ori_freeformcurve.Cylinderstt.All_Base_onCyl_Pos ;
            
            %             BaseCC_OnCyl = zeros(obj.Ori_freeformcurve.Cylinderstt.N,1 ) ;
            for Bi = 1:  length(obj.MagicDNA_design.containBundle)
                %                 Bi
                Bundle = obj.MagicDNA_design.containBundle{Bi} ;
                InplaneXY = (Bundle.CylInplanePosition-mean(Bundle.CylInplanePosition))*2.5/2 ; % spacing to 2.5nm
                
                NeutralCurve_t = obj.Ori_freeformcurve.Bundlett.New_tarr(Bi:Bi+1) ;
                cc= 1 ;
                NBaseOnCylinders = zeros( size(InplaneXY,1) ,1) ;
                for cyl_k = 1: size(InplaneXY,1)
                    [~,iid] = ismember([Bi cyl_k], Rtable(:,1:2), 'rows'  ) ;
                    Cyl_C5 = Rtable(iid,5) ;
                    RefOnJoint = obj.Ori_freeformcurve.Cylinderstt.DotOnCylinders(cyl_k,Bi:Bi+1,:) ;
                    RefOnJoint= reshape(RefOnJoint,2,3)  ;
                    
                    Query_arr =  Bundle.Zbase1(cyl_k): Bundle.Zbase2(cyl_k) ;
                    if ~isempty( Bundle.skipPosition)
                    Query_arr= setdiff(Query_arr  , Bundle.skipPosition(Bundle.skipPosition(:,2)==cyl_k,4 )) ;
                    end
                    if ~isempty( Bundle.InsertPosition)
                    Query_arr=sort([Query_arr,Bundle.InsertPosition(Bundle.InsertPosition(:,2)==cyl_k,4 )']);
                    end
                                        
                    
                    [~,IndOnCcurve]= ismember(RefOnJoint , All_Base_onCyl_Pos{cyl_k},'rows') ;
                    RefOnJoint2=All_Base_onCyl_Pos{cyl_k}(IndOnCcurve(1):IndOnCcurve(2) ,:) ;
                    
                    [~,IndOnCcurve2]= ismember(RefOnJoint2 , All_Base_onCyl_Pos{cyl_k},'rows') ;
                    
                    QueryBases= Query_arr-Bundle.Zbase1(cyl_k)+IndOnCcurve2(1) ;
                    QueryBases=linspace(QueryBases(1)+0.5,QueryBases(end)-0.5, length(QueryBases)) ;
                    % avoid repeat base at boundary, and interpolate again.
                    BasePos = interp1( IndOnCcurve2, RefOnJoint2 , QueryBases ,'pchip','extrap')  ;
                    
                    %                     if sum(sum(isnan(BasePos)))>0
                    %                         sdfs=3
                    %                     end
                    
                    %---------------------------- TVec
                    t_base_arr = linspace(NeutralCurve_t(1), NeutralCurve_t(2), length(Query_arr) ) ;
                    Tvec = fnval(fnder(SpineCurve), t_base_arr) ;
                    Tvec=vec3norm_CM(Tvec') ;
                    %---------
                    
%                     NVec = interp1( NeutralCurve_t,obj.Ori_freeformcurve.Bundlett.NVec, t_base_arr  )
                    
                    %                     vq = interp1(x,v,xq)
%                     sdfsdf=3
                    %----------
                    
                    
                    All_C5_Cyl_traj(ct_b:ct_b+ length(Query_arr)-1 ,1:2) =[Cyl_C5*ones(1,length(Query_arr))  ; Query_arr]' ;
                    All_C5_Cyl_traj(ct_b:ct_b+ length(Query_arr)-1 ,3:5) =BasePos ;
                    All_C5_Cyl_traj(ct_b:ct_b+ length(Query_arr)-1 ,6:8) =Tvec ;
                    
                    ct_b=ct_b+length(Query_arr) ;
                    %----------
 
                end
            end
            All_C5_Cyl_traj=All_C5_Cyl_traj(1:ct_b-1 ,:);
            %       scatter3(All_C5_Cyl_traj(:,3) ,All_C5_Cyl_traj(:,4), All_C5_Cyl_traj(:,5),'.' ) ;
            
            
            Inds = 1:10:ct_b-1 ;
            Ez =All_C5_Cyl_traj(:,3:8)' ;
            % q1 = quiver3(Ez(1,Inds),Ez(2,Inds),Ez(3,Inds),Ez(4,Inds),Ez(5,Inds),Ez(6,Inds) ) ;
            
            %       quiver3(All_C5_Cyl_traj(:,3) ,All_C5_Cyl_traj(:,4), All_C5_Cyl_traj(:,5),'.' ) ;
            
            %-----------
            %    BaseC5_inoxDNA  ;
            %             %------return result
            obj.MatrixForObjFunc =  All_C5_Cyl_traj ; % double-stranded C5 and quivers,
            %             set(gcf,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) ) ;
            xlabel('X') ; ylabel('Y') ; zlabel('Z');
            title(' This is vector field.')            ;axis equal ;
            
            
        end
        
        
        function extrudePerfectModel(obj)
            obj.Scaled_FFC.VisualizeFFC_noGUI ;
            obj.f_align = obj.Scaled_FFC.fH2 ;
            Rtable =obj.MagicDNA_design.RelateTable ;
            figure(obj.f_align) ; hold on ; title(' This is vector field.')
            %---------  MagicDNA_design.ScafAllBase
            BaseC5_inoxDNA = [-1,-1 ,-1,-1] ; % [C5 , base , strand , baseInStrand]
            strand =1 ;
            for k =1:length(  obj.MagicDNA_design.ScafAllBase )
                Ez=obj.MagicDNA_design.ScafAllBase{k} ;
                NewLines =  [Ez, strand*ones(size(Ez,1) ,1),[1:size(Ez,1)]' ] ;
                BaseC5_inoxDNA= [ BaseC5_inoxDNA ;  NewLines];
                strand=strand+1 ;
            end
            ScafInd = 1: size(BaseC5_inoxDNA,1)-1 ;
            for k =1:length(  obj.MagicDNA_design.StapAllBase )
                Ez=obj.MagicDNA_design.StapAllBase{k} ;
                NewLines =  [Ez, strand*ones(size(Ez,1) ,1),[1:size(Ez,1)]' ] ;
                BaseC5_inoxDNA= [ BaseC5_inoxDNA ;  NewLines];
                strand=strand+1 ;
            end
            BaseC5_inoxDNA ;
            BaseC5_inoxDNA = BaseC5_inoxDNA(2:end ,:) ;
            %-----------------
            SpineCurve = obj.Scaled_FFC.SplineCurve{1} ;
            
            AllBundle_Mid_ei = zeros(3, length(obj.MagicDNA_design.containBundle) ) ;
            for k=1 :  size(AllBundle_Mid_ei ,2)
                AllBundle_Mid_ei(:,k) =   diff(obj.MagicDNA_design.containBundle{k}.LocalCoordinatFromLineModel(:,1:2) ,1,2) ;
                AllBundle_Mid_ei(:,k)=AllBundle_Mid_ei(:,k)/norm(AllBundle_Mid_ei(:,k)) ;
            end
            
            %----allocate larger array,
            All_C5_Cyl_traj = zeros(round(1.5* size(BaseC5_inoxDNA,1)) , 8 ) ;
            %             All_ek_OnTraj= zeros(round(1.5* size(BaseC5_inoxDNA,1)) , 7 ) ;
            ct_b = 1 ;
            
            for Bi = 1:  length(obj.MagicDNA_design.containBundle)
                Bundle = obj.MagicDNA_design.containBundle{Bi} ;
                InplaneXY = (Bundle.CylInplanePosition-mean(Bundle.CylInplanePosition))*2.5/2 ; % spacing to 2.5nm
                
                StartP =  fnval( SpineCurve , SpineCurve.breaks(Bi) ) ;
                EndP =   fnval( SpineCurve , SpineCurve.breaks(Bi+1) ) ;
                %                 e_k_Mid =  fnval(fnder(SpineCurve), mean(SpineCurve.breaks(Bi:Bi+1))) ; e_k_Mid=e_k_Mid/norm(e_k_Mid) ;
                
                StartEnd_all = [ SpineCurve.breaks(Bi:Bi+1)' , [ round(mean(Bundle.Zbase1)); round(mean(Bundle.Zbase2))] ]  ; %neutral curve
                %                 [t ,Base ]
                N_slice  =length(StartEnd_all(1,2):StartEnd_all(2,2)) ;
                Cyl_traj = zeros(3, (N_slice+1)* size(InplaneXY,1) ) ; SaveC5Rep = zeros(2, (N_slice+1)* size(InplaneXY,1) ) ;
                ek_OnTraj = zeros(3, (N_slice+1)* size(InplaneXY,1) ) ;
                cc= 1 ;
                NBaseOnCylinders = zeros( size(InplaneXY,1) ,1) ;
                for cyl_k = 1: 1 %size(InplaneXY,1)
                    [~,iid] = ismember([Bi cyl_k], Rtable(:,1:2), 'rows'  ) ;
                    Cyl_C5 = Rtable(iid,5) ;
                    
                    for bs =  Bundle.Zbase1(cyl_k): Bundle.Zbase2(cyl_k)   %              StartEnd_all(1,2):StartEnd_all(2,2)
                        t_Atslice =interp1([ Bundle.Zbase1(cyl_k)-1; Bundle.Zbase2(cyl_k)+1] , StartEnd_all(:,1)  , bs ) ;
                        P_Atslice =  fnval(SpineCurve, t_Atslice) ;
                        
                        e_k =  fnval(fnder(SpineCurve), t_Atslice)  ; e_k=e_k/norm(e_k) ;
                        e_i =fnval(fnder(fnder(SpineCurve)), t_Atslice) ;
                        e_i=e_i/norm(e_i)  ;  %e_i=e_i';
                        
                        e_j = cross(e_k,e_i) ; e_i= cross(e_j ,  e_k) ;
                        
                        Cyl_traj(: ,cc) = P_Atslice+ InplaneXY(cyl_k,1)*e_i+ InplaneXY(cyl_k,2)*e_j ;
                        ek_OnTraj(:,cc) = e_k ;
                        SaveC5Rep(:,cc) = [Cyl_C5 ; bs];
                        cc=cc+1 ; NBaseOnCylinders(cyl_k)=Bundle.Zbase2(cyl_k) -Bundle.Zbase1(cyl_k) +1 ;
                    end
                    Cyl_traj(: ,cc) =[nan;nan;nan] ;
                    ek_OnTraj(: ,cc) =[nan;nan;nan] ;
                    SaveC5Rep(: ,cc) =[nan;nan] ;
                    cc=cc+1 ;
                end
                Cyl_traj=Cyl_traj(:, 1:cc-1) ;
                ek_OnTraj=ek_OnTraj(:, 1:cc-1) ;
                SaveC5Rep=SaveC5Rep(:, 1:cc-1) ;
                plot3(Cyl_traj(1,:),Cyl_traj(2,:),Cyl_traj(3,:),'-') ;
                %               Inds = ~isnan(Cyl_traj(1,:)) ;
                %                                 quiver3(Cyl_traj(1,:),Cyl_traj(2,:),Cyl_traj(3,:) ,ek_OnTraj(1,:),ek_OnTraj(2,:),ek_OnTraj(3,:))
                %                 quiver3(x,y,z,u,v,w)
                Merge = [SaveC5Rep ; Cyl_traj; ek_OnTraj ] ;
                Inds = ~isnan(Cyl_traj(1,:)) ;
                Merge=Merge(:,Inds)  ; nn= size(Merge,2) ;
                
                All_C5_Cyl_traj(ct_b:ct_b+nn-1 ,:) =Merge' ;
                ct_b=ct_b+nn ;
            end
            All_C5_Cyl_traj=All_C5_Cyl_traj(1:ct_b-1,:) ;
            BaseC5_inoxDNA ;
            %-----------
            %--------relax cylinders of multiple bundles.
            %             fprintf('Warning!! Only works for constant cross-sections!!\n')
            %             for Cxi =1: size(InplaneXY,1)
            %                 [~, All_Cyl_inds]= ismember( Rtable(:,2) ,Cxi)  ;
            %                 C5cylinders =  Rtable(  find( All_Cyl_inds),5)  ;
            %                 [~,   AllbaseInd_in_table ] =ismember( All_C5_Cyl_traj(:,1) ,C5cylinders ) ;
            %                 Inds = find(AllbaseInd_in_table) ;
            %
            %                 Sub_alongCyls =  All_C5_Cyl_traj(Inds ,:) ;
            %                 Ds = diff(Sub_alongCyls(:,3:5) ) ;
            %                 Base_Base_Dist = sqrt(sum(Ds.^2 ,2) );
            %                 X_before = cumsum([0;Base_Base_Dist] ) ; X_before=X_before/max(X_before) ;
            %                 %              figure;plot(Base_Base_Dist) ;
            %                 V_after = interp1(X_before, Sub_alongCyls(:,3:5) , [1:size(Sub_alongCyls,1)]/size(Sub_alongCyls,1) ) ;
            %
            %                 All_C5_Cyl_traj(Inds ,3:5) = V_after ;
            %                 %              Ds2 = diff(V_after)
            %                 %              Base_Base_Dist2 = sqrt(sum(Ds2.^2 ,2) );
            %                 %              figure;plot(Base_Base_Dist2) ;
            %
            %             end
            %
            
            %------return result
            obj.BaseC5_inoxDNA_sBase =BaseC5_inoxDNA ; % single base C5 and [strand base] in ini. oxDNA
            obj.MatrixForObjFunc =  All_C5_Cyl_traj ; % double-stranded C5 and quivers,
            %             set(gcf,'KeyPressFcn',@(src,evn)ChangeAxesLimit(src,evn) ) ;
            xlabel('X') ; ylabel('Y') ; zlabel('Z');
            title(' This is vector field.')            ;axis equal ;
            
        end
        
        function XYZ=collect_pH(obj, Strand_Base, pH)
            XYZ= zeros(size(Strand_Base,1) ,3) ;
            for k=1: size(XYZ,1)
                %                 if Strand_Base(k,2) ==39
                %                     sdfsf=3
                %                 end
                XYZ(k,1) =   pH{Strand_Base(k,1)}.XData(Strand_Base(k,2)) ;
                XYZ(k,2) =   pH{Strand_Base(k,1)}.YData(Strand_Base(k,2)) ;
                XYZ(k,3) =   pH{Strand_Base(k,1)}.ZData(Strand_Base(k,2)) ;
            end
        end
        
        function Deploy_pH(obj, Strand_Base, pH,NewXYZ)
            for k=1: size(NewXYZ,1)
                pH{Strand_Base(k,1)}.XData(Strand_Base(k,2))=  NewXYZ(k,1);
                pH{Strand_Base(k,1)}.YData(Strand_Base(k,2))=  NewXYZ(k,2);
                pH{Strand_Base(k,1)}.ZData(Strand_Base(k,2))=  NewXYZ(k,3);
            end
        end
        
        
    end
    
end

