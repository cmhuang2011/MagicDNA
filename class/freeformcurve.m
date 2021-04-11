classdef freeformcurve <handle
    % This class handles the free-form sketch tool in the very beginning of
    % MagicDNA. Users can specify spline/straight-line objects to
    % MagicDNA?through this GUI.
    %   Detailed explanation goes here
    
    properties
        % -------cell array
        PointOnCurve  ;  % blue dot as point on the curve
        SplineCurve  ; %  curve object passing PointOnCurve as free form curve
        CurveParas;   % red curve as xyz and t
        
        Cline; % line object of CurveParas
        Cscatter ; % scatter object of PointOnCurve
        %-----------
        TempAddPoint  ;
        ls=1  ;% value of popupH2
        %--- below coming from GUI_continuouSlicing
        PointOnScaleCurve = [];  %joints between bundles in straight line model.
        Scaled_Curve =[];
        PointOnScaleCurve_afterSlicing = [];  %joints between bundles.
        t_arrs_fine = [];   %resolution is about half of bps
        PTN_Vec_onCurve = []; 
        Bundlett =[];
        Cylinderstt =[];
        
    end
    
    properties (Hidden = true)
        popupH ; 
        TRtxt;%translate or rotate interval
        popupH2 ; % select one of lines
        check_keyControl ;
        check_SWsplineStraight ;
        UIgroup_XYZ_AlignProj;
        btn1_projXY ;  btn2_projXZ ;   btn3_projYZ ;
        
        UIgroup_Other;
        btn4_Save ;  btn5_Load ;
        btn6_Export ; btn7_Add;
        
        fH ; fH2 ;
        ss_STEP ;  % tab handle to MagicDNA
        editH;
        XYZlineOnPoint=struct('x',[],'y',[],'z',[] ) ; % xyz local
    end
    
    methods
        function  exCommand(obj)
            FFC  = freeformcurve(2, xyz')   ;
            fnval( obj.Scaled_FFC.SplineCurve{1} , 4.9405) ;
            
            
            ss_STEP= findobj(0,'Tag','ss_STEP') ;
            FFC = ss_STEP.UserData.FFC ;
            FFC.get_Equi_Distance_Point(1,50) ;
        end
        
        
        function get_Equi_Distance_Point(obj , CurveInd , N_point)            
            spCurve = obj.SplineCurve{CurveInd} ;
            fnprime = fnder(spCurve);    Lfun = @(s) sqrt(sum(fnval(fnprime,s).^2,1));
            s_array = linspace(min(spCurve.breaks) ,max(spCurve.breaks) , 1000) ;
            Ori_XYZ_fine =fnval(spCurve , s_array) ;
            x= Ori_XYZ_fine(1,:) ;   y= Ori_XYZ_fine(2,:) ;   z= Ori_XYZ_fine(3,:) ;
            
            ds = zeros(1, length(s_array )-1 ) ;
            for k =1: length(ds)
                ds(k) =    integral(Lfun, s_array(k), s_array(k+1)) ;
            end
            Cs =  cumsum(ds) ;
            TargetS = linspace(0, max(Cs) , N_point) ;
            IndWeightL =  zeros(size(TargetS) ) ;
            for k =1:length(IndWeightL)
                IndWeightL(k) = find(  min(abs(TargetS(k) -Cs))==abs(TargetS(k) -Cs)) ;
            end
            x=x(IndWeightL) ; y =y(IndWeightL) ; z =z(IndWeightL) ;
            obj.PointOnCurve{CurveInd} =[x;y;z]' ;
            
        end
        
        function obj = freeformcurve(type,varargin)
            if type ==1 % No input sketch from blank
                LoadXYZ(obj,[],[],[]) ;  % load Hi.xyz
                if nargin>2
                    obj.ss_STEP =   varargin{2} ;
                end
            elseif type ==2 % xyz input
                obj.PointOnCurve{obj.ls} = varargin{obj.ls} ;
                obj.SplineCurve{obj.ls} =  cscvn( obj.PointOnCurve{obj.ls}(1:end,: )') ;
                nargin;
                if nargin>2
                    obj.ss_STEP =   varargin{2} ;
                end
                
            end
        end
        function VisualizeFFC_noGUI(obj)
            obj.fH2 = figure(2425) ; clf ;
            %             ax= gca ; ax.Position = [0.05 0.05 0.75 0.9 ] ;
            for k = 1 :length(obj.PointOnCurve)
                obj.Cscatter{k} = scatter3(obj.PointOnCurve{k}(:,1),obj.PointOnCurve{k}(:,2) ,obj.PointOnCurve{k}(:,3),58,'filled'  );
                obj.Cscatter{k}.CData= repmat( [0 0.45 0.75 ],length(obj.PointOnCurve{k}(:,1)),1  );
                %             obj.Cscatter.ButtonDownFcn
                
                hold on ;
                [obj.CurveParas.p{k},obj.CurveParas.t{k}] = fnplt(obj.SplineCurve{k}) ;
                obj.Cline{k}=  plot3( obj.CurveParas.p{k}(1,:),obj.CurveParas.p{k}(2,:),obj.CurveParas.p{k}(3,:),'LineWidth',3 );
                obj.Cline{k}.HitTest='off' ;
            end
            
            %---------export spline to Chimera
            %             prompt = {'Enter a value of scale '};
            %             dlgtitle = 'Scale Value';
            %             definput = {'1'};
            %             dims = [1 40];
            %             opts.Interpreter = 'tex';
            %             answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
            %             Rs = str2double(answer{1});
            
            Rs=1 ;
            Thicker = 0.3;
            file_name='Spline_extrusion'
            Radius =0.15*Thicker;
            fileID = fopen([pwd filesep file_name '.bild'],'w');
            for Bi = 1:length( obj.Cline )
                fprintf(fileID ,'\n' );
                %                         fprintf(fileID ,'.transparency 0.5\n' );
                fprintf(fileID ,'.comment (If need colors for bundle %i) .color  r g b\n',Bi );
                s_array = linspace( obj.SplineCurve{Bi}.breaks(1),obj.SplineCurve{Bi}.breaks(end) , 1000 ) ;
                XYZ_onCyl =    fnval(obj.SplineCurve{Bi} , s_array) ;
                XYZ_onCyl=XYZ_onCyl';
                
                for Cj= 1:size(XYZ_onCyl,1)-1
                    fprintf(fileID , '.cylinder %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f open\n',XYZ_onCyl(Cj,:),XYZ_onCyl(Cj+1,:),Radius )    ;
                    fprintf(fileID , '.sphere %8.6f %8.6f %8.6f  %8.6f \n',XYZ_onCyl(Cj,:),Radius )    ;
                end
                fprintf(fileID , '.sphere %8.6f %8.6f %8.6f  %8.6f \n',XYZ_onCyl(Cj+1,:),Radius )    ;
            end
            fclose(fileID);
            %--------------
            file_name='Spline_BreakPoint' 
            filename2 = 'Spline_BPconnect'
            Radius =0.3*Thicker;
            fileID = fopen([pwd filesep file_name '.bild'],'w');
            fileID2 = fopen([pwd filesep filename2 '.bild'],'w');
            for Bi = 1:length( obj.Cline )
                fprintf(fileID ,'\n' );
                %                         fprintf(fileID ,'.transparency 0.5\n' );
                fprintf(fileID ,'.comment (If need colors for bundle %i) .color  r g b\n',Bi );
                
                fprintf(fileID2 ,'.comment (If need colors for bundle %i) .color  r g b\n',Bi );
                fprintf(fileID2 ,'.color %4.2f %4.2f %4.2f \n' , rand(1,3)) ;
                XYZ_onCyl = Rs*[obj.Cscatter{Bi}.XData ;obj.Cscatter{Bi}.YData;obj.Cscatter{Bi}.ZData]' ;
                
%                 s_array = linspace( obj.SplineCurve{1}.breaks(1),obj.SplineCurve{1}.breaks(end) , 1000 ) ;
%                 XYZ_onCyl =    fnval(obj.SplineCurve{Bi} , s_array) ;
                
                for Cj= 1:size(XYZ_onCyl,1)
                    fprintf(fileID , '.sphere %8.6f %8.6f %8.6f  %8.6f \n',XYZ_onCyl(Cj,:),Radius )    ;
                    if Cj~=1
                       fprintf(fileID2 , '.cylinder %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',XYZ_onCyl(Cj-1,:),XYZ_onCyl(Cj,:),Radius*0.7 )    ; 
                    end                    
                end
            end
            fclose(fileID);            
            fclose(fileID2);      
        end
        
        function Cont_FFC_data=VisualizeFFC_noGUI_wNormalVector(obj,InplaneXY,VecAll, L_bundle,XsecLattice)
            %Don't Use.
            % This function was origainally for sweeping. However, the
            % codes are implemented in other external function.
            fprintf('Start Wnormal \n')
            f2426=figure(2426) ; cla ;   hold on ;
            nCylinder =size(InplaneXY ,1) ;  nBundle =size(obj.PointOnCurve{1},1)-1 ;
            
            Bundle1Length_in_nm = L_bundle(1) *0.34 ;
            Ori_L_curve = norm( diff(obj.PointOnCurve{1}(1:2,:)) ) ;
            
            obj.PointOnScaleCurve = obj.PointOnCurve{1}*Bundle1Length_in_nm/Ori_L_curve ; POSC =obj.PointOnScaleCurve;
            obj.Scaled_Curve =    cscvn(  obj.PointOnScaleCurve') ; Scale_C =obj.Scaled_Curve ;
            %             ax= gca ; ax.Position = [0.05 0.05 0.75 0.9 ] ;
            for k = 1 :1 %length(obj.PointOnCurve)
                AA = scatter3(POSC(:,1),POSC(:,2) ,POSC(:,3),58,'filled'  );
                AA.CData= repmat( [0 0.45 0.75 ],length(POSC(:,1)),1  );
                [WW.p,obj.CurveParas.t] = fnplt(Scale_C) ;
                BB=  plot3( WW.p(1,:),WW.p(2,:),WW.p(3,:),'LineWidth',3 );
                BB.HitTest='off' ;
            end
            ax=gca;
            %---------
            N_slice= 2*nBundle+1 ; % 2n+1 , or N_slice =1000
            InitOrt_X = VecAll(1,7:9)';  InitOrt_X=InitOrt_X/norm(InitOrt_X) ;  % initial orientation of x from sketch
            for k = 1 :1
                t_arrs= sort( [Scale_C.breaks ,movmean(Scale_C.breaks , 2,'Endpoints','discard')] )  ;
                Position_TanVec_NormalVec = zeros(9,N_slice ) ;
                for  pointi = 1: N_slice
                    Position_TanVec_NormalVec(1:3, pointi) = fnval(Scale_C, t_arrs(pointi))  ;% position
                end
            end
            
            QQ= gradient(Position_TanVec_NormalVec(1:3,: ) ) ; % tangent direction
            Position_TanVec_NormalVec(4:6,:)=  vec3norm_CM( QQ' )' ;
            Position_TanVec_NormalVec(7:9,1) =InitOrt_X ;
            %------Make continuous
            %             Formers =QQ2(:,1:end-1 ); Laters =QQ2(:,2:end );
            for k =2:    N_slice
                NVec =  Position_TanVec_NormalVec(7:9,k-1) ;
                CurTVec = Position_TanVec_NormalVec(4:6,k) ;
                PrevTVec = Position_TanVec_NormalVec(4:6,k-1) ;
                if dot(CurTVec,PrevTVec)>0  % angle smaller than 90
                    Ytemp= cross(CurTVec, NVec) ;Ytemp=Ytemp/norm(Ytemp) ;
                else
                    Ytemp= cross(-CurTVec, NVec) ;Ytemp=Ytemp/norm(Ytemp)  ;
                end
                New_NVec = cross(Ytemp,Position_TanVec_NormalVec(4:6,k)) ;
                Position_TanVec_NormalVec(7:9,k) = New_NVec ;
            end
            %------------local orientation of each slice
            InplaneXY=InplaneXY-mean(InplaneXY) ;
            InplaneXY=InplaneXY*1.25 ;
            
            S_PTNC =  Position_TanVec_NormalVec ;
            N_FineSlice = round(2* sum(L_bundle) );  % Don't need to be too dense
            Fine_tarr = linspace(t_arrs(1),t_arrs(end), N_FineSlice) ;
            PTN_vec_onSpline = interp1(t_arrs,S_PTNC', Fine_tarr  ,'PCHIP'  ) ;
            
            SumFFC_unitBPS = sum(L_bundle) ;
            
            PTN_vec_onSpline(:,1:3) = fnval(Scale_C, Fine_tarr)'  ;% position
            
            QQ = fnval( fnder(Scale_C), Fine_tarr)'  ; % tangent direction
            PTN_vec_onSpline(:,4:6)= vec3norm_CM( QQ ) ;   %Tan Vec
            
            Temp= cross(  PTN_vec_onSpline(:,7:9 ), PTN_vec_onSpline(:,4:6) ) ;
            QQ2 =cross(PTN_vec_onSpline(:,4:6 ), Temp ) ;
            PTN_vec_onSpline(:,7:9)=  vec3norm_CM( QQ2 ) ;
            
            PTN_vec_onSpline(: , 10:12) =cross(   PTN_vec_onSpline(: , 4:6) ,   PTN_vec_onSpline(: , 7:9)) ;
            % [N_FineSlice x 12 ], [ P ; Tang ; Nor(x) ; y ]
            
            AllXYZ = zeros(3*nCylinder , N_FineSlice ) ;
            for clyi=1:nCylinder
                LocalXY = InplaneXY(clyi ,:) ;
                P_onCyl = PTN_vec_onSpline(:,1:3) + LocalXY(1)*PTN_vec_onSpline(:,7:9) ++ LocalXY(2)*PTN_vec_onSpline(:,10:12) ;
                AllXYZ(3*clyi-2:3*clyi ,:) =P_onCyl' ;
            end
            
            Cyl_h = gobjects(nCylinder,1) ;
            for k =  1: nCylinder
                Cyl_h(k) =plot3(AllXYZ(3*k-2 ,:),AllXYZ(3*k-1 ,:),AllXYZ(3*k ,:)  ,'-') ;
                Cyl_h(k).ButtonDownFcn = @(src,evn)ShowInd(obj,src,evn) ;
                Cyl_h(k).HitTest='off' ;
            end
            %----------
            set(Cyl_h,'Visible','off') ;
            set(Cyl_h(end),'Visible','on') ;
            set(Cyl_h(1),'Visible','on') ;
            %--------Visualize
            
            simplifyLine_byGtheory(obj,Cyl_h) ;
            
            Inds = 1:2:N_FineSlice ;   Inds = 1:1:N_slice ;
            Ez =PTN_vec_onSpline' ;     Ez =Position_TanVec_NormalVec ;
            
            q1 = quiver3(Ez(1,Inds),Ez(2,Inds),Ez(3,Inds),Ez(4,Inds),Ez(5,Inds),Ez(6,Inds) ) ;
            q2= quiver3(Ez(1,Inds),Ez(2,Inds),Ez(3,Inds),Ez(7,Inds),Ez(8,Inds),Ez(9,Inds) ) ;
            q1.HitTest='off' ;    q2.HitTest='off' ;
            
            axis equal ;
            xlabel('X'); ylabel('Y') ;zlabel('Z');
            set(figure(2426),'KeyPressFcn',@(src,evn)ChangeAxesLimit(obj,src,evn) )  ;
            
            
            fprintf('finish Wnormal \n')
            
            %             return
            %----------------decide how many bundles on this curve
            Rough_split_ForEachSlice = 20 ;  %
            sH=scatter3(AllXYZ(1,1),AllXYZ(1,2) ,AllXYZ(1,3),64 ,'filled'  );    sH.CData=[0 1 1];%initialize and some point
            sH.CData=[0 1 1]; sH.UserData.tarrs= t_arrs ; sH.UserData.Scale_C=Scale_C ;
            sH.ButtonDownFcn= @(src,evn)ClickSH(obj,src,evn) ;
            
            
            %              for k=1:5
            %             while 1
            [New_tarr, NewNodes]=SlicingCoarseCut(obj ,Rough_split_ForEachSlice,Scale_C ,N_slice ,t_arrs) ;
            sH.XData= NewNodes(:,1)' ; sH.YData= NewNodes(:,2)' ; sH.ZData= NewNodes(:,3)' ;  % update graphic
            
            %                 Result = AskGoodOrNot(Rough_split_ForEachSlice)   ;
            %                 Rough_split_ForEachSlice= Result.N ;
            %
            %                 switch Result.Info
            %                     case 'Close'
            %                         break
            %                     case 'Assign'
            %                         Rough_split_ForEachSlice=Result.N ;
            %                 end
            % %             end
            %----------change nodes on neutral ;
            
            %             pause()
            %---------------
            [DotOnCylinders,ds_Cyli_ForBundle ]=SlicingCoarseCut_part2(obj,nCylinder,NewNodes,Cyl_h) ;
            
            %--------
            AllJoints = reshape(DotOnCylinders,size(DotOnCylinders,1)*size(DotOnCylinders,2),3) ;
            scatter3(ax,AllJoints(:,1),AllJoints(:,2),AllJoints(:,3),'r.') ;
            
            
            %--------------------------------------------------------------
            %            % Exprot Hyperbundle ------------------
            %           Cont_FFC_data = [] ;
            %           Cont_FFC_data
            Center_t =  0.5*(New_tarr(1:end-1) +New_tarr(2:end));
            %          Bundle.Center =     fnval(Scale_C , Center_t)' ;
            Bundle.Center =   0.5*( NewNodes(1:end-1,:)+ NewNodes(2:end,:) ) ;
            %          NewNodes
            
            %          Bundle.TanVec =     fnval(fnprime, Center_t )'  ;
            Bundle.TanVec =NewNodes(2:end,:)-NewNodes(1:end-1,:) ;
            Bundle.NVec= interp1(Fine_tarr,PTN_vec_onSpline(:,7:9), Center_t  ,'PCHIP'  ) ;
            Cont_FFC_data.Bundle=Bundle ;
            Cylinders.N = nCylinder ;
            Cylinders.Nbundles = size(NewNodes,1)-1 ;
            Cylinders.DotOnCylinders=DotOnCylinders ;
            Cylinders.ds_Cyli_ForBundle=ds_Cyli_ForBundle;
            
            Cont_FFC_data.Bundle=Bundle ;
            Cont_FFC_data.Cylinders=Cylinders ;
            
            figure(325);clf;
            subplot(3,1,1) ; hold on ;
            NBaseOncyl= zeros(nCylinder,1) ;
            for k =  1: nCylinder
                XYZ = [ Cyl_h(k).XData ; Cyl_h(k).YData ; Cyl_h(k).ZData ]' ;
                dx = diff(XYZ) ;  qq = sqrt(sum(dx.^2 ,2)) ;  plot(qq)
                NBaseOncyl(k)=length(Cyl_h(k).XData);
            end
            title('Base to base distance on cylinders.'); xlabel('base');ylabel('distance (nm)')
            subplot(3,1,2) ; hold on ;
            for k =  1: nCylinder
                plot(ds_Cyli_ForBundle(k,:) ,'--');
            end
            plot(mean(ds_Cyli_ForBundle ) ,'LineWidth',2)
            title('Cylinder lengths of each section, Bold for neutral.'); xlabel('bundle');ylabel('Length (nm)');
            
            subplot(3,1,3) ; hold on ;
            bar(NBaseOncyl)
            title('Entire lengths of cylinders'); xlabel('cylinder');ylabel('N (bps)')
            
            drawnow ;
            %             sdfs=3
        end
        function  [New_tarr, NewNodes]=SlicingCoarseCut(obj ,Total_split_ForEachSlice,Scale_C ,N_slice ,t_arrs)
%             figure(2426);
%             Total_split_ForEachSlice = 20 ;  % number of bundles in Assembly, will be round up due to non-zeros case            
            
%             Scale_C;Position_TanVec_NormalVec;N_slice; t_arrs;
            fnprime = fnder(Scale_C);    Lfun = @(s) sqrt(sum(fnval(fnprime,s).^2,1));
            AvgCurvature = zeros(N_slice-1 ,1)  ;  %Weight for slitting
            for k =1:N_slice-1
                dT=fnval(fnprime, t_arrs(k:k+1) ) ;
                dT=dT(:,2)-dT(:,1) ;
                ds =    integral(Lfun,t_arrs(k),t_arrs(k+1)) ;
                AvgCurvature(k)=norm(dT)/ds ;
            end
            
            Nb = ceil(AvgCurvature/sum(AvgCurvature)*Total_split_ForEachSlice ) ;
            New_tarr = [] ;
            for k=1:N_slice-1
                New_tarr= union(New_tarr, linspace(t_arrs(k),t_arrs(k+1),Nb(k)+1 )   )  ;
            end
            NewNodes = fnval(Scale_C , New_tarr ) ;
            NewNodes=NewNodes';  % use this to detect closest point on cylinders
%             sH=scatter3(NewNodes(:,1),NewNodes(:,2) ,NewNodes(:,3),64 ,'filled'  );
%             sH.CData=[0 1 1];
            %-----------

        end
        function  [DotOnCylinders,ds_Cyli_ForBundle ] =SlicingCoarseCut_part2(obj,nCylinder,NewNodes,Cyl_h)
            
            DotOnCylinders = zeros(nCylinder , size(NewNodes,1) ,3 ) ;
            ds_Cyli_ForBundle=  zeros(nCylinder , size(NewNodes,1)-1  ) ;  % unit: nm
            for k_cyl =  1: nCylinder
                XYZ_onCyl = [Cyl_h(k_cyl).XData ;Cyl_h(k_cyl).YData; Cyl_h(k_cyl).ZData]' ;
                for j_spl = 1 :size(NewNodes,1)
                    dXYZ = XYZ_onCyl- ones(size(XYZ_onCyl,1),1)*NewNodes(j_spl,:) ;
                    d =  sqrt(sum(dXYZ.^2 ,2)) ;
                    Ind = find(d==min(d)) ;Ind=Ind(1) ;
                    DotOnCylinders(k_cyl, j_spl ,1) =  XYZ_onCyl(Ind,1) ;
                    DotOnCylinders(k_cyl, j_spl ,2) =  XYZ_onCyl(Ind,2) ;
                    DotOnCylinders(k_cyl, j_spl ,3) =  XYZ_onCyl(Ind,3) ;
                    if j_spl~=1  %numerical integral on the cylinder curve, instead of end-to-end distance
                        ds = XYZ_onCyl(PrevInd+1:Ind,:) -  XYZ_onCyl(PrevInd:Ind-1,:) ;
                        Total_ds= sum( sqrt(sum(ds.^2 ,2)) ) ;
                        ds_Cyli_ForBundle(k_cyl, j_spl-1) =Total_ds ;
                    end
                    PrevInd = Ind ;
                end
            end
%             AllJoints = reshape(DotOnCylinders,size(DotOnCylinders,1)*size(DotOnCylinders,2),3) ;
%             scatter3(AllJoints(:,1),AllJoints(:,2),AllJoints(:,3),'r.') ;
        end
        function simplifyLine_byGtheory(obj, Cyl_h)
            InterVal_D = 0.34  ; % nm
            for k =1: length(Cyl_h)
                XYZ = [Cyl_h(k).XData; Cyl_h(k).YData; Cyl_h(k).ZData ]';                                size(XYZ);
                Dxyz = diff(XYZ) ;   Dist= sqrt(sum(Dxyz.^2 ,2)) ;
                cumS_D = cumsum(Dist) ;
                IntegerInds = ceil(cumS_D/InterVal_D) ;
                EquPre= IntegerInds(1:end-1)~=IntegerInds(2:end) ;
                XYZ=XYZ([EquPre;false],:) ;   size(XYZ)   ;
                %--------------
                %
                N_movAvg = 20 ;
                XYZ=[ repmat(XYZ(1,:),N_movAvg-1, 1) ;XYZ ;  repmat(XYZ(end,:),N_movAvg-1, 1) ] ;
                
                XYZ= movmean( XYZ ,20  ,'Endpoints','shrink') ;
                XYZ= [XYZ(1,:);  XYZ(N_movAvg+1: end-N_movAvg,:); XYZ(end,:)] ;
                %-----------
%                 if  k==1  %|| k==length(Cyl_h)

% SaveXYZ =XYZ ;
                    cyl_spline=  cscvn( XYZ') ;
                    %                 cyl_spline = csape({XYZ(:,1),XYZ(:,2),XYZ(:,3)},'second') ;
                    
                    fnprime = fnder(cyl_spline);
                    Lfun = @(s) sqrt(sum(fnval(fnprime,s).^2,1));
                    L = integral(Lfun,cyl_spline.breaks(1),cyl_spline.breaks(end)) ; % nm
                    
                    N_base  = round(L/InterVal_D );
                    s_i = linspace(cyl_spline.breaks(1),cyl_spline.breaks(end) , N_base) ;
%                     size(s_i)
                    for rp=1:5
                    dxyz_ini =diff(fnval(cyl_spline , s_i)'  ,1) ;
                    ds_ini =  sqrt(sum(dxyz_ini.^2 ,2) ) ;
                    Cum= cumsum(ds_ini) ; 
                    Cum2 = [ 0;Cum ]; 
                    EqSum = linspace(Cum2(1) ,Cum2(end) , length(Cum2)) ;                    
                    S_i_eq = interp1(Cum2,s_i(1:end),EqSum) ;
                    s_i= 0.5*(S_i_eq +s_i) ;
%                     size(s_i)
                    end
                    
%                     dxyz_eq =diff(fnval(cyl_spline , [ S_i_eq])'  ,1) ;
%                     ds_eq =  sqrt(sum(dxyz_eq.^2 ,2) ) ;
                    
%                     vq = interp1(x,v,xq)
                    XYZ=fnval(cyl_spline , S_i_eq)' ;
%                 end
                %
                Cyl_h(k).XData=XYZ(:,1)'; Cyl_h(k).YData=XYZ(:,2)';  Cyl_h(k).ZData=XYZ(:,3)';
            end
            %-------------
        end
        
        function ShowInd(obj,src,evn)
            
            XYZ =[src.XData ; src.YData ;src.ZData]' ;
            ClicP = evn.IntersectionPoint ;
            d = XYZ- ones(length(src.XData),1)*ClicP ;
            dd = sum(d.^2 ,2) ;
            Ind = find(dd==min(dd)) ;
            %             sdfs=3
        end
        
        
        function VisualizeFFC(obj)
            % show points and curve, allow user control points
            obj.fH = figure(26634) ;
            obj.fH.Name='Spline sketch GUI' ; set( obj.fH,'NumberTitle','off');
            obj.fH.Units='normalized'; obj.fH.OuterPosition=[0 0 1 1];  clf;

            ax= gca ; ax.Position = [0.05 0.1 0.75 0.8 ] ;
            obj.CurveParas=[] ;
            for k = 1 :length(obj.PointOnCurve)
                obj.Cscatter{k} = scatter3(obj.PointOnCurve{k}(:,1),obj.PointOnCurve{k}(:,2) ,obj.PointOnCurve{k}(:,3),58,'filled'  );
                obj.Cscatter{k}.CData= repmat( [0 0.45 0.75 ],length(obj.PointOnCurve{k}(:,1)),1  );
                %             obj.Cscatter.ButtonDownFcn
                hold on ;
                [obj.CurveParas.p{k},obj.CurveParas.t{k}] = fnplt(obj.SplineCurve{k}) ;
                obj.Cline{k}=  plot3( obj.CurveParas.p{k}(1,:),obj.CurveParas.p{k}(2,:),obj.CurveParas.p{k}(3,:) );
                obj.Cline{k}.UserData.splinePlotXYZ = obj.CurveParas.p{k};
                obj.Cline{k}.HitTest='off' ;
            end
            
            %---------GUI components
            obj.editH = uicontrol(obj.fH,'Style', 'edit','String', '1','Unit','normalized','Position', [0.85 0.92 0.04 0.05]);  
            
            obj.TRtxt = uicontrol(obj.fH,'Style','text','Unit','normalized','FontSize',14,'Position', [0.91 0.92 0.06 0.05],....
                'String','Interval','HorizontalAlignment','left');  % [0.92 0.85 0.05 0.05]

            
            Str= { 'Move node' , 'Insert Node','Delete Node' ,'Break' ,'Rotate node'}  ; % option
            obj.popupH = uicontrol(obj.fH,'Style', 'popupmenu',...
                'String', Str,'Unit','normalized','Position', [0.85 0.8 0.1 0.05] ,'FontSize',12);    %            
            
            obj.popupH.Callback= @(src,evn)SwitchMode(obj,src,evn) ;
            obj.popupH2 = uicontrol(obj.fH,'Style', 'popupmenu',...
                'String', 'Spline: 1','Unit','normalized','Position', [0.85 0.7 0.1 0.05] ,'FontSize',12 );    %
            obj.popupH2.Callback= @(src,evn)SelectLine(obj,src,evn) ;
            obj.check_keyControl = uicontrol(obj.fH,'Style', 'checkbox',...
                'String', 'Move point/View','Value',1,'Unit','normalized','Position', [0.85 0.65 0.1 0.05] ,'FontSize',12);    %
            
            obj.check_keyControl.Callback = @(src,evn)ChangeKeyControl(obj,src,evn )  ;
            obj.check_SWsplineStraight = uicontrol(obj.fH,'Style', 'checkbox',...
                'String', 'Switch Spline or Straight','Value',1,'Unit','normalized','Position', [0.85 0.68 0.1 0.05] ,'FontSize',12);    %
            obj.check_SWsplineStraight.Callback = @(src,evn)SwitchRepresentation(obj,src,evn )  ;
            
            align([obj.editH obj.popupH obj.popupH2 obj.check_keyControl obj.check_SWsplineStraight],'Left','distribute');
            
            
            obj.UIgroup_XYZ_AlignProj = uibuttongroup(obj.fH,'Position',[0.85 0.38 0.12 0.24],'Title','Align / Projection','FontSize',12,'BorderWidth',2) ;
            obj.btn1_projXY = uicontrol(obj.UIgroup_XYZ_AlignProj,'Style', 'pushbutton',...
                'String', 'XY','Unit','normalized','Position', [0.1 0.05 0.8 0.25] ,'FontSize',12);    %
            obj.btn2_projXZ = uicontrol(obj.UIgroup_XYZ_AlignProj,'Style', 'pushbutton',...
                'String', 'XZ','Unit','normalized','Position', [0.1 0.5 0.8 0.25] ,'FontSize',12);    %
            obj.btn3_projYZ = uicontrol(obj.UIgroup_XYZ_AlignProj,'Style', 'pushbutton',...
                'String', 'YZ','Unit','normalized','Position', [0.1 0.65 0.8 0.25] ,'FontSize',12 );    %            
%             obj.btn1_projXY = uicontrol(obj.fH,'Style', 'pushbutton',...
%                 'String', 'XY','Unit','normalized','Position', [0.9 0.05 0.05 0.05] );    %
%             obj.btn2_projXZ = uicontrol(obj.fH,'Style', 'pushbutton',...
%                 'String', 'XZ','Unit','normalized','Position', [0.9 0.12 0.05 0.05] );    %
%             obj.btn3_projYZ = uicontrol(obj.fH,'Style', 'pushbutton',...
%                 'String', 'YZ','Unit','normalized','Position', [0.9 0.19 0.05 0.05] );    %
            obj.btn1_projXY.Callback= @(src,evn)Proj_plane(obj,src,evn) ;
            obj.btn2_projXZ.Callback= @(src,evn)Proj_plane(obj,src,evn) ;
            obj.btn3_projYZ.Callback= @(src,evn)Proj_plane(obj,src,evn) ;
            
            AssignIcon(  obj.btn1_projXY,'XY.jpg' ) ;  obj.btn1_projXY.TooltipString='Align or project to XY plane.' ;
            AssignIcon(  obj.btn2_projXZ,'XZ.jpg' ) ;  obj.btn2_projXZ.TooltipString='Align or project to XZ plane.' ;
            AssignIcon(  obj.btn3_projYZ,'YZ.jpg' ) ;  obj.btn3_projYZ.TooltipString='Align or project to YZ plane.' ;
            obj.btn1_projXY.UserData.String='XY' ; obj.btn2_projXZ.UserData.String='XZ' ; obj.btn3_projYZ.UserData.String='YZ' ;
            
            align([obj.btn3_projYZ obj.btn2_projXZ obj.btn1_projXY ],'Left','distribute');
            

            obj.UIgroup_Other = uibuttongroup(obj.fH,'Position',[0.85 0.05 0.12 0.3],'Title','Spline Operations','FontSize',12,'BorderWidth',2) ;
            obj.btn4_Save = uicontrol(obj.UIgroup_Other,'Style', 'pushbutton','String', 'Save',...
                'Unit','normalized','Position', [0.1 0.25 0.8 0.2], 'Callback',@(src,evn)SaveXYZ(obj,src,evn) ,'FontSize',12  );    %
            obj.btn5_Load = uicontrol(obj.UIgroup_Other,'Style', 'pushbutton','String', 'Load',...
                'Unit','normalized','Position', [0.1 0.45 0.8 0.2], 'Callback',@(src,evn)LoadXYZ(obj,src,evn) ,'FontSize',12   );    %
            obj.btn6_Export = uicontrol(obj.UIgroup_Other,'Style', 'pushbutton','String', 'Export',...
                'Unit','normalized','Position', [0.1 0.05 0.8 0.2], 'Callback',@(src,evn)ExportXYZ(obj,src,evn) ,'FontSize',12   );    %
            obj.btn7_Add = uicontrol(obj.UIgroup_Other,'Style', 'pushbutton','String', 'Add',...
                'Unit','normalized','Position', [0.1 0.75 0.8 0.2], 'Callback',@(src,evn)AddSpline(obj,src,evn) ,'FontSize',12   );    %
            AssignIcon(  obj.btn4_Save,'SpSave.jpg' ) ;  obj.btn4_Save.TooltipString='Save splines(s).' ;
            AssignIcon(  obj.btn5_Load,'SpLoad.jpg' ) ;  obj.btn5_Load.TooltipString='Load splines(s).' ;
            AssignIcon(  obj.btn6_Export,'SpExport.jpg' ) ;  obj.btn6_Export.TooltipString='Export paths for extrude or sweep.' ;
            AssignIcon(  obj.btn7_Add,'SpAdd.jpg' ) ;  obj.btn7_Add.TooltipString='Add Spline(s).' ;
            align([obj.btn7_Add obj.btn6_Export obj.btn5_Load obj.btn4_Save ],'Left','distribute');
            
%             obj.btn4_Save = uicontrol(obj.fH,'Style', 'pushbutton','String', 'Save',...
%                 'Unit','normalized','Position', [0.9 0.30 0.05 0.05], 'Callback',@(src,evn)SaveXYZ(obj,src,evn)   );    %
%             obj.btn5_Load = uicontrol(obj.fH,'Style', 'pushbutton','String', 'Load',...
%                 'Unit','normalized','Position', [0.9 0.37 0.05 0.05], 'Callback',@(src,evn)LoadXYZ(obj,src,evn)    );    %
%             obj.btn6_Export = uicontrol(obj.fH,'Style', 'pushbutton','String', 'Export',...
%                 'Unit','normalized','Position', [0.9 0.44 0.05 0.05], 'Callback',@(src,evn)ExportXYZ(obj,src,evn)    );    %
%             obj.btn7_Add = uicontrol(obj.fH,'Style', 'pushbutton','String', 'Add',...
%                 'Unit','normalized','Position', [0.9 0.51 0.05 0.05], 'Callback',@(src,evn)AddSpline(obj,src,evn)    );    %        
            
            %---------add legend as instructions, single row
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            hLg= legend(ForLegend,'Click me for instructions','Location','northwest') ;

            hLg.String={'\bf{Click me for instructions}'};
            hLg.Interpreter='tex';  %latex
            hLg.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_Spline( src,evn,ax );

            hLg.Units='normalized'; hLg.AutoUpdate ='off';
            hLg.Position=[0.0063 0.9728 0.1569 0.025];
            %------------------



            for k = 1 :length(obj.PointOnCurve)
                obj.Cscatter{k}.ButtonDownFcn =@(src,evn)MoveCscatter(obj,src,evn) ;
            end
            obj.fH.WindowKeyPressFcn =@(src,evn)windowKeyIn(obj,src,evn) ;
            axis equal ;
            xlabel('X'); ylabel('Y') ;zlabel('Z');
            
            drawnow ;
            SwitchMode(obj,obj.popupH,[]) ;
            UpdateData(obj) ;
            
        end
        
        function AddSpline(obj,src,evn)
            
            [FileName,PathName] = uigetfile('*.txt;*.xyz','Select the file with XYZ points');
            fffid=fopen(strcat(PathName,FileName));
            Oneline = fgetl(fffid) ;   Node= [];   WithHeader= 0;
            if  contains(Oneline,'HeadNode' )
                Node = str2num(Oneline(9:end) ) ;     WithHeader=1 ;
            else  % no head, pure XYZ, reset reading
                fclose(fffid);   fffid=fopen(strcat(PathName,FileName));
            end
            %------
            count=0;
            while 1
                Oneline = fgetl(fffid) ;   count=count+1 ;
                if isempty(Oneline)  % ||strcmp(Oneline, '-1')  ||count>20
                    break
                end
                if isa(Oneline,'double')
                    if Oneline== -1
                        break
                    end
                end
            end
            fclose(fffid);
            %-------------------------
            XYZ = zeros( count-1,3) ;
            fffid2=fopen(strcat(PathName,FileName)); count=0 ;
            if WithHeader==1;  Oneline = fgetl(fffid2) ;   end
            while 1
                Oneline = fgetl(fffid2) ;    count=count+1 ;
                if ~isa(Oneline,'double')
                    XYZ(count ,: ) = str2num( Oneline) ;
                end
                if isempty(Oneline)  % ||strcmp(Oneline, '-1')  ||count>20
                    break
                end
                if isa(Oneline,'double')
                    if Oneline== -1
                        break
                    end
                end
            end
            fclose(fffid2);
            %---------
            obj.ls=1 ;
            if ~isempty(Node)
                TailNode = Node(2:end)-1 ; TailNode=[TailNode , size(XYZ,1)] ;
                for k =1: length(Node)
                    obj.PointOnCurve{end+1} = XYZ(Node(k):TailNode(k) ,:  ) ;
                    obj.SplineCurve{end+1} =  cscvn( obj.PointOnCurve{end}(1:end,: )') ;
                end
            else
                obj.PointOnCurve{end+1} = XYZ ;
                obj.SplineCurve{end+1} =  cscvn( XYZ') ;
            end
            obj.VisualizeFFC ;
        end
        
        function  ChangeKeyControl(obj,src,evn )
            if src.Value==1
                obj.fH.WindowKeyPressFcn= @(src1,evn1)windowKeyIn(obj,src1,evn1) ;
            else
                obj.fH.WindowKeyPressFcn= @(src2,evn2)ChangeAxesLimit(obj,src2,evn2) ;
            end
        end
        function tableBtnDown(obj,src, evn)
            Ls = evn.Indices(1) ;
            for k=1: length(obj.Cscatter)
                obj.Cscatter{k}.HitTest='off' ;   obj.Cline{k}.HitTest='off' ;    obj.Cline{k}.LineStyle=':' ;
            end
            obj.Cline{Ls}.LineStyle='-' ;
            XYZ = [obj.Cline{Ls}.XData([1 end]) ; obj.Cline{Ls}.YData([1 end]);obj.Cline{Ls}.ZData([1 end]) ];
            if ~isfield(src.UserData, 'hasplot')
                src.UserData.pp = plot3(findobj(obj.fH,'type','Axes'),XYZ(1,:),XYZ(2,:),XYZ(3,:),'-k', 'LineWidth',2 ) ;
                src.UserData.hasplot =1 ;
            else
                src.UserData.pp.XData = XYZ(1,:) ;     src.UserData.pp.YData = XYZ(2,:) ;      src.UserData.pp.ZData = XYZ(3,:) ;
            end
        end
        
        function AskWhichOption(obj) 
            % decision.pair=2; 'menubar','none', pixels
%             decision.whichBundle=1;'menubar','none',
%             obj.fH.Units='normalized';
            fh = figure('units','pixels','position',[300 300 600 250],'menubar','none',...
                'name','Advanced options','Units','normalized',...
                'numbertitle','off', 'resize','off');
            movegui(fh,'center') ;
            
            str = {'Mirror to plane:' ,'Axial-symmetric to axis:' ,'Insert break nodes', 'Equi-distance nodes','None'} ;
            Op.pop1 = uicontrol(fh,'style','pop','Value',1,...
                'units','normalized','position',[0.3 0.6 0.3 0.1],    'String',str);
            
            str2 = {'XY-plane / Z-axis' ,'YZ-plane / X-axis' , 'XZ-plane/ Y-axis'} ;
            Op.pop2 = uicontrol(fh,'style','pop','Value',1, ...
                'units','normalized','position',[0.3 0.2 0.3 0.1],    'String',str2);
            
            UIgroup1 = uibuttongroup(fh,'Position',[0.62 0.2 0.35 0.5],'Title','Applied to ','FontSize',12,'BorderWidth',2) ;

            Op.Radio_SingleSpline=uicontrol('Style','radiobutton','Parent',UIgroup1,'Units','normalized','Position',[0.05 0.2 0.8 0.2],'String','One Spline'  );
            Op.Radio_AllSpline=uicontrol('Style','radiobutton','Parent',UIgroup1,'Units','normalized','Position',[0.05 0.7 0.8 0.2],'String','All Splines'  );
            
%             txtH2 = uicontrol(fh,'Style','text','FontSize',12,'Units','normalized',...
%             'Position', [0.3 0.5 0.2 0.2], 'String','Select options');
            
            btn = uicontrol(fh,'style','pushbutton','Parent',fh,'Units','normalized',...
                'unit','normalized',    'position',[0.02 0.3 0.2 0.2],...
                'string','OK',  'callback',@(src,evn)DecideOption(obj,src,evn ,fh,Op)    );
            
            h = findobj(fh,'type','uicontrol');  set(h, 'FontSize' ,12) ;
            uiwait(fh);
        end
        function AskEquiDistance(obj)
            fh = figure('units','pixels','position',[300 300 300 250],'menubar','none',...
                'name','Advanced options','Units','normalized', 'numbertitle','off', 'resize','off');        
            movegui(fh,'center') ;
        
            Str={'--'};
            for k=1:length( obj.Cscatter) ;  Str{k}=strcat('Spline : ',num2str(k));    end
            pop_Whichspline = uicontrol(fh,'style','pop','Value',1,...
                'units','normalized','position',[0.2 0.6 0.3 0.1],    'String',Str);
            Edit_Whichspline = uicontrol(fh,'style','edit',...
                'units','normalized','position',[0.55 0.6 0.3 0.1],    'String','30');
            btn_OK = uicontrol(fh,'style','pushbutton','Parent',fh,'Units','normalized',...
                'unit','normalized',    'position',[0.3 0.2 0.4 0.2],...
                'string','OK',  'callback',@(src,evn)DecideEquiNode(obj,src,evn,pop_Whichspline,Edit_Whichspline,fh )    );
%                     obj.editH = uicontrol(obj.fH,'Style', 'edit','String', '1','Unit','normalized','Position', [0.85 0.92 0.04 0.05]);  
            h = findobj(fh,'type','uicontrol');  set(h, 'FontSize' ,12) ;
            uiwait(fh);
        end
        function DecideEquiNode(obj,src,evn,pop_Whichspline,Edit_Whichspline ,fh)
            Ind_Spline =   pop_Whichspline.Value  ;
            N_point = str2num(Edit_Whichspline.String )  ;
            obj.get_Equi_Distance_Point(Ind_Spline,N_point) ;
            obj.VisualizeFFC ;  close(fh) ;
        end
        
        
        function DecideOption(obj,src,evn , fOption,Op)
            if Op.Radio_SingleSpline.Value ==1
                SingleSpline =true ;
            else
                SingleSpline =false ; 
            end
            
            if Op.pop1.Value==1
                Sym ='Mirror' ;
            elseif Op.pop1.Value==2
                Sym ='Assymetric' ;
            elseif Op.pop1.Value==3
                obj.InsertBreakNodeBetween(SingleSpline) ;
                obj.UpdateData ;
                return ;
            elseif Op.pop1.Value==4
                obj.AskEquiDistance ;
                obj.UpdateData ;
                close(fOption) ;
                return ;
            else
                close(fOption) ;
                return ;
            end
%             str2 = {'XY-plane / Z-axis' ,'YZ-plane / X-axis' , 'XZ-plane/ Y-axis'} ;
            if Op.pop2.Value==1
                strXYZ = {'ZData','XData','YData'} ;
            elseif Op.pop2.Value==2
                strXYZ = {'XData','YData','ZData'} ;
            else
                strXYZ = {'YData','XData','ZData'} ;
            end
            
            obj.ApplySymmetric(Sym, strXYZ ,SingleSpline ) ;
            obj.UpdateData ;
%             close(fOption) ;
        end
        
        function InsertBreakNodeBetween(obj,SingleSpline)            
            for j =1: length(obj.Cscatter)
                if SingleSpline==1 && j ~=obj.ls
                    continue ;
                end
                Ori_Sarr = obj.SplineCurve{j}.breaks ;
                InsertMid = movmean(Ori_Sarr ,2 ,'Endpoints','discard') ;                
                New_Sarr = sort([Ori_Sarr ,InsertMid ] ) ;
                NewXYZ = fnval(obj.SplineCurve{j} , New_Sarr) ;
                obj.Cscatter{j}.XData = NewXYZ(1,:) ;
                obj.Cscatter{j}.YData = NewXYZ(2,:) ;
                obj.Cscatter{j}.ZData = NewXYZ(3,:) ;
            end
        end
        
        
        function ApplySymmetric(obj , MirrorOrAssymetric , XYZcase,SingleSpline)
            
%             MirrorOrAssymetric ='Mirror' ;
            %                     MirrorOrAssymetric ='Assymetric' ;
            
            strXYZ=XYZcase ;
%             strXYZ = {'XData','YData','ZData'} ;
            %                      strXYZ = {'ZData','YData','XData'} ;
            
            for j =1: length(obj.Cscatter)
                if SingleSpline==1 && j ~=obj.ls
                    continue ;
                end
                for k = 1: round(length(obj.Cscatter{j}.XData)/2)
                    
                    if strcmp(MirrorOrAssymetric , 'Mirror' )
                        meanABSX = 0.5* (abs(obj.Cscatter{j}.(strXYZ{1})(k)) +abs(obj.Cscatter{j}.(strXYZ{1})(end-k+1)) ) ;
                        if obj.Cscatter{j}.(strXYZ{1})(k)<0
                            obj.Cscatter{j}.(strXYZ{1})(k)= -meanABSX ;
                            obj.Cscatter{j}.(strXYZ{1})(end-k+1) =meanABSX ;
                        else
                            obj.Cscatter{j}.(strXYZ{1})(k)= meanABSX ;
                            obj.Cscatter{j}.(strXYZ{1})(end-k+1) = -meanABSX ;
                        end
                    else
                         meanABSX = 0.5* (abs(obj.Cscatter{j}.(strXYZ{1})(k)) +abs(obj.Cscatter{j}.(strXYZ{1})(end-k+1)) ) ;
                         obj.Cscatter{j}.(strXYZ{1})(k)= meanABSX ;obj.Cscatter{j}.(strXYZ{1})(end-k+1) =meanABSX ;
                        
                    end
                    
                    if strcmp(MirrorOrAssymetric , 'Mirror' )
                        meanZ = 0.5* ((obj.Cscatter{j}.(strXYZ{3})(k)) +(obj.Cscatter{j}.(strXYZ{3})(end-k+1)) ) ;
                        obj.Cscatter{j}.(strXYZ{3})(k)= meanZ ;  obj.Cscatter{j}.(strXYZ{3})(end-k+1) =meanZ ;
                        
                        meanY = 0.5* ((obj.Cscatter{j}.(strXYZ{2})(k)) +(obj.Cscatter{j}.(strXYZ{2})(end-k+1)) ) ;
                        obj.Cscatter{j}.(strXYZ{2})(k)= meanY ;  obj.Cscatter{j}.(strXYZ{2})(end-k+1) =meanY ;
                    else
                        meanZ = 0.5* (abs(obj.Cscatter{j}.(strXYZ{3})(k)) + abs(obj.Cscatter{j}.(strXYZ{3})(end-k+1)) ) ;
                        meanY = 0.5* (abs(obj.Cscatter{j}.(strXYZ{2})(k)) + abs(obj.Cscatter{j}.(strXYZ{2})(end-k+1)) ) ;
                        if obj.Cscatter{j}.(strXYZ{3})(k)<0
                            obj.Cscatter{j}.(strXYZ{3})(k)= -meanZ ;  obj.Cscatter{j}.(strXYZ{3})(end-k+1) =meanZ ;
                        else
                            obj.Cscatter{j}.(strXYZ{3})(k)= meanZ ;   obj.Cscatter{j}.(strXYZ{3})(end-k+1) = -meanZ ;
                        end
                        if obj.Cscatter{j}.(strXYZ{2})(k)<0
                            obj.Cscatter{j}.(strXYZ{2})(k)= -meanY ;  obj.Cscatter{j}.(strXYZ{2})(end-k+1) =meanY ;                            
                        else
                            obj.Cscatter{j}.(strXYZ{2})(k)= meanY ;   obj.Cscatter{j}.(strXYZ{2})(end-k+1) = -meanY ;
                        end
                        
                    end
                end
                if mod(length(obj.Cscatter{j}.XData) ,2 )==1
                    obj.Cscatter{j}.(strXYZ{1})(k) =0 ;
                    if ~strcmp(MirrorOrAssymetric , 'Mirror' )
                        obj.Cscatter{j}.(strXYZ{3})(k) =0 ;
                        obj.Cscatter{j}.(strXYZ{2})(k) =0 ;
                    end
                end
            end
        end
        
        
        
        function ExportTableChange(obj,src,evn)
            Ls = evn.Indices(1) ;
            if src.Data{Ls,2}==1
                src.UserData.pp.LineStyle='-' ;
            else
                src.UserData.pp.LineStyle=':' ;
            end
        end        
        function DecideConnectLine(obj,src,evn,t,btn,NewFig)
            btn.UserData.ConnectDataInTable = t.Data ;
            close(NewFig) ;
        end
        
        function ExportXYZ(obj,src,evn)
            NewFig = dialog('units','pixels', 'position',[300 300 400 500],...
                'menubar','none',  'name','Connect head to end or not',...
                'numbertitle','off', 'resize','off');
            btn =src ;
            movegui(NewFig,'center')
            N=length(obj.PointOnCurve);
            
            Mat=[ (1:N)'  ,zeros(N,1)]  ;
            Datat=mat2cell(Mat, ones(1,N ) , [1,1]  );
            for rev=1:size(Datat,1)
                Datat{rev,1}= num2str(Datat{rev,1});        Datat{rev,2}=false;
            end
            
            t = uitable(NewFig,'Units','normalized','Position',[0.05 0.2 0.9 0.75],....
                'ColumnWidth','auto','ColumnFormat',({[] [] }),...
                'ColumnEditable', true,'Data',Datat,....
                'ColumnName',{'Spline'; 'Connect/disconnect' });
            btn_ExportNewHB = uicontrol(NewFig,'Style', 'pushbutton', 'String', '	OK ','Unit','normalized', 'Position', [0.5 0.05 0.4 0.1] );
            t.CellSelectionCallback =@(src,evn)tableBtnDown(obj,src, evn) ;
            t.CellEditCallback=  @(src,evn)ExportTableChange(obj,src, evn);
            btn_ExportNewHB.Callback = @(src,evn)DecideConnectLine(obj,src,evn,t,btn,NewFig ) ;
            uiwait(NewFig) ;
            %             finishAsking =1
            if isempty(obj.ss_STEP)
                return
            end
            
%             src.UserData.ConnectDataInTable;
            %             return
            [A,~] = cellfun(@size,obj.PointOnCurve) ;
            HeadNodes= cumsum([1 ,A] ) ; HeadNodes=HeadNodes(1:end-1) ;
            TailNode = HeadNodes(2:end)-1 ; TailNode=[TailNode , sum(A)] ;
            
            XYZ = zeros(sum(A),3); c=1;
            for k =1:length(obj.PointOnCurve)
                XYZ( c:A(k)+c-1 ,: )   = obj.PointOnCurve{k} ; c=c+A(k);
            end
            obj.ss_STEP.UserData.UsePoints.pXYZ=XYZ;
            
            QQ =[] ;
            for k=1:length(A)
                OneLine = [    HeadNodes(k):TailNode(k)-1 ; HeadNodes(k)+1:TailNode(k)]' ;
                QQ=[QQ;OneLine ] ;
                
                if src.UserData.ConnectDataInTable{k,2} ==1
                    QQ=[QQ;OneLine(1) OneLine(end) ] ;    % looped spline
                end
            end
            obj.ss_STEP.UserData.UsePoints.Edges= QQ ;
            obj.ss_STEP.UserData.FFC = obj;
            delete(obj.fH) ;
        end
        
        function Proj_plane(obj,src,evn)
            IndSelect = obj.Cscatter{obj.ls}.UserData.Ind  ;
            if length(IndSelect)==1 || isempty(IndSelect)
                % if only select one or nothing, use button to project
                % ALL points to the planes (xy, xz, or yz planes).
                switch src.UserData.String
                    case 'XY'
                        obj.Cscatter{obj.ls}.ZData=zeros(size( obj.Cscatter{obj.ls}.ZData));
                    case 'XZ'
                        obj.Cscatter{obj.ls}.YData=zeros(size( obj.Cscatter{obj.ls}.YData));
                    case 'YZ'
                        obj.Cscatter{obj.ls}.XData=zeros(size( obj.Cscatter{obj.ls}.XData));
                end
            else
                % if select multiple, use buttons to align selected nodes.
                switch src.UserData.String
                    case 'XY'
                     obj.Cscatter{obj.ls}.ZData(IndSelect)= mean( obj.Cscatter{obj.ls}.ZData(IndSelect))*ones(size(IndSelect))   ;
                    case 'XZ'
                     obj.Cscatter{obj.ls}.YData(IndSelect)= mean( obj.Cscatter{obj.ls}.YData(IndSelect))*ones(size(IndSelect))   ;
                    case 'YZ'
                     obj.Cscatter{obj.ls}.XData(IndSelect)= mean( obj.Cscatter{obj.ls}.XData(IndSelect))*ones(size(IndSelect))   ;
                end
            end
            obj.UpdateData ;
        end
        
        function SelectLine(obj,src,evn)
            % callback of popupH2
            obj.ls =src.Value ;
            %             for k= 1: length(obj.Cscatter)
            %               obj.Cline{k}.LineStyle=':' ;
            %             end
            %             obj.Cline{obj.ls}.LineStyle='-' ;
            uistack(obj.Cscatter{obj.ls},'top') ;   uistack(obj.Cline{obj.ls},'top') ;
            SwitchMode(obj,obj.popupH,evn) ;
        end
        
        function SwitchMode(obj,src,evn)
            %             FFC.popupH.ButtonDownFcn
            obj.fH.WindowKeyPressFcn= @(src1,evn1)windowKeyIn(obj,src1,evn1)  ;
            switch src.Value
                case 1  % 'Move node'
                    for k=1: length(obj.Cscatter)
                        obj.Cscatter{k}.HitTest='off' ;
                        obj.Cline{k}.HitTest='off' ; obj.Cline{k}.LineStyle=':' ;
                    end
                    obj.Cscatter{obj.ls}.HitTest='on' ;  obj.Cline{obj.ls}.LineStyle='-' ;
                    delete(obj.TempAddPoint) ;
                     set(obj.TRtxt ,'String','Interval');
                case 2 %  'Insert Node'
                    for k=1: length(obj.Cscatter)
                        obj.Cscatter{k}.HitTest='off' ;
                        obj.Cline{k}.HitTest='off' ; obj.Cline{k}.LineStyle=':' ;
                        obj.Cline{k}.ButtonDownFcn = [];
                    end
                    obj.Cscatter{obj.ls}.HitTest='off' ;
                    obj.Cline{obj.ls}.HitTest='on' ; obj.Cline{obj.ls}.LineStyle='-' ;
                    obj.Cline{obj.ls}.ButtonDownFcn = @(src,evn)InsertNode(obj,src, evn) ;
                case 3 %  'Delete Node'
                    for k=1: length(obj.Cscatter)
                        obj.Cscatter{k}.HitTest='off' ;
                        obj.Cline{k}.HitTest='off' ; obj.Cline{k}.LineStyle=':' ;
                        obj.Cline{k}.ButtonDownFcn = [];
                    end
                    
                    obj.Cscatter{obj.ls}.HitTest='on' ; obj.Cline{obj.ls}.LineStyle='-' ;
                    obj.Cline{obj.ls}.HitTest='off' ;
                    delete(obj.TempAddPoint) ;
                case 4 % 'Break'
                    for k = 1:length( obj.Cscatter)
                        obj.Cscatter{k}.HitTest='off' ;
                        obj.Cline{k}.HitTest='off' ;  obj.Cline{k}.LineStyle=':' ;
                    end
                    obj.Cline{obj.ls}.HitTest='off' ;  obj.Cline{obj.ls}.LineStyle='-' ;
                    obj.Cscatter{obj.ls}.HitTest='on' ;
                    
                    delete(obj.TempAddPoint) ;
                case 5 % 'Rotate'
                    obj.fH.WindowKeyPressFcn= @(src2,evn2)RotateNode(obj,src2,evn2) ;
                    
                    set(obj.TRtxt ,'String','Interval (deg)');
            end
        end
        
        function InsertNode(obj,src, evn)
            d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2  + (src.ZData-evn.IntersectionPoint(3)).^2 ;
            Ind = find(d==min(d))  ;
            %             obj.CurveParas.t{obj.ls}(Ind) ;
            if ~isgraphics(obj.TempAddPoint)
                obj.TempAddPoint =scatter3(src.XData(Ind),src.YData(Ind),src.ZData(Ind) , 56 ,'k') ;
                obj.TempAddPoint.UserData.t =  obj.CurveParas.t{obj.ls}(Ind) ;
            else
                obj.TempAddPoint.XData=src.XData(Ind) ; obj.TempAddPoint.YData=src.YData(Ind) ; obj.TempAddPoint.ZData=src.ZData(Ind) ;
                obj.TempAddPoint.UserData.t =  obj.CurveParas.t{obj.ls}(Ind) ;
            end
            %            g= fnrfn(obj.SplineCurve ,[src.XData(Ind) ,src.YData(Ind), src.ZData(Ind) ])
        end
        
        
        function MoveCscatter(obj,src,evn)
            % CScatter buttondown
            d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2  + (src.ZData-evn.IntersectionPoint(3)).^2 ;
            Ind = find(d==min(d))  ;
            
            allSc = findobj(gca,'type','Scatter') ;
            for k=1:length(allSc)
                if ~isequal(src, allSc(k))  % only clear the non-selected splines
                    allSc(k).CData= repmat( [0 0.45 0.75 ],length( allSc(k).XData),1  );
                    allSc(k).UserData.Ind= [] ;
                end
            end
            
            if evn.Button== 1
                Ind=Ind(1); 
                src.UserData.Ind= Ind ;
            elseif evn.Button== 3
                src.UserData.Ind= union(Ind , src.UserData.Ind ) ;
            elseif  evn.Button== 2
                src.UserData.Ind= union(Ind , src.UserData.Ind ) ;
                src.UserData.Ind = min( src.UserData.Ind):max( src.UserData.Ind) ;
            end
            %              src.UserData.Ind
            src.CData= repmat( [0 0.45 0.75 ],length(obj.PointOnCurve{obj.ls}(:,1)),1  );
            src.CData( src.UserData.Ind,:) = repmat([1,0,0],length(src.UserData.Ind) ,1)  ;
            
            XYZ = [src.XData(Ind) ,src.YData(Ind) ,src.ZData(Ind) ];
            if ~isgraphics(obj.XYZlineOnPoint.x)
                obj.XYZlineOnPoint.x = plot3( [XYZ(1),XYZ(1)+1] ,[XYZ(2),XYZ(2)], [XYZ(3),XYZ(3)] ,'r' ,'HitTest','off') ;
                obj.XYZlineOnPoint.y = plot3( [XYZ(1),XYZ(1)] ,[XYZ(2),XYZ(2)+1], [XYZ(3),XYZ(3)] ,'g','HitTest','off');
                obj.XYZlineOnPoint.z = plot3( [XYZ(1),XYZ(1)] ,[XYZ(2),XYZ(2)], [XYZ(3),XYZ(3)+1] ,'b','HitTest','off')  ;
            else
                obj.XYZlineOnPoint.x.XData=[XYZ(1),XYZ(1)+1]; obj.XYZlineOnPoint.x.YData=[XYZ(2),XYZ(2)]; obj.XYZlineOnPoint.x.ZData=[XYZ(3),XYZ(3)];
                obj.XYZlineOnPoint.y.XData=[XYZ(1),XYZ(1)]; obj.XYZlineOnPoint.y.YData=[XYZ(2),XYZ(2)+1]; obj.XYZlineOnPoint.y.ZData=[XYZ(3),XYZ(3)];
                obj.XYZlineOnPoint.z.XData=[XYZ(1),XYZ(1)]; obj.XYZlineOnPoint.z.YData=[XYZ(2),XYZ(2)]; obj.XYZlineOnPoint.z.ZData=[XYZ(3),XYZ(3)+1];
            end
        end
        
        function RotateNode(obj,src,evn)
            theta = str2num(obj.editH.String) ; theta=theta*pi/180;
            Ind=obj.Cscatter{obj.ls}.UserData.Ind ;
            OriP1 = [obj.Cscatter{obj.ls}.XData(Ind) ;obj.Cscatter{obj.ls}.YData(Ind);obj.Cscatter{obj.ls}.ZData(Ind) ]' ;
            Gcenter= mean(OriP1 ,1) ;
            
%             Gcenter=[0 0 0 ]; % hard code
            %             sdfsdf=3
            switch evn.Key
                case 'q'
                    RMat=[1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
                case 'a'
                    RMat=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
                case 'w'
                    RMat=[cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
                case 's'
                    RMat=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
                case 'e'
                    RMat=[cos(theta) sin(theta) 0; -sin(theta)  cos(theta)  0; 0 0 1];
                case 'd'
                    RMat=[cos(theta) -sin(theta) 0; sin(theta)  cos(theta)  0; 0 0 1];
                otherwise
                    return
            end
            
            OriP1=OriP1 - ones(size(OriP1,1),1)*Gcenter ;
            OriP1=OriP1* RMat ;
            OriP1=OriP1 + ones(size(OriP1,1),1)*Gcenter ;
            
            obj.Cscatter{obj.ls}.XData(Ind) = OriP1(:,1)' ;
            obj.Cscatter{obj.ls}.YData(Ind) = OriP1(:,2)' ;
            obj.Cscatter{obj.ls}.ZData(Ind) = OriP1(:,3)' ;
            
            obj.UpdateData ;
        end
        
        
        function windowKeyIn(obj,src,evn)
            StepSize = str2num(obj.editH.String) ;
            %             [aa,Ind] = ismember([1,0,0],obj.Cscatter.CData  ,'rows') ;
            %             Ind
            %                         evn.Key
            SelectSpline = obj.popupH2.Value ;
            switch evn.Key
                
                case 'o'
                    fprintf('Take the mean of the first and last as the new first node. Remove the last node. \n') ;
                    HeadAndTail_x = mean(obj.Cscatter{SelectSpline}.XData([1 end])) ;
                    HeadAndTail_y = mean(obj.Cscatter{SelectSpline}.YData([1 end])) ;
                    HeadAndTail_z = mean(obj.Cscatter{SelectSpline}.ZData([1 end])) ;
                    obj.Cscatter{SelectSpline}.XData(1)=HeadAndTail_x ;
                    obj.Cscatter{SelectSpline}.YData(1)=HeadAndTail_y ;
                    obj.Cscatter{SelectSpline}.ZData(1)=HeadAndTail_z ;
                    obj.Cscatter{SelectSpline}.XData(end) =[];
                    obj.Cscatter{SelectSpline}.YData(end) =[];
                    obj.Cscatter{SelectSpline}.ZData(end) =[];
                case 'control'
                    
                    obj.AskWhichOption ;
%                     obj.UpdateData ;
                    return
                    %------mirror to XZ plane
                     MirrorOrAssymetric ='Mirror' ;
%                     MirrorOrAssymetric ='Assymetric' ;
                    strXYZ = {'XData','YData','ZData'} ;
                    %                      strXYZ = {'ZData','YData','XData'} ;
                    %                     for j =1 : 1 %length(obj.Cscatter)
                    j= obj.ls;
                    for k = 1: round(length(obj.Cscatter{j}.XData)/2)
                        meanABSX = 0.5* (abs(obj.Cscatter{j}.(strXYZ{1})(k)) +abs(obj.Cscatter{j}.(strXYZ{1})(end-k+1)) ) ;
                        if obj.Cscatter{j}.(strXYZ{1})(k)<0
                            obj.Cscatter{j}.(strXYZ{1})(k)= -meanABSX ;
                            obj.Cscatter{j}.(strXYZ{1})(end-k+1) =meanABSX ;
                        else
                            obj.Cscatter{j}.(strXYZ{1})(k)= meanABSX ;
                            obj.Cscatter{j}.(strXYZ{1})(end-k+1) = -meanABSX ;
                        end
                        
                        if strcmp(MirrorOrAssymetric , 'Mirror' )
                            meanZ = 0.5* ((obj.Cscatter{j}.(strXYZ{3})(k)) +(obj.Cscatter{j}.(strXYZ{3})(end-k+1)) ) ;
                            obj.Cscatter{j}.(strXYZ{3})(k)= meanZ ;
                            obj.Cscatter{j}.(strXYZ{3})(end-k+1) =meanZ ;
                        else
                            meanZ = 0.5* (abs(obj.Cscatter{j}.(strXYZ{3})(k)) + abs(obj.Cscatter{j}.(strXYZ{3})(end-k+1)) ) ;
                            if obj.Cscatter{j}.(strXYZ{3})(k)<0
                                obj.Cscatter{j}.(strXYZ{3})(k)= -meanZ ;
                                obj.Cscatter{j}.(strXYZ{3})(end-k+1) =meanZ ;
                            else
                                obj.Cscatter{j}.(strXYZ{3})(k)= meanZ ;
                                obj.Cscatter{j}.(strXYZ{3})(end-k+1) = -meanZ ;
                            end
                            
                        end
                    end
                    if mod(length(obj.Cscatter{j}.XData) ,2 )==1
                        obj.Cscatter{j}.(strXYZ{1})(k) =0 ;
                        if ~strcmp(MirrorOrAssymetric , 'Mirror' )
                            obj.Cscatter{j}.(strXYZ{3})(k) =0 ;
                        end
                    end
                    %-------
                    
                    %                 case 'x'
                    %                     ratio = 1.01 ;
                    %                     for k =1 : length(obj.Cscatter)
                    %                         obj.Cscatter{k}.XData =   obj.Cscatter{k}.XData*ratio ;
                    %                     end
                    %                     [obj.Cscatter{obj.ls}.XData ; obj.Cscatter{obj.ls}.YData ;obj.Cscatter{obj.ls}.ZData] ;
                    
                case 'm'
                    %                 Yarr = linspace(0, 15 ,13) ;
                    AllXYZ =[0 0 0] ;
                    for k =1 : length(obj.Cscatter)
                        AllXYZ= [ AllXYZ ;   [ obj.Cscatter{k}.XData ;obj.Cscatter{k}.YData;obj.Cscatter{k}.ZData      ]'] ;
                    end
                    mXYZ = mean(AllXYZ) ;
                    for k =1 : length(obj.Cscatter)
                        obj.Cscatter{k}.XData =   obj.Cscatter{k}.XData-mXYZ(1) ;
                        obj.Cscatter{k}.YData =   obj.Cscatter{k}.YData-mXYZ(2) ;
                        obj.Cscatter{k}.ZData =   obj.Cscatter{k}.ZData-mXYZ(3) ;
                    end
                    %                     SDSF=3
                case 't' % print
                    for linei = 1: length(obj.Cscatter)
                        XYZ = [obj.Cscatter{linei}.XData;obj.Cscatter{linei}.YData;obj.Cscatter{linei}.ZData];
                        fprintf('For Spline %i , \n' ,linei)
                        for k=1: length( obj.Cscatter{linei}.XData)
                            fprintf('%i %4.2f %4.2f %4.2f \n',k,XYZ(1:3,k) );
                        end
                        fprintf('\n') ;
                    end
                case 'q'
                    Ind=obj.Cscatter{obj.ls}.UserData.Ind ;
                    obj.Cscatter{obj.ls}.XData(Ind) =  obj.Cscatter{obj.ls}.XData(Ind) + StepSize ;
                case 'a'
                    Ind=obj.Cscatter{obj.ls}.UserData.Ind ;
                    obj.Cscatter{obj.ls}.XData(Ind) =  obj.Cscatter{obj.ls}.XData(Ind) - StepSize ;
                case 'w'
                    Ind=obj.Cscatter{obj.ls}.UserData.Ind ;
                    obj.Cscatter{obj.ls}.YData(Ind) =  obj.Cscatter{obj.ls}.YData(Ind) + StepSize ;
                case 's'
                    Ind=obj.Cscatter{obj.ls}.UserData.Ind ;
                    obj.Cscatter{obj.ls}.YData(Ind) =  obj.Cscatter{obj.ls}.YData(Ind) - StepSize ;
                case 'e'
                    Ind=obj.Cscatter{obj.ls}.UserData.Ind ;
                    obj.Cscatter{obj.ls}.ZData(Ind) =  obj.Cscatter{obj.ls}.ZData(Ind) + StepSize ;
                case 'd'
                    Ind=obj.Cscatter{obj.ls}.UserData.Ind ;
                    obj.Cscatter{obj.ls}.ZData(Ind) =  obj.Cscatter{obj.ls}.ZData(Ind) - StepSize ;
                case 'b'  %break
                    Ind=obj.Cscatter{obj.ls}.UserData.Ind ;
                    PointOnCurveOri = [obj.Cscatter{obj.ls}.XData; obj.Cscatter{obj.ls}.YData; obj.Cscatter{obj.ls}.ZData ];
                    PointsOnOldLine = 1:min(Ind) ;
                    ToNewLine = min(Ind)+1: length(obj.Cscatter{obj.ls}.XData) ;
                    
                    obj.Cscatter{obj.ls}.XData=PointOnCurveOri(1,PointsOnOldLine) ;
                    obj.Cscatter{obj.ls}.YData=PointOnCurveOri(2,PointsOnOldLine) ;
                    obj.Cscatter{obj.ls}.ZData=PointOnCurveOri(3,PointsOnOldLine) ;
                    obj.Cscatter{obj.ls}.CData= repmat( [0 0.45 0.75 ],length(PointOnCurveOri(1,PointsOnOldLine)),1  );
                    
                    obj.Cscatter{end+1}=scatter3(PointOnCurveOri(1,ToNewLine),PointOnCurveOri(2,ToNewLine) ,PointOnCurveOri(3,ToNewLine),58,'filled'  );
                    obj.Cscatter{end}.CData= repmat( [0 0.45 0.75 ],length(PointOnCurveOri(1,ToNewLine)),1  );
                    obj.Cscatter{end}.ButtonDownFcn =@(src,evn)MoveCscatter(obj,src,evn) ;
                    
                    obj.PointOnCurve{end+1} = PointOnCurveOri(:,ToNewLine)' ;
                    obj.SplineCurve{end+1} =  cscvn( obj.PointOnCurve{end}(1:end,: )') ;
                    
                    [obj.CurveParas.p{end+1},obj.CurveParas.t{end+1}] = fnplt(obj.SplineCurve{end}) ;
                    obj.Cline{end+1}=  plot3( obj.CurveParas.p{end}(1,:),obj.CurveParas.p{end}(2,:),obj.CurveParas.p{end}(3,:) );
                    obj.Cline{end}.HitTest='off' ;
                    
                case 'return' % add node
                    if isgraphics(obj.TempAddPoint)
                        NewBreaks = sort( [ obj.SplineCurve{obj.ls}.breaks, obj.TempAddPoint.UserData.t] ) ;
                        NewPointOnCurve = fnval(obj.SplineCurve{obj.ls} , NewBreaks) ;
                        obj.Cscatter{obj.ls}.XData=NewPointOnCurve(1,:) ;obj.Cscatter{obj.ls}.YData=NewPointOnCurve(2,:) ;obj.Cscatter{obj.ls}.ZData=NewPointOnCurve(3,:) ;
                        obj.Cscatter{obj.ls}.UserData.Ind = [] ;
                    end
                    
                case 'backspace' %delete node
                    if obj.popupH.Value==3 && ~isempty( obj.Cscatter{obj.ls}.UserData.Ind )
                        RemainNode = setdiff(1:length(obj.Cscatter{obj.ls}.XData) , obj.Cscatter{obj.ls}.UserData.Ind) ;
                        NewPointOnCurve = [obj.Cscatter{obj.ls}.XData(RemainNode); obj.Cscatter{obj.ls}.YData(RemainNode); obj.Cscatter{obj.ls}.ZData(RemainNode) ];
                        obj.Cscatter{obj.ls}.XData=NewPointOnCurve(1,:) ;obj.Cscatter{obj.ls}.YData=NewPointOnCurve(2,:) ;obj.Cscatter{obj.ls}.ZData=NewPointOnCurve(3,:) ;
                        obj.Cscatter{obj.ls}.UserData.Ind = [] ;
                        obj.PointOnCurve{obj.ls} = [  obj.Cscatter{obj.ls}.XData' ,obj.Cscatter{obj.ls}.YData',obj.Cscatter{obj.ls}.ZData' ] ;
                    end
                    
                    [A,~] = cellfun(@size,obj.PointOnCurve) ;                    
                    if sum(A<=1)>0  % no node or one node
                        IndRR =A>1  ;
                        obj.PointOnCurve=obj.PointOnCurve(IndRR) ;
                        obj.SplineCurve=obj.SplineCurve(IndRR)  ;
                        obj.CurveParas.p= obj.CurveParas.p(IndRR) ;
                        obj.CurveParas.t= obj.CurveParas.t(IndRR) ;
                        obj.Cline= obj.Cline(IndRR) ;
                        obj.Cscatter= obj.Cscatter(IndRR) ;
                        obj.ls =1 ;
                        obj.VisualizeFFC ;
                        obj.popupH2.Value =1;
                    end
                    
                case 'l'  % connect two splines
                    %----query option GUI
                    d = figure('Position',[300 300 500 200],'Name','Show how to connect two lines from ends','numbertitle','off');
                    movegui(d,'center') ;
                    Str =num2str([1:length( obj.Cscatter) ]') ;
                    popup1 = uicontrol('Parent',d,    'Style','popup','Units','normalized',...
                        'Position',[ 0.1750 0.5 0.3 0.25],'FontSize',12,    'String',Str );
                    
                    popup2 = uicontrol('Parent',d,    'Style','popup','Units','normalized',...
                        'Position',[0.525 0.5 0.3 0.25],'FontSize',12,    'String',Str );
                                        
                    Str2 = {'H1-H2','H1-T2' , 'T1-H2' ,'T1-T2'  } ;
                    popupChoice= uicontrol(d,'Style', 'popupmenu',...
                        'String', Str2,'Unit','normalized','Position', [0.35 0.3 0.3 0.2] );    %
                    popupChoice.Callback= @(src,evn)showHowToConnect(obj,src,evn,popup1,popup2) ;
                    btnDecide = uicontrol('Parent',d,...
                        'Units','normalized','Position',[0.25 0.1 0.5 0.1],...
                        'String','Connect','FontSize',12,... ;%...
                        'Callback',{@(src,evn)Connect(obj,src,evn,d,popup1,popup2,popupChoice  )});  %
                case 'p' % duplicate a spline
                    d = figure('Position',[300 300 500 200],'Name','Copy/duplicate spline','numbertitle','off');
                    movegui(d,'center') ;
                    popup1 = uicontrol('Parent',d,    'Style','popup','Units','normalized',...
                        'Position',[ 0.25 0.5 0.5 0.25],'FontSize',12,  'String',num2str([1:length( obj.Cscatter) ]') );
                    popup1.Callback=@(src,evn)duplicateSpline(obj,src,evn , popup1 ) ;
                    btnDecide = uicontrol('Parent',d, 'Units','normalized','Position',[0.25 0.1 0.5 0.1],...
                        'String','Duplicate','FontSize',12 ) ;% , 'Callback',@(src,evn)duplicateSpline(obj,src,evn , popup1) );  %
                    btnDecide.Callback=@(src,evn)DecideDuplicate(obj,src,evn , popup1,d) ;
            end
            obj.UpdateData ;
        end
        
        function duplicateSpline(obj,src,evn,popup1)
            Ls = popup1.Value ;
            for k=1: length(obj.Cscatter)
                obj.Cscatter{k}.HitTest='off' ;   obj.Cline{k}.HitTest='off' ;
                obj.Cline{k}.LineStyle=':' ;
            end
            obj.Cline{Ls}.LineStyle='-' ;
        end
        function DecideDuplicate(obj,src,evn , popup1,d)
            Ls = popup1.Value ;
            PointOnCurveOri = [obj.Cscatter{Ls}.XData; obj.Cscatter{Ls}.YData; obj.Cscatter{Ls}.ZData ];
            ax= findobj(obj.fH,'type','Axes');
            
            obj.Cscatter{end+1}=scatter3(ax,PointOnCurveOri(1,:),PointOnCurveOri(2,:) ,PointOnCurveOri(3,:),58,'filled'  );
            obj.Cscatter{end}.CData= repmat( [0 0.45 0.75 ],length(PointOnCurveOri(1,:)),1  );
            obj.Cscatter{end}.ButtonDownFcn =@(src,evn)MoveCscatter(obj,src,evn) ;
            
            obj.PointOnCurve{end+1} = PointOnCurveOri' ;
            obj.SplineCurve{end+1} =  cscvn( obj.PointOnCurve{end}(1:end,: )') ;
            
            [obj.CurveParas.p{end+1},obj.CurveParas.t{end+1}] = fnplt(obj.SplineCurve{end}) ;
            obj.Cline{end+1}=  plot3(ax , obj.CurveParas.p{end}(1,:),obj.CurveParas.p{end}(2,:),obj.CurveParas.p{end}(3,:) );
            obj.Cline{end}.HitTest='off' ;
            close(d); obj.VisualizeFFC ;
        end
        function showHowToConnect(obj,src,evn,popup1,popup2)
            L1 = popup1.Value;
            L2 = popup2.Value;
            if L1 ==L2 ; return ;end %
            
            for k=1: length(obj.Cscatter)
                obj.Cscatter{k}.HitTest='off' ;   obj.Cline{k}.HitTest='off' ;
                obj.Cline{k}.LineStyle=':' ;
            end
            %            obj.Cscatter{obj.ls}.HitTest='on' ;
            obj.Cline{L1}.LineStyle='-' ;   obj.Cline{L2}.LineStyle='-' ;
            
            d_subUI = src.Parent ;
            LineHead = [ obj.Cline{L1}.XData(1),obj.Cline{L1}.YData(1), obj.Cline{L1}.ZData(1) ;...
                obj.Cline{L2}.XData(1),obj.Cline{L2}.YData(1), obj.Cline{L2}.ZData(1)] ;
            LineTail = [ obj.Cline{L1}.XData(end),obj.Cline{L1}.YData(end), obj.Cline{L1}.ZData(end) ;...
                obj.Cline{L2}.XData(end),obj.Cline{L2}.YData(end), obj.Cline{L2}.ZData(end)] ;
            switch src.Value
                case 1  %'H1-H2'
                    XYZ = [LineHead(1,:) ; LineHead(2,:) ]  ;
                case 2  %'H1-T2'
                    XYZ = [LineHead(1,:) ; LineTail(2,:) ]  ;
                case 3  %'T1-H2'
                    XYZ = [LineTail(1,:) ; LineHead(2,:) ]  ;
                case 4 %'T1-T2'
                    XYZ = [LineTail(1,:) ; LineTail(2,:) ]  ;
            end
            if ~isfield(d_subUI.UserData,'AlreadyPlot' )
                d_subUI.UserData.ShowLine = plot3(findobj(obj.fH,'type','Axes') , XYZ(:,1),XYZ(:,2),XYZ(:,3) ,'-k','LineWidth', 2 );
                d_subUI.UserData.AlreadyPlot =1 ;
            else
                d_subUI.UserData.ShowLine.XData= XYZ(:,1); d_subUI.UserData.ShowLine.YData= XYZ(:,2); d_subUI.UserData.ShowLine.ZData= XYZ(:,3);
            end
        end
        
        function Connect(obj,src,evn,d,popup1,popup2,popupChoice)
            L1 = popup1.Value;           L2 = popup2.Value;
            Choice = popupChoice.Value ;
            switch Choice
                case 1  %'H1-H2'
                    XYZ = [flip(obj.PointOnCurve{L1} ,1) ; obj.PointOnCurve{L2} ]  ;
                case 2  %'H1-T2'
                    XYZ = [flip(obj.PointOnCurve{L1}) ; flip(obj.PointOnCurve{L2} ,1) ]  ;
                case 3  %'T1-H2'
                    XYZ = [obj.PointOnCurve{L1} ; obj.PointOnCurve{L2} ]  ;
                case 4 %'T1-T2'
                    XYZ = [obj.PointOnCurve{L1} ; flip(obj.PointOnCurve{L2} ,1) ]  ;
            end
%             Choice
%             XYZ 
            
            obj.PointOnCurve{L1} = XYZ ;
            obj.SplineCurve{L1} =  cscvn( obj.PointOnCurve{L1}(1:end,: )') ;
            obj.PointOnCurve{L2} =[];
            obj.SplineCurve{L2} =[];
            delete(obj.Cline{L2}) ;
            delete(obj.Cscatter{L2}) ;
            
            IndRR =~cellfun('isempty',obj.PointOnCurve) ;
            obj.PointOnCurve=obj.PointOnCurve(~cellfun('isempty',obj.PointOnCurve)) ;
            obj.SplineCurve=obj.SplineCurve(~cellfun('isempty',obj.SplineCurve))  ;
            obj.CurveParas.p= obj.CurveParas.p(IndRR) ;
            obj.CurveParas.t= obj.CurveParas.t(IndRR) ;
            obj.Cline= obj.Cline(IndRR) ;
            obj.Cscatter= obj.Cscatter(IndRR) ;
            obj.VisualizeFFC ;
            
            close(d) ;
            %---Ref
            %             if ~isempty(Node)
            %                 TailNode = Node(2:end)-1 ; TailNode=[TailNode , size(XYZ,1)] ;
            %                 for k =1: length(Node)
            %                     obj.PointOnCurve{k} = XYZ(Node(k):TailNode(k) ,:  ) ;
            %                     obj.SplineCurve{k} =  cscvn( obj.PointOnCurve{k}(1:end,: )') ;
            %                 end
            %             else
            %                 obj.PointOnCurve{1} = XYZ ;
            %                 obj.SplineCurve{1} =  cscvn( obj.PointOnCurve{obj.ls}(1:end,: )') ;
            %             end
            %             obj.VisualizeFFC ;
            %-----
        end
        
        function UpdateData(obj)
            %             length( obj.Cscatter)
            RemoveEmpty= cellfun(@isempty, obj.Cscatter) ;
            obj.Cscatter=  obj.Cscatter(~RemoveEmpty);
            
            for upi = 1: length( obj.Cscatter)
                if ~isgraphics(obj.Cscatter{upi})
                    continue
                end
                
                obj.PointOnCurve{upi} = [  obj.Cscatter{upi}.XData' ,obj.Cscatter{upi}.YData',obj.Cscatter{upi}.ZData' ] ;
                obj.SplineCurve{upi} =  cscvn( obj.PointOnCurve{upi}(1:end,: )') ;
                [obj.CurveParas.p{upi},obj.CurveParas.t{upi}] = fnplt(obj.SplineCurve{upi}) ;
                if isgraphics(obj.Cline{upi})
                    obj.Cline{upi}.XData = obj.CurveParas.p{upi}(1,:);
                    obj.Cline{upi}.YData = obj.CurveParas.p{upi}(2,:);
                    obj.Cline{upi}.ZData = obj.CurveParas.p{upi}(3,:);
                    obj.Cline{upi}.UserData.splinePlotXYZ=obj.CurveParas.p{upi} ;
                    obj.Cline{upi}.LineWidth =3;
                    %                     sdfsdf=2
                end
                if isgraphics(obj.Cscatter{upi})
                    obj.Cscatter{upi}.XData = obj.PointOnCurve{upi}(:,1);
                    obj.Cscatter{upi}.YData = obj.PointOnCurve{upi}(:,2);
                    obj.Cscatter{upi}.ZData = obj.PointOnCurve{upi}(:,3);
                    obj.Cscatter{upi}.CData= repmat( [0 0.45 0.75 ],length(obj.PointOnCurve{upi}(:,1)),1  );
                    obj.Cscatter{upi}.SizeData =96 ;
                    uistack(  obj.Cscatter{upi},'top') ;
                end
                
                if isfield(obj.Cscatter{upi}.UserData ,'Ind')
                    if ~isempty( obj.Cscatter{upi}.UserData.Ind)
                        obj.Cscatter{upi}.CData(obj.Cscatter{upi}.UserData.Ind,:) = repmat([1,0,0],length(obj.Cscatter{upi}.UserData.Ind) ,1)  ;
                    end
                end
                
            end
            obj.SwitchRepresentation(obj.check_SWsplineStraight ,[] ) ;
            
            Str={'--'};
            for k=1:length( obj.Cscatter) ;  Str{k}=strcat('Spline : ',num2str(k));    end

            obj.popupH2.String= Str ;
            h= findobj(obj.fH,'type','Line','LineWidth',2) ;
            if ~isempty(h); delete(h) ;end
        end
        
        function SaveXYZ(obj,src,evn)
            [A,~] = cellfun(@size,obj.PointOnCurve) ;
            HeadNodes= cumsum([1 ,A] ) ; HeadNodes=HeadNodes(1:end-1) ;
            
            prompt = {'Enter File name:'};  dlg_title = 'Input';
            num_lines = 1;   defaultans = {'XYZ_1'};
            file3_name = inputdlg(prompt,dlg_title,num_lines,defaultans);
            
            fileID2 = fopen(strcat(file3_name{1},'.xyz')  ,'w') ;
            fprintf(fileID2,'HeadNode %s\n'  ,  num2str(HeadNodes) );  % header
            
            %             XYZ = obj.PointOnCurve{obj.ls} ;
            XYZ = zeros(sum(A),3); c=1;
            for k =1:length(obj.PointOnCurve)
                XYZ( c:A(k)+c-1 ,: )   = obj.PointOnCurve{k} ; c=c+A(k);
            end
            
            for k = 1:size(XYZ ,1)
                fprintf(fileID2,'%6.2f %6.2f %6.2f \n'  , XYZ(k,1),XYZ(k,2),XYZ(k,3) );  % header
            end
            fclose(fileID2);
        end
        
        function LoadXYZ(obj,src,evn,varargin)
            if nargin <4
                [FileName,PathName] = uigetfile('*.txt;*.xyz','Select the file with XYZ points');
            else
                FileName='Hi.xyz'   ; PathName=[pwd filesep];
            end
            
            fffid=fopen(strcat(PathName,FileName));
            Oneline = fgetl(fffid) ;
            Node= [];
            WithHeader= 0;
            if  contains(Oneline,'HeadNode' )
                Node = str2num(Oneline(9:end) ) ;
                WithHeader=1 ;
            else  % no head, pure XYZ, reset reading
                fclose(fffid);
                fffid=fopen(strcat(PathName,FileName));
            end
            
            %              sdfsf=3 ;
            %------
            count=0;
            while 1
                Oneline = fgetl(fffid) ;
                count=count+1 ;
                if isempty(Oneline)  % ||strcmp(Oneline, '-1')  ||count>20
                    break
                end
                if isa(Oneline,'double')
                    if Oneline== -1
                        break
                    end
                end
            end
            fclose(fffid);
            %-------------------------
            XYZ = zeros( count-1,3) ;
            fffid2=fopen(strcat(PathName,FileName)); count=0 ;
            if WithHeader==1
                Oneline = fgetl(fffid2) ;
            end
            
            while 1
                Oneline = fgetl(fffid2) ;
                count=count+1 ;
                if ~isa(Oneline,'double')
                    XYZ(count ,: ) = str2num( Oneline) ;
                end
                if isempty(Oneline)  % ||strcmp(Oneline, '-1')  ||count>20
                    break
                end
                if isa(Oneline,'double')
                    if Oneline== -1
                        break
                    end
                end
            end
            fclose(fffid2);
            %             cc=1;
            obj.PointOnCurve =[];
            obj.SplineCurve=[] ; obj.ls=1 ;
            obj.CurveParas=[];
            obj.Cline =[];
            obj.Cscatter =[];
            
            
            if ~isempty(Node)
                TailNode = Node(2:end)-1 ; TailNode=[TailNode , size(XYZ,1)] ;
                for k =1: length(Node)
                    obj.PointOnCurve{k} = XYZ(Node(k):TailNode(k) ,:  ) ;
                    obj.SplineCurve{k} =  cscvn( obj.PointOnCurve{k}(1:end,: )') ;
                end
            else
                obj.PointOnCurve{1} = XYZ ;
                obj.SplineCurve{1} =  cscvn( obj.PointOnCurve{obj.ls}(1:end,: )') ;
            end
            obj.VisualizeFFC ;
        end
        
        function ChangeAxesLimit(obj,src,evn ,sH)
            ax=gca;
%             evn
            interval =5; R=1.1;
%             dfsf=2
%             switch evn.Key
%                 case 'rightarrow'
% %                     sH.UserData
%                     sH
%                     [qq,indRed ]= ismember([1 0 0],sH.CData,'rows');
%                     if qq==0 ;      return;       end
%                     if indRed==1 || indRed==length(sH.XData) ;    return;       end
%                    
%                     Intv = 0.01*(max( sH.UserData.tarrs(indRed))-min( sH.UserData.tarrs(indRed))) ;
%                     t_current = sH.UserData.tarrs(indRed) ;
%                     t_new= t_current+Intv ; 
%                     NewXYZ = fnval( sH.UserData.Scale_C , t_new) ;
%                     sH.XData(indRed)=NewXYZ(1); sH.YData(indRed)=NewXYZ(2); sH.ZData(indRed)=NewXYZ(3); 
% 
%                      return
%                 case 'leftarrow'    
%                     sH
%                      [qq,indRed ]= ismember([1 0 0],sH.CData,'rows');
%                     if qq==0 ;      return;       end
%                     if indRed==1 || indRed==length(sH.XData) ;    return;       end
%                    
%                     Intv = 0.01*(max( sH.UserData.tarrs(indRed))-min( sH.UserData.tarrs(indRed))) ;
%                     t_current = sH.UserData.tarrs(indRed) ;
%                     t_new= t_current-Intv; 
%                     NewXYZ = fnval( sH.UserData.Scale_C , t_new) ;
%                     sH.XData(indRed)=NewXYZ(1); sH.YData(indRed)=NewXYZ(2); sH.ZData(indRed)=NewXYZ(3); 
% 
%                     return
%             end
            
            
            
            switch  evn.Character
                case 'q'
                    ax.XLim= ax.XLim + interval ;
                case 'Q'
                    ax.XLim= R*(ax.XLim -mean(ax.XLim))+ mean(ax.XLim)   ;
                case 'a'
                    ax.XLim= ax.XLim - interval ;
                case 'A'
                    ax.XLim= (ax.XLim -mean(ax.XLim))/R + mean(ax.XLim)   ;
                    %-----------
                case 'w'
                    ax.YLim= ax.YLim + interval ;
                case 's'
                    ax.YLim= ax.YLim - interval ;
                case 'W'
                    ax.YLim= R*(ax.YLim -mean(ax.YLim))+ mean(ax.YLim)   ;
                case 'S'
                    ax.YLim= (ax.YLim -mean(ax.YLim))/R + mean(ax.YLim)   ;
                    %------------
                case 'e'
                    ax.ZLim= ax.ZLim + interval ;
                case 'd'
                    ax.ZLim= ax.ZLim - interval ;
                case 'E'
                    ax.ZLim= R*(ax.ZLim -mean(ax.ZLim))+ mean(ax.ZLim)   ;
                case 'D'
                    ax.ZLim= (ax.ZLim -mean(ax.ZLim))/R + mean(ax.ZLim)   ;
                    
                case 'x'
                    axis equal;
                    fprintf('set axis equal \n')
                case 'X'
                    axis auto
                    fprintf('set axis auto \n')
                case 'n'
                    axis normal
                    fprintf('set axis normal \n')
                case 'N'
                    axis normal
                    fprintf('set axis normal \n')
            end
            
        end
        function SwitchRepresentation(obj,src,~ ) 
            if src.Value==1 %spline
                for k = 1:length(obj.Cline)
                   obj.Cline{k}.XData = obj.Cline{k}.UserData.splinePlotXYZ(1,:) ;
                   obj.Cline{k}.YData = obj.Cline{k}.UserData.splinePlotXYZ(2,:) ;
                   obj.Cline{k}.ZData = obj.Cline{k}.UserData.splinePlotXYZ(3,:) ;
                end
            else  % straight-line
                for k = 1:length(obj.Cline)
                     obj.Cline{k}.XData = obj.Cscatter{k}.XData ;
                     obj.Cline{k}.YData = obj.Cscatter{k}.YData ;
                     obj.Cline{k}.ZData = obj.Cscatter{k}.ZData ;
                    
                end
            end
            
            
%             dssf=3
%             
%             sdfs=3
        end
        
        
    end
    
end

