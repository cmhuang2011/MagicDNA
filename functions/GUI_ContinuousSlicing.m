function Cont_FFC_data = GUI_ContinuousSlicing( FFC,InplaneXY,VecAll, L_bundle,XsecLattice , state,Prev_Newtarr ,AdjLines)
% This is the GUI for assign details about ss connetion and tuning the
% orientation
%-


if isempty(state)
    Cont_FFC_data=[];
    % if ~isempty(Prev_ssDNA) ; output= Prev_ssDNA ;   end
    
    d = figure('Position', get( groot, 'Screensize' ),'Name','Slicing  options');
    d.Units='normalized';d.OuterPosition=[0 0 1 1];  clf; hold on ;
    aH=gca;
    aH.Position=[0.1,0.1,0.75,0.8]; axis equal;
    
    
    % str1='Length (unit: nt) of extended scaffold ends to form loops for preventing base stacking ';
    % txt1 = uicontrol('Parent',d,'Style','text','Units','normalized',...
    %     'Position',[0.85 0.9 0.1 0.05], 'String',str1 ,'Enable','Inactive','HorizontalAlignment','left');
    
    str = cellfun(@num2str,num2cell([10 20 30 40 50 60 70]),'uniformoutput',0);  % single strand options
    popup1 = uicontrol('Parent',d , 'Style','popup','Units','normalized',...
        'Position',[0.9 0.85 0.08 0.05], 'String',str);
    popup1.Callback= @(src,evn)TargetNSlice(src,evn, FFC,InplaneXY,VecAll, L_bundle,XsecLattice  ,AdjLines) ;
    
    btn1 = uicontrol('Parent',d,'Units','normalized', 'Position',[0.9 0.05 0.08 0.1],   'String','OK');
    % btn1.Position=[0.88 0.05 0.09 0.1] ;  % ok
    % btn1.Callback=@(src,evn)Exportresult(src,evn,t,popup1,ppH) ;
    btn2 = uicontrol('Parent',d,'Units','normalized', 'Position',[0.9 0.25 0.08 0.1],   'String','Preview');
    
    
    AssignIcon(  btn1,'Sweep_ok.jpg' ) ;  btn1.TooltipString='Export to Assembly.' ;
    AssignIcon(  btn2,'Sweep_view.jpg' ) ;btn2.TooltipString='Check Cylinder lengths in each segment.' ;
    
    Rough_split_ForEachSlice = 10 ;  %
    d.UserData.popup1 =popup1;
else
    d= gcf ;
    Cont_FFC_data=[];
    Rough_split_ForEachSlice =str2num(d.UserData.popup1.String{d.UserData.popup1.Value}) ;  %
    
    cla ;
end
%------------------------
InplaneXY ;
nCylinder =size(InplaneXY ,1) ; 
Bundle1Length_in_nm = L_bundle(1) *0.34 ;
Ori_L_curve = norm( diff(FFC.PointOnCurve{1}(1:2,:)) ) ;
% AdjLines
% size(AdjLines)
% size(FFC.PointOnCurve{1} ,1)
if AdjLines(1, size(AdjLines ,2)) == 0 % open curve, Feb 28 2021
    
    FFC.PointOnScaleCurve = FFC.PointOnCurve{1}*Bundle1Length_in_nm/Ori_L_curve ;
    nBundle =size(FFC.PointOnCurve{1},1)-1 ;
    LoopCase =false ;
else
    
    FFC.PointOnScaleCurve = [FFC.PointOnCurve{1} ;FFC.PointOnCurve{1}(1,:) ] *Bundle1Length_in_nm/Ori_L_curve ;
    nBundle =size(FFC.PointOnCurve{1},1) ;
    LoopCase =true ;
end


POSC =FFC.PointOnScaleCurve;
FFC.Scaled_Curve =    cscvn(  FFC.PointOnScaleCurve') ; Scale_C =FFC.Scaled_Curve ;

sdsfs=3;

%             ax= gca ; ax.Position = [0.05 0.05 0.75 0.9 ] ;
for k = 1 :1 %length(FFC.PointOnCurve)
    AA = scatter3(POSC(:,1),POSC(:,2) ,POSC(:,3),58,'filled'  );
    AA.CData= repmat( [0 0.45 0.75 ],length(POSC(:,1)),1  );
    AA.HitTest='off' ;AA.Visible='off';
    [WW.p,FFC.CurveParas.t] = fnplt(Scale_C) ;
    BB=  plot3( WW.p(1,:),WW.p(2,:),WW.p(3,:),'LineWidth',3 );
    BB.HitTest='off' ;
end
aH=gca;
%---------
N_slice= 20*nBundle+1 ;
%  N_slice =1000 ;
% 2n+1 , or N_slice =1000
% number of slices before finer slicing. For robot-like multi-bundles, the
% orientation are not so continuous. Slicing more may have different end
% orietations of bundles. If this number is not small, when the local
% coordinate is not parallel to tangential direction, unfeasible twisting
% may happen.

InitOrt_X = VecAll(1,7:9)';  InitOrt_X=InitOrt_X/norm(InitOrt_X) ;

% initial orientation of x from sketch
% t_arrs= sort( [Scale_C.breaks ,movmean(Scale_C.breaks , 2,'Endpoints','discard')] )  ;
t_arrs=linspace(Scale_C.breaks(1),Scale_C.breaks(end),N_slice ) ;
Position_TanVec_NormalVec = zeros(9,N_slice ) ;
for  pointi = 1: N_slice
    Position_TanVec_NormalVec(1:3, pointi) = fnval(Scale_C, t_arrs(pointi))  ;% position
end

QQ= gradient(Position_TanVec_NormalVec(1:3,: ) ) ; % tangent direction
Position_TanVec_NormalVec(4:6,:)=  vec3norm_CM( QQ' )' ;
InitOrt_X= cross(Position_TanVec_NormalVec(4:6,1)   ,cross(InitOrt_X,Position_TanVec_NormalVec(4:6,1))) ;

if LoopCase
%     MeanTan = mean(Position_TanVec_NormalVec(4:6,[1 end])' ) ;
%     MeanTan = vec3norm_CM( MeanTan ) ;    
%     Position_TanVec_NormalVec(4:6 ,1) = MeanTan' ;
%     Position_TanVec_NormalVec(4:6 ,end) = MeanTan' ;   
%     Position_TanVec_NormalVec(:,end) =   Position_TanVec_NormalVec(:,1) ;
    Position_TanVec_NormalVec(:,end+1) =   Position_TanVec_NormalVec(:,2) ;    
end


Position_TanVec_NormalVec(7:9,1) =InitOrt_X ;
%------Make continuous
%             Formers =QQ2(:,1:end-1 ); Laters =QQ2(:,2:end );
for k =2:    size(Position_TanVec_NormalVec,2)
    NVec =  Position_TanVec_NormalVec(7:9,k-1) ;
    CurTVec = Position_TanVec_NormalVec(4:6,k) ;
    PrevTVec = Position_TanVec_NormalVec(4:6,k-1) ;
    if dot(CurTVec,PrevTVec)>0  % angle smaller than 90
        Ytemp= cross(CurTVec, NVec) ;Ytemp=Ytemp/norm(Ytemp) ;
    else
        Ytemp= cross(-CurTVec, NVec) ;Ytemp=Ytemp/norm(Ytemp)  ;
    end
    New_NVec = cross(Ytemp,Position_TanVec_NormalVec(4:6,k)) ;
    
    %     dT =Position_TanVec_NormalVec(4:6,k)-Position_TanVec_NormalVec(4:6,k-1);
    %     New_NVec= dT/norm(dT) ;
    
    if dot(New_NVec, Position_TanVec_NormalVec(7:9,k-1))<0
        New_NVec=-New_NVec ;
    end
    
    Position_TanVec_NormalVec(7:9,k) = New_NVec ;
end

if LoopCase
    Position_TanVec_NormalVec(4:6 ,end-1 ) = Position_TanVec_NormalVec(4:6 ,1 ) ; % force tangent equalat head and tail
    Position_TanVec_NormalVec=Position_TanVec_NormalVec(:,1:end-1) ;
    Zdir = Position_TanVec_NormalVec(4:6 ,1) ;
    X1 = Position_TanVec_NormalVec(7:9 ,1) ;
    X2 = Position_TanVec_NormalVec(7:9 ,end) ;    
    Y1 = cross( Zdir , X1) ;
    
    tan2x_onfirst = dot(X2 ,  X1 )  ;
    tan2y_onfirst = dot(X2 ,  Y1 )  ;
    closedLoop_twist = atan2d(tan2y_onfirst,tan2x_onfirst) ;
    IndivTwist =  closedLoop_twist/(size(Position_TanVec_NormalVec,2)-1 )  ;
    CumTwistAng = IndivTwist ;
    for k =2 : size(Position_TanVec_NormalVec,2)
        RMat_twist = RotationAxisToRotaionMatrix(   Position_TanVec_NormalVec( 4:6,k), -CumTwistAng )   ;
        Position_TanVec_NormalVec( 7:9,k) = RMat_twist* Position_TanVec_NormalVec( 7:9,k) ;
        CumTwistAng=CumTwistAng+IndivTwist ;
    end
end

%------------local orientation of each slice
InplaneXY=InplaneXY-mean(InplaneXY) ;
InplaneXY=InplaneXY*1.25 ;

S_PTNC =  Position_TanVec_NormalVec ;
N_FineSlice = round(2* sum(L_bundle) ) ; % Don't need to be too dense
Fine_tarr = linspace(t_arrs(1),t_arrs(end), N_FineSlice) ;
PTN_vec_onSpline = interp1(t_arrs,S_PTNC', Fine_tarr  ,'PCHIP'  ) ;

% SumFFC_unitBPS = sum(L_bundle) ;

PTN_vec_onSpline(:,1:3) = fnval(Scale_C, Fine_tarr)'  ;% position
QQ = fnval( fnder(Scale_C), Fine_tarr)'  ; % tangent direction
PTN_vec_onSpline(:,4:6)= vec3norm_CM( QQ ) ;   %Tan Vec

% PTN_vec_onSpline([1 end],:)

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
    Cyl_h(k) =plot3(AllXYZ(3*k-2 ,:),AllXYZ(3*k-1 ,:),AllXYZ(3*k ,:)  ,'-k','LineWidth',0.5) ;
    Cyl_h(k).ButtonDownFcn = @(src,evn)ShowInd(FFC,src,evn) ;
    Cyl_h(k).HitTest='off' ;
end
%----------
% set(Cyl_h,'Visible','off') ;
% set(Cyl_h(end),'Visible','on') ;
% set(Cyl_h(1),'Visible','on') ;
%--------Visualize

simplifyLine_byGtheory(Cyl_h) ;

Inds = 1:10:N_FineSlice ;  % Inds = 1:1:N_slice ;
Ez =PTN_vec_onSpline' ;   %  Ez =Position_TanVec_NormalVec ;

% q1 = quiver3(Ez(1,Inds),Ez(2,Inds),Ez(3,Inds),Ez(4,Inds),Ez(5,Inds),Ez(6,Inds) ,'Color',[0.9290 0.6940 0.1250] ) ;
% q2=  quiver3(Ez(1,Inds),Ez(2,Inds),Ez(3,Inds),Ez(7,Inds),Ez(8,Inds),Ez(9,Inds),'Color',[0.4940 0.1840 0.5560] ) ;
% q1.HitTest='off' ;    q2.HitTest='off' ;
% save 'SweepData.mat' Ez;


axis equal ;
xlabel('X'); ylabel('Y') ;zlabel('Z');
% set(figure(2426),'KeyPressFcn',@(src,evn)ChangeAxesLimit(FFC,src,evn) )  ;
%             return
%----------------decide how many bundles on this curve
sH=scatter3(AllXYZ(1,1),AllXYZ(1,2) ,AllXYZ(1,3),96 ,'filled'  );    sH.CData=[0 1 1];%initialize and some point
sH.CData=[0 1 1];
sH.ButtonDownFcn= @(src,evn)ClickSH(src,evn) ;


[New_tarr, NewNodes]=SlicingCoarseCut(Rough_split_ForEachSlice,Scale_C ,N_slice ,t_arrs) ;
if ~isempty(Prev_Newtarr) && isempty(state)
    Prev_NewNodes   = fnval(Scale_C , Prev_Newtarr ) ;
    NewNodes= Prev_NewNodes' ;
    New_tarr=Prev_Newtarr ;
end
% ssdf=4
sH.XData= NewNodes(:,1)' ; sH.YData= NewNodes(:,2)' ; sH.ZData= NewNodes(:,3)' ;  % update graphic
sH.UserData.New_tarr= New_tarr ; sH.UserData.Scale_C=Scale_C ;

[DotOnCylinders,ds_Cyli_ForBundle ,All_Base_onCyl_Pos]=SlicingCoarseCut_part2(nCylinder,sH,Cyl_h,LoopCase) ;

%----------------------------
AllJoints = reshape(DotOnCylinders,size(DotOnCylinders,1)*size(DotOnCylinders,2),3) ;
sH_alljoint = scatter3(aH,AllJoints(:,1),AllJoints(:,2),AllJoints(:,3),'r.') ;
sH.UserData.ds_Cyli_ForBundle =ds_Cyli_ForBundle ;


%--------------------------------------------------------------
%            % Exprot Hyperbundle ------------------
%           Cont_FFC_data = [] ;
%           Cont_FFC_data
Center_t =  0.5*(New_tarr(1:end-1) +New_tarr(2:end));
Bundle.Center =   0.5*( NewNodes(1:end-1,:)+ NewNodes(2:end,:) ) ;
Bundle.TanVec =NewNodes(2:end,:)-NewNodes(1:end-1,:) ;
Bundle.NVec= interp1(Fine_tarr,PTN_vec_onSpline(:,7:9), Center_t  ,'PCHIP'  ) ;
Bundle.New_tarr=New_tarr;

Cylinders.N = nCylinder ;
Cylinders.Nbundles = size(NewNodes,1)-1 ;
Cylinders.DotOnCylinders=DotOnCylinders ;
Cylinders.ds_Cyli_ForBundle=ds_Cyli_ForBundle;
Cylinders.All_Base_onCyl_Pos =All_Base_onCyl_Pos ;

FFC.t_arrs_fine = Fine_tarr;
FFC.PTN_Vec_onCurve = PTN_vec_onSpline;

Cont_FFC_data.Bundle=Bundle ;
Cont_FFC_data.Cylinders=Cylinders ;

% figure(325);clf;
% subplot(3,1,1) ; hold on ;
% NBaseOncyl= zeros(nCylinder,1) ;
% for k =  1: nCylinder
%     XYZ = [ Cyl_h(k).XData ; Cyl_h(k).YData ; Cyl_h(k).ZData ]' ;
%     dx = diff(XYZ) ;  qq = sqrt(sum(dx.^2 ,2)) ;  plot(qq)
%     NBaseOncyl(k)=length(Cyl_h(k).XData);
% end
% title('Base to base distance on cylinders.'); xlabel('base');ylabel('distance (nm)')
% subplot(3,1,2) ; hold on ;
% for k =  1: nCylinder
%     plot(ds_Cyli_ForBundle(k,:) ,'--');
% end
% plot(mean(ds_Cyli_ForBundle ) ,'LineWidth',2)
% title('Cylinder lengths of each section, Bold for neutral.'); xlabel('bundle');ylabel('Length (nm)');
%
% subplot(3,1,3) ; hold on ;
% bar(NBaseOncyl);
% title('Entire lengths of cylinders'); xlabel('cylinder');ylabel('N (bps)')
% drawnow ;

d.UserData.sH=sH ;
d.UserData.Cyl_h=Cyl_h ;

d.WindowKeyPressFcn=@(src,evn)ChangeAxesLimitLocal(src,evn,sH,sH_alljoint,nCylinder,Cyl_h ,LoopCase) ;
btn2.Callback=@(src,evn)Preview(src,evn,nCylinder,Cyl_h,sH) ;
btn1.Callback=[];
btn1.Callback=@(src,evn)Exportresult(src,evn,nCylinder,Cyl_h,sH,sH_alljoint,PTN_vec_onSpline,Fine_tarr,FFC) ;
% sH
uistack(sH,'top') ;
%---------------

uiwait(d);


    function Exportresult(src,evn,nCylinder,~,~,sH_alljoint,PTN_vec_onSpline,Fine_tarr,FFC)
        d=gcf;
        
        sH=d.UserData.sH;
        Cyl_h=d.UserData.Cyl_h;
        New_tarrtt=sH.UserData.New_tarr ;
        NewNodestt= [sH.XData ; sH.YData ;sH.ZData ]' ;
        
        
        Center_ttt =  0.5*(New_tarrtt(1:end-1) +New_tarrtt(2:end));
        Bundlett.Center =   0.5*( NewNodestt(1:end-1,:)+ NewNodestt(2:end,:) ) ;
        Bundlett.TanVec =NewNodestt(2:end,:)-NewNodestt(1:end-1,:) ;
        Bundlett.NVec= interp1(Fine_tarr,PTN_vec_onSpline(:,7:9), Center_ttt  ,'PCHIP'  ) ;
        Bundlett.New_tarr=New_tarrtt ;
        
        [DotOnCylinderstt,ds_Cyli_ForBundlett,All_Base_onCyl_Postt ]=SlicingCoarseCut_part2(nCylinder,sH,Cyl_h ,LoopCase) ;
        
        Cylinderstt.N = nCylinder ;
        Cylinderstt.Nbundles = size(NewNodestt,1)-1 ;
        Cylinderstt.DotOnCylinders=DotOnCylinderstt ;
        Cylinderstt.ds_Cyli_ForBundle=ds_Cyli_ForBundlett;
        Cylinderstt.All_Base_onCyl_Pos =All_Base_onCyl_Postt ;
        
        
        Cont_FFC_data.Bundle=Bundlett ;
        Cont_FFC_data.Cylinders=Cylinderstt ;
        
        FFC.PointOnScaleCurve_afterSlicing= NewNodestt ;
        FFC.Bundlett =Bundlett;
        FFC.Cylinderstt =Cylinderstt;
        
        
        delete(gcf)
    end % end fcn Exportresult






end
% sdf=3
%-------------------------------------------------------End main
function Preview(src,evn,nCylinder,~,~)
d=gcf;
sH=d.UserData.sH;
Cyl_h=d.UserData.Cyl_h;



figure(3251);clf;
subplot(3,1,1) ; hold on ;
NBaseOncyl= zeros(nCylinder,1) ;

SaveXYZ = cell(nCylinder ,1) ; % use for transformed routing
for k =  1: nCylinder
    XYZ = [ Cyl_h(k).XData ; Cyl_h(k).YData ; Cyl_h(k).ZData ]' ;
    dx = diff(XYZ) ;  qq = sqrt(sum(dx.^2 ,2)) ;  plot(qq)
    NBaseOncyl(k)=length(Cyl_h(k).XData);
    
    SaveXYZ{k} = XYZ ;    
end
% save ExtrudedXYZ.mat SaveXYZ

title('Base to base distance on cylinders.'); xlabel('base');ylabel('distance (nm)')
subplot(3,1,2) ; hold on ;
for k =  1: nCylinder
    plot(sH.UserData.ds_Cyli_ForBundle(k,:) ,'--');
end
plot(mean(sH.UserData.ds_Cyli_ForBundle ) ,'LineWidth',2) ;
plot([0,size(sH.UserData.ds_Cyli_ForBundle,2) ],[10 10 ],'--k' ) ;
title('Cylinder lengths of each section, Bold for neutral.'); xlabel('bundle');ylabel('Length (nm)');

subplot(3,1,3) ; hold on ;
bar(NBaseOncyl);
title('Entire lengths of cylinders'); xlabel('cylinder');ylabel('N (bps)')
drawnow ;


end

function TargetNSlice(src,evn, FFC,InplaneXY,VecAll, L_bundle,XsecLattice ,AdjLines)

GUI_ContinuousSlicing( FFC,InplaneXY,VecAll, L_bundle,XsecLattice , 'Repeat',[] ,AdjLines) ;

end


function ClickSH(src,evn)
%             size(src.CData)
evn.Button;
if  evn.Button ==1
    src.CData = repmat([0 1 1] ,length(src.XData),1 )  ;
    
    XYZ = [src.XData;src.YData;src.ZData ]' ;
    dxyz =  XYZ-  repmat(evn.IntersectionPoint,length(src.XData),1) ;
    d = sqrt( sum(dxyz.^2 ,2) ) ;
    ind = find(d==min(d)) ; ind=ind(1) ;
    src.CData(ind,:) = [1 0 0] ;
    % elseif evn.Button ==1
    %     [qq,indRed ]= ismember([1 0 0],src.CData,'rows')
    %     if qq==0 ;      return;       end
    %     if indRed==1 || indRed==length(src.XData) ;    return;       end
    %
    %     Intv = 0.05*(max( src.UserData.tarrs)-min( src.UserData.tarrs)) ;
    %     t_current = src.UserData.tarrs(indRed) ;
    %     t_new= t_current-Intv;
    %     NewXYZ = fnval( src.UserData.Scale_C , t_new)
    %     src.XData(indRed)=NewXYZ(1)   ;
    %     src.YData(indRed)=NewXYZ(2)   ;
    %     src.ZData(indRed)=NewXYZ(3)   ;
    %     src.UserData.tarrs(indRed) =t_new ;
    % elseif evn.Button ==3
    %     [qq,indRed ]= ismember([1 0 0],src.CData,'rows')
    %     if qq==0 ;      return;       end
    %     if indRed==1 || indRed==length(src.XData) ;    return;       end
    %
    %     Intv = 0.05*(max( src.UserData.tarrs(indRed))-min( src.UserData.tarrs(indRed))) ;
    %     t_current = src.UserData.tarrs(indRed) ;
    %     t_new= t_current+Intv;
    %     NewXYZ = fnval( src.UserData.Scale_C , t_new)
    %     src.XData(indRed)=NewXYZ(1)   ;
    %     src.YData(indRed)=NewXYZ(2)   ;
    %     src.ZData(indRed)=NewXYZ(3)   ;
    %     src.UserData.tarrs(indRed) =t_new ;
    %
    
end

end

function  [New_tarr, NewNodes]=SlicingCoarseCut(Total_split_ForEachSlice,Scale_C ,N_slice ,t_arrs)

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

if Total_split_ForEachSlice<N_slice
    [aa,bb] =sort(AvgCurvature,'descend' ) ;
    %     Nb =zeros(size(Nb) ) ;
    %     New_tarr= sort( [t_arrs(1), t_arrs((bb(1:Total_split_ForEachSlice))), t_arrs(end)] ) ;
    
    New_tarr= linspace(t_arrs(1), t_arrs(end) ,Total_split_ForEachSlice  );
%     sdfs=3
    
end

NewNodes = fnval(Scale_C , New_tarr ) ;
NewNodes=NewNodes';  % use this to detect closest point on cylinders
%             sH=scatter3(NewNodes(:,1),NewNodes(:,2) ,NewNodes(:,3),64 ,'filled'  );
%             sH.CData=[0 1 1];
%-----------

end


function  [DotOnCylinders,ds_Cyli_ForBundle,All_Base_onCyl_Pos ] =SlicingCoarseCut_part2(nCylinder,sH,Cyl_h ,LoopCase)
sHXYZ =[sH.XData ;sH.YData ;sH.ZData ]' ;
% NewNodes=sHXYZ;
% DotOnCylinders = zeros(nCylinder , size(NewNodes,1) ,3 ) ;
if LoopCase % closed-loop case
NewNodes=sHXYZ;
DotOnCylinders = zeros(nCylinder , size(NewNodes,1) ,3 ) ;   
ds_Cyli_ForBundle=  zeros(nCylinder , size(NewNodes,1)-1  ) ;  % unit: nm
else % open-chain case
NewNodes=sHXYZ;
DotOnCylinders = zeros(nCylinder , size(NewNodes,1) ,3 ) ;  
ds_Cyli_ForBundle=  zeros(nCylinder , size(NewNodes,1)-1  ) ;  % unit: nm
end

All_Base_onCyl_Pos = cell(nCylinder,1) ;
for k_cyl =  1: nCylinder
    XYZ_onCyl = [Cyl_h(k_cyl).XData ;Cyl_h(k_cyl).YData; Cyl_h(k_cyl).ZData]' ;
    
%     k = dsearchn(P,PQ)
   Indk  =   dsearchn(XYZ_onCyl,NewNodes)     ;
   Indk(1)=1 ; Indk(end) =size(XYZ_onCyl,1) ;  % forced to the head and tail
   PointOnEachCyl = XYZ_onCyl(Indk ,:) ;

   DotOnCylinders(k_cyl,:,1) = PointOnEachCyl(:,1) ;
   DotOnCylinders(k_cyl,:,2) = PointOnEachCyl(:,2) ;
   DotOnCylinders(k_cyl,:,3) = PointOnEachCyl(:,3) ;
   
   for j = 1:length(Indk)-1
      ds = XYZ_onCyl(Indk(j)+1:Indk(j+1),:) -  XYZ_onCyl(Indk(j):Indk(j+1)-1,:) ;
      Total_ds= sum( sqrt(sum(ds.^2 ,2)) )  ;
     ds_Cyli_ForBundle(k_cyl, j) =Total_ds ;
   end
   
   All_Base_onCyl_Pos{k_cyl}= XYZ_onCyl ;
% sdfsdf=2
% 
%     for j_spl = 1 :size(ds_Cyli_ForBundle,2)
%         dXYZ = XYZ_onCyl- ones(size(XYZ_onCyl,1),1)*NewNodes(j_spl,:) ;
%         d =  sqrt(sum(dXYZ.^2 ,2)) ;
%         if j_spl==1
%         Ind = find(d==min(d)) ;Ind=Ind(1) ;
%         else
%         Ind = find(d==min(d)) ;Ind=max(Ind) ;
%         end
%         
%         DotOnCylinders(k_cyl, j_spl ,1) =  XYZ_onCyl(Ind,1) ;
%         DotOnCylinders(k_cyl, j_spl ,2) =  XYZ_onCyl(Ind,2) ;
%         DotOnCylinders(k_cyl, j_spl ,3) =  XYZ_onCyl(Ind,3) ;
%         if j_spl~=1  %numerical integral on the cylinder curve, instead of end-to-end distance
%             ds = XYZ_onCyl(PrevInd+1:Ind,:) -  XYZ_onCyl(PrevInd:Ind-1,:) ;
%             Total_ds= sum( sqrt(sum(ds.^2 ,2)) ) 
%             ds_Cyli_ForBundle(k_cyl, j_spl-1) =Total_ds ;
%         end
%         PrevInd = Ind ;
%     end
%     All_Base_onCyl_Pos{k_cyl}= XYZ_onCyl ;
end
end

% function  [DotOnCylinders,ds_Cyli_ForBundle,All_Base_onCyl_Pos ] =SlicingCoarseCut_part2(nCylinder,sH,Cyl_h ,LoopCase)
% sHXYZ =[sH.XData ;sH.YData ;sH.ZData ]' ;
% % NewNodes=sHXYZ;
% % DotOnCylinders = zeros(nCylinder , size(NewNodes,1) ,3 ) ;
% if LoopCase
% NewNodes=sHXYZ;
% DotOnCylinders = zeros(nCylinder , size(NewNodes,1) ,3 ) ;   
% ds_Cyli_ForBundle=  zeros(nCylinder , size(NewNodes,1)  ) ;  % unit: nm
% else
% NewNodes=sHXYZ;
% DotOnCylinders = zeros(nCylinder , size(NewNodes,1) ,3 ) ;
%    
% ds_Cyli_ForBundle=  zeros(nCylinder , size(NewNodes,1)-1  ) ;  % unit: nm
% end
% 
% All_Base_onCyl_Pos = cell(nCylinder,1) ;
% for k_cyl =  1: nCylinder
%     XYZ_onCyl = [Cyl_h(k_cyl).XData ;Cyl_h(k_cyl).YData; Cyl_h(k_cyl).ZData]' ;
%     for j_spl = 1 :size(ds_Cyli_ForBundle,2)
%         dXYZ = XYZ_onCyl- ones(size(XYZ_onCyl,1),1)*NewNodes(j_spl,:) ;
%         d =  sqrt(sum(dXYZ.^2 ,2)) ;
%         if j_spl==1
%         Ind = find(d==min(d)) ;Ind=Ind(1) ;
%         else
%         Ind = find(d==min(d)) ;Ind=max(Ind) ;
%         end
%         
%         DotOnCylinders(k_cyl, j_spl ,1) =  XYZ_onCyl(Ind,1) ;
%         DotOnCylinders(k_cyl, j_spl ,2) =  XYZ_onCyl(Ind,2) ;
%         DotOnCylinders(k_cyl, j_spl ,3) =  XYZ_onCyl(Ind,3) ;
%         if j_spl~=1  %numerical integral on the cylinder curve, instead of end-to-end distance
%             ds = XYZ_onCyl(PrevInd+1:Ind,:) -  XYZ_onCyl(PrevInd:Ind-1,:) ;
%             Total_ds= sum( sqrt(sum(ds.^2 ,2)) ) 
%             ds_Cyli_ForBundle(k_cyl, j_spl-1) =Total_ds ;
%         end
%         PrevInd = Ind ;
%     end
%     All_Base_onCyl_Pos{k_cyl}= XYZ_onCyl ;
% end
% end

function simplifyLine_byGtheory( Cyl_h)
InterVal_D = 0.34  ; % nm
for k =1: length(Cyl_h)
    XYZ = [Cyl_h(k).XData; Cyl_h(k).YData; Cyl_h(k).ZData ]';                                size(XYZ);
    Dxyz = diff(XYZ) ;   Dist= sqrt(sum(Dxyz.^2 ,2)) ;
    cumS_D = cumsum(Dist) ;
    IntegerInds = ceil(cumS_D/InterVal_D) ;
    EquPre= IntegerInds(1:end-1)~=IntegerInds(2:end) ;
    XYZ=XYZ([EquPre;false],:) ;   size(XYZ)   ;
    %--------------
    % % %                 sdsf=3
    XYZ2 =XYZ ; nNode =size(XYZ2,1) ;
    ST =zeros(nNode,2 ) ; c=1 ;
    for bi =1 : nNode
        for bj =  bi+1:nNode
            if norm(XYZ(bi,:)-XYZ(bj,:) ) < 1*InterVal_D
                %                                                                 scatter3(XYZ2([bi bj],1),XYZ2([bi bj],2),XYZ2([bi bj],3) )
                ST(c,:) =[bi,bj ] ;c=c+1 ;
            end
        end
    end
    ST=ST(1:c-1,:) ;   ST2 =[ST ;[ 1:nNode-1 ;2:nNode]'] ;
    ST2= unique(ST2,'rows') ;
    
    G = digraph(ST2(:,1),ST2(:,2)) ;
    Pind = shortestpath(G,1,nNode ,'Method','positive') ;
    XYZ= XYZ(Pind,:) ;                size(XYZ);
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


function ChangeAxesLimitLocal(src,evn,sH,sH_alljoint,nCylinder,Cyl_h  ,LoopCase)
ax=gca;
interval =5; R=1.1;
% evn
MovesNodeUse =0;
switch evn.Key
    case 'c'  % export curved cylinder model to Chimera
            fprintf('printing extrusion models for Chimera.\n')
   
% %         Cyl_h
%             prompt = {'Enter a value of scale '};
%             dlgtitle = 'Scale Value';
%             definput = {'1'};
%             dims = [1 40];
%             opts.Interpreter = 'tex';
%             answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
%             Rs = str2double(answer{1});
            Rs =1 ;
            
            file_name='Curved_extrusion'    
            Radius =1.28;
            fileID = fopen([pwd filesep file_name '.bild'],'w');
            for Bi = 1:length(Cyl_h )
                fprintf(fileID ,'\n' );
                %                         fprintf(fileID ,'.transparency 0.5\n' );
                fprintf(fileID ,'.comment (If need colors for bundle %i) .color  r g b\n',Bi );
                 XYZ_onCyl = Rs*[Cyl_h(Bi).XData ;Cyl_h(Bi).YData; Cyl_h(Bi).ZData]' ;
                 
                 % hard code for extract certain segments
%                  BP=[sH_alljoint.XData ; sH_alljoint.YData ;sH_alljoint.ZData ]' ;
%                  IndCut = ismember( XYZ_onCyl , BP ,'rows' ) ;
%                  Ind = find(IndCut) ;  Ind=Ind(1:3) ;
%                  XYZ_onCyl=XYZ_onCyl(Ind(1):Ind(3),:) ;
%                                   
                for Cj= 1:size(XYZ_onCyl,1)-1
                    fprintf(fileID , '.cylinder %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f open\n',XYZ_onCyl(Cj,:),XYZ_onCyl(Cj+1,:),Radius )    ;
                     fprintf(fileID , '.sphere %8.6f %8.6f %8.6f  %8.6f \n',XYZ_onCyl(Cj,:),Radius )    ;
                end
                  fprintf(fileID , '.sphere %8.6f %8.6f %8.6f  %8.6f \n',XYZ_onCyl(Cj+1,:),Radius )    ;
            end
            fclose(fileID);

            
%             sH_alljoint;
            file_name2='Curved_extrusion_CyleDivid'   
            Radius =1.02*Radius;
            fileID2 = fopen([pwd filesep file_name2 '.bild'],'w');
             fprintf(fileID ,'.comment (If need colors for bundle %i) .color  r g b\n',Bi );
            for k =1: length(sH_alljoint.XData)
                XYZi =Rs*[sH_alljoint.XData(k),sH_alljoint.YData(k),sH_alljoint.ZData(k) ] ;
              fprintf(fileID , '.sphere %8.6f %8.6f %8.6f  %8.6f \n',XYZi,Radius )    ;   
            end
            
             fclose(fileID2);
        
    case 'backspace'
        [qq,indRed ]= ismember([1 0 0],sH.CData,'rows') ;
        if qq==0 ;      return;       end
        if indRed==1 || indRed==length(sH.XData) ;    return;       end
        sH.XData(indRed)=[]   ;    sH.YData(indRed)=[]   ;    sH.ZData(indRed)=[]   ;
        sH.UserData.New_tarr(indRed) =[] ; sH.CData(indRed,:) =[];
        
        [DotOnCylinders,ds_Cyli_ForBundle ]=SlicingCoarseCut_part2(nCylinder,sH,Cyl_h ,LoopCase) ;
        AllJoints = reshape(DotOnCylinders,size(DotOnCylinders,1)*size(DotOnCylinders,2),3) ;
        %     sH_alljoint = scatter3(aH,AllJoints(:,1),AllJoints(:,2),AllJoints(:,3),'r.') ;
        sH_alljoint.XData =AllJoints(:,1) ;
        sH_alljoint.YData =AllJoints(:,2) ;
        sH_alljoint.ZData =AllJoints(:,3) ;
        sH.UserData.ds_Cyli_ForBundle =ds_Cyli_ForBundle ;
        %            simplifyLine_byGtheory( Cyl_h) ;
        return ;
    case 'leftarrow'
        [qq,indRed ]= ismember([1 0 0],sH.CData,'rows') ;
        if qq==0 ;      return;       end
        if indRed==1 || indRed==length(sH.XData) ;    return;       end
        Intv = 0.002*(max( sH.UserData.New_tarr)-min( sH.UserData.New_tarr)) ;
        MovesNodeUse=1 ;
    case 'rightarrow'
        [qq,indRed ]= ismember([1 0 0],sH.CData,'rows') ;
        if qq==0 ;      return;       end
        if indRed==1 || indRed==length(sH.XData) ;    return;       end
        Intv = -0.002*(max( sH.UserData.New_tarr)-min( sH.UserData.New_tarr)) ;
        MovesNodeUse=1 ;
end
if MovesNodeUse==1
    t_current = sH.UserData.New_tarr(indRed) ;
    t_new= t_current-Intv;
    if t_new<sH.UserData.New_tarr(indRed-1) || t_new>sH.UserData.New_tarr(indRed+1)
        return;
    end
    
    NewXYZ = fnval( sH.UserData.Scale_C , t_new) ;
    sH.XData(indRed)=NewXYZ(1)   ;    sH.YData(indRed)=NewXYZ(2)   ;    sH.ZData(indRed)=NewXYZ(3)   ;
    sH.UserData.New_tarr(indRed) =t_new ;
    
    [DotOnCylinders,ds_Cyli_ForBundle ]=SlicingCoarseCut_part2(nCylinder,sH,Cyl_h,LoopCase) ;
    AllJoints = reshape(DotOnCylinders,size(DotOnCylinders,1)*size(DotOnCylinders,2),3) ;
    %     sH_alljoint = scatter3(aH,AllJoints(:,1),AllJoints(:,2),AllJoints(:,3),'r.') ;
    sH_alljoint.XData =AllJoints(:,1) ;
    sH_alljoint.YData =AllJoints(:,2) ;
    sH_alljoint.ZData =AllJoints(:,3) ;
    sH.UserData.ds_Cyli_ForBundle =ds_Cyli_ForBundle ;
    %     simplifyLine_byGtheory( Cyl_h) ;
    return;
end




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
    case 'P'
        fprintf('Capturing current figure \n') ;
        print('-r100',gcf,'SnapShotFromSSDNA_r100','-djpeg')  ;
    case 'p'
        fprintf('Capturing current figure \n') ;
        print('-r100',gcf,'SnapShotFromSSDNA_r100','-djpeg')  ;
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


function CheckFcnshowssDNAL(src,evn,t,ppH)

if src.Value ==1
    for k_ss=1:length(ppH)
        ppH{k_ss}.UserData.labelL.Visible='on';
    end
else
    for k_ss=1:length(ppH)
        ppH{k_ss}.UserData.labelL.Visible='off';
    end
end

end  % end fcn CheckFcnshowssDNAL












% end  %end of all


