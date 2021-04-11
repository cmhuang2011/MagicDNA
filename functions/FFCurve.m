function FFCurve(src,evn,allGraphics,ss_STEP,allbtn)
% this function is used to rapidly generate a mechanism with many
% bundles, select the STEP created from Solidworks and follow the
% instrcutions.

% addpath(genpath(pwd));
clc;

% ax= gca;
h1=findobj(ss_STEP.Children,'tag','SS') ;
h2=findobj(ss_STEP.Children) ;
Removed=setdiff(h2,h1) ;
for k=1:length(Removed)
    delete(Removed(k));
end
cla(allGraphics.allaxes.STEP);

isequal(allbtn.StartSTEP ,src); % use STEP as start

isequal(allbtn.StartPoints ,src); % use XYZ points as start

if isequal(allbtn.StartSTEP ,src)
    %     [FileName,PathName] = uigetfile('*.STEP','Select the STEP file');
    %     fffid=fopen(strcat(PathName,FileName));
    %     Oneline = fgetl(fffid);    count=1;
    %     nn_header=0;
    %     while 1
    %         Oneline = fgetl(fffid); count=count+1;
    %         nn_header=nn_header+1;
    %         if isempty(Oneline) ||  contains(Oneline,{'ENDSEC;'})
    %             break
    %         end
    %     end
    %     nn_rest=nn_header+1;
    %
    %     while 1
    %         Oneline = fgetl(fffid); count=count+1;
    %         if ~isempty(Oneline)
    %             nn_rest=nn_rest+1;
    %         end
    %         if strcmp(Oneline,'DATA;')
    %             DataStart =count ;
    %         end
    %
    %         if strcmp(Oneline,'ENDSEC;')
    %             break
    %         end
    %     end
    %     fclose(fffid);
    %     %  sdfsdf=3
    %     %--------------------------
    %     fffid2=fopen(strcat(PathName,FileName));
    %     for k=1:DataStart %header_____not interesting
    %         Oneline = fgetl(fffid2)   ;
    %     end
    %     field1='number';
    %     field2='type' ;
    %     field3='content' ;
    %     field3='DeleteLine' ;
    %     s(length(nn_header+3 : nn_rest-1)) = struct(field1,[],field2,[],field3,[]);
    %     CriticalWords={'LINE','CARTESIAN_POINT','VECTOR','DIRECTION','TRIMMED_CURVE'};
    %     Exclude={'SPLINE_CURVE'};
    %     localindex=1;
    %     UU=[];
    %     LinesInd=[]; DeleteLine=[]; LineToAdd= [];
    %     table=zeros(length(nn_header+3 : nn_rest-1),2); ToSkip =0;
    %     for iCare=nn_header+3 : nn_rest-1
    %         if ToSkip==0
    %             Oneline = fgetl(fffid2) ;
    %         else
    %             ToSkip=0;
    %             localindex= localindex+1;
    %             continue;
    %         end
    %         poundsing=strfind(Oneline,'#')   ;
    %         if isempty(poundsing)
    %             localindex= localindex+1;
    %             continue;
    %         end
    %         %          if poundsing(1)~=1
    %         %              sdfsf=2
    %         %          end
    %
    %         eqaulsign=strfind(Oneline,'=')  ;
    %         s(localindex).number=  str2num(Oneline(poundsing(1)+1:eqaulsign-1)) ;
    %         if ~isempty( str2num(Oneline(2:eqaulsign-1)))
    %             table(localindex,:)=[localindex,s(localindex).number];
    %         else
    %             table(localindex,:)=[localindex,0];
    %         end
    %         TF = contains(Oneline,CriticalWords) ;   %work for 2016b and later
    %         if TF
    %             UU=union(UU,localindex);
    %             if contains(Oneline,CriticalWords(1))  &&   ~contains(Oneline,Exclude(1))
    %                 %LINE type
    %                 LinesInd=union(LinesInd,localindex);
    %                 s(localindex).type='LINE';
    %                 hash= strfind(Oneline,'#');
    %                 for getnumberi=2:3
    %                     cc=hash(getnumberi)+1;
    %                     while ~isempty( str2num(Oneline(cc)) )
    %                         cc=cc+1   ;
    %                     end
    %                     addnumber=str2num(Oneline( hash(getnumberi)+1:cc-1) );
    %                     s(localindex).content(getnumberi-1)= addnumber ;
    %                 end
    %             elseif contains(Oneline,CriticalWords(2))
    %                 %CARTESIAN_POINT type
    %                 s(localindex).type='CARTESIAN_POINT';
    %                 if strcmp(Oneline(end),';')
    %                     Secondline=[];
    %                 else
    %                     Secondline = fgetl(fffid2) ; ToSkip=1;
    %                 end
    %                 WholeLine =[Oneline Secondline]  ;
    %
    %                 leftP=  strfind(WholeLine,'(') ; leftP=leftP(2);
    %                 rightP=  strfind(WholeLine,')') ; rightP=rightP(1);
    %                 comaP= strfind(WholeLine,',') ;
    %                 comaP(comaP<leftP)=[];
    %                 comaP(comaP>rightP)=[];
    %                 x=str2num(WholeLine(leftP+1:comaP(1)-1));
    %                 y=str2num(WholeLine(comaP(1)+1:comaP(2)-1  ));
    %                 z=str2num(WholeLine(comaP(2)+1:rightP-1  ));
    %                 s(localindex).content=[x,y,z];
    %             elseif contains(Oneline,CriticalWords(3))
    %                 %VECTOR type
    %                 s(localindex).type='VECTOR';
    %                 hash= strfind(Oneline,'#');
    %                 cc=hash(2)+1;
    %                 while ~isempty( str2num(Oneline(cc)) )
    %                     cc=cc+1;
    %                 end
    %                 addnumber=str2num(Oneline( hash(2)+1:cc-1) );
    %                 comaP= strfind(Oneline,',');
    %                 rightP=  strfind(Oneline,')') ;
    %                 Lenghs= str2num(Oneline(comaP(2)+1:rightP-1));
    %                 s(localindex).content=[ addnumber,Lenghs];
    %             elseif contains(Oneline,CriticalWords(4))
    %                 %DIRECTION type
    %                 s(localindex).type='DIRECTION';
    %                 if strcmp(Oneline(end),';')
    %                     Secondline=[];
    %                 else
    %                     Secondline = fgetl(fffid2) ; ToSkip=1;
    %                 end
    %                 WholeLine =[Oneline Secondline]  ;
    %
    %                 %                  Oneline
    %                 leftP=  strfind(WholeLine,'(') ; leftP=leftP(2);
    %                 rightP=  strfind(WholeLine,')') ; rightP=rightP(1);
    %
    %                 comaP= strfind(WholeLine,',') ;
    %                 comaP(comaP<leftP)=[];
    %                 comaP(comaP>rightP)=[];
    %                 x=str2num(WholeLine(leftP+1:comaP(1)-1));
    %                 y=str2num(WholeLine(comaP(1)+1:comaP(2)-1  ));
    %                 z=str2num(WholeLine(comaP(2)+1:rightP-1  ));
    %                 s(localindex).content=[x,y,z];
    %             elseif contains(Oneline,CriticalWords(5))   % Trimmed curve for freeCAD step files
    %                 s(localindex).type='TRIMMED_CURVE';
    %                 if strcmp(Oneline(end),';')
    %                     Secondline=[];
    %                 else
    %                     Secondline = fgetl(fffid2) ; ToSkip=1;
    %                 end
    %                 WholeLine =[Oneline Secondline]  ;
    %                 poundsing=strfind(WholeLine,'#')  ;
    %                 comaP= strfind(WholeLine,',') ;
    %                 s(localindex).DeleteLine= str2num(WholeLine(poundsing(2)+1:comaP(2)-1   )) ;
    %                 %                   mixedSymbol={',' ,')', ' '} ;
    %                 IndsMix = union(union(strfind(WholeLine,',') ,strfind(WholeLine,')') ), strfind(WholeLine,' '));
    %
    %                 Next1=sort( IndsMix(IndsMix>poundsing(3)) );
    %                 Next2=sort( IndsMix(IndsMix>poundsing(4)) );
    %                 Point1 =  str2num(WholeLine(poundsing(3)+1:Next1(1)-1   ))  ;
    %                 Point2 =  str2num(WholeLine(poundsing(4)+1:Next2(1)-1   )) ;
    %                 leftP=  strfind(WholeLine,'(') ;
    %                 rightP=  strfind(WholeLine,')') ;
    %                 LLength= str2num(WholeLine(leftP(5)+1:rightP(3)-1   )) ;
    %                 s(localindex).content=[Point1,Point2,LLength];
    %
    %                 DeleteLine=union(DeleteLine,str2num(WholeLine(poundsing(2)+1:comaP(2)-1   )) ) ;  % ind of object
    %                 LineToAdd=union(LineToAdd, localindex) ;
    %             end
    %         end
    %         localindex= localindex+1;
    %     end
    %     fclose(fffid2);
    %
    %     f22=figure(22);  clf;   %later will be close
    %     hold on;
    %     T2=table(:,2);
    %     T1=table(:,1);
    %     LineToDelete =table(ismember(table(:,2),DeleteLine  ) ,1) ;
    %
    %     LinesInd=setdiff(LinesInd, LineToDelete) ;
    %     coordinates=zeros(2*length(LinesInd)+2*length(LineToAdd)  ,3);
    %     % % % %      coordinates=zeros(2*length(LinesInd),3);
    %     for Li=1:length(LinesInd)
    %         LL= LinesInd(Li);
    %
    %         CartInd=s(LL).content(1);
    %         CorreInd=T1( T2==CartInd);
    %         PA= s(CorreInd).content;
    %
    %         VecInd=s(LL).content(2);
    %         CV=T1( T2==VecInd);
    %         Lenghts=s(CV).content(2) ;
    %
    %         DirInd=s(CV).content(1);
    %         CorrDir=T1( T2==DirInd);
    %         Dir=s(CorrDir).content;
    %
    %         PB=PA+Lenghts*Dir;       PP=[PA;PB];
    %         PP= circshift(PP,1,2);
    %         plot3(PP(:,1),PP(:,2),PP(:,3));
    %         coordinates(2*Li-1:2*Li,:)=PP;
    %     end
    %
    %     for Lj=length(LinesInd)+1: length(LineToAdd)
    %         LL2 = LineToAdd(Lj-length(LinesInd)) ;
    %
    %         PA_objInd =  s(LL2).content(1) ;
    %         PA_CorreInd=T1( T2==PA_objInd);
    %         PA= s(PA_CorreInd).content;
    %
    %         PB_objInd =  s(LL2).content(2) ;
    %         PB_CorreInd=T1( T2==PB_objInd);
    %         PB= s(PB_CorreInd).content;
    %
    %         PP=[PA;PB];
    %         PP= circshift(PP,1,2);
    %         plot3(PP(:,1),PP(:,2),PP(:,3));
    %         coordinates(2*Lj-1:2*Lj,:)=PP;
    %     end
    %
    %     xlabel('x');ylabel('y');zlabel('z') ;%axis equal;
    %     axis auto;
    %
    %
    %     %       clearvars -except  coordinates edges  f22 allGraphics ss_STEP
    %     %       Li=
    %     % % 04082019
    %
    %
    %     %-------------------------------
    %     [C,~,IC]=uniquetol(coordinates,'ByRows',true);
    %     shp_test = alphaShape(C(:,1),C(:,2),C(:,3))  ;   [tri_test,P_test] =  boundaryFacets(shp_test);
    %     if isempty(P_test) % planar points
    %         %          extraP=10*rand(1,3) ;  % don't use stochastic for consistence
    %         [coeff,score,~] = pca(coordinates) ;
    %         extraP=mean(score) + std(std(score))*coeff(:,3)';
    %         C=[C;extraP];
    %     end
    %     if isempty(Li)
    %         Li=0 ;
    %         edges=[   2*Li+1:2:2*Li+2*Lj-1  ; 2*Li+2:2:2*Li+2*Lj ]' ;   %
    %
    %     else
    %         edges=[ 1:2:2*Li-1  ; 2:2:2*Li ;  2*Li+1:2:2*Li+2*Lj-1  ; 2*Li+2:2:2*Li+2*Lj ]' ;   %
    %
    %     end
    %
    %     edges2=zeros(size(edges));
    %     for k=1:size(edges2,1)
    %         edges2(k,1)= IC(edges(k,1));
    %         edges2(k,2)= IC(edges(k,2));
    %     end
    %     shp = alphaShape(C(:,1),C(:,2),C(:,3))  ;
    %     %       figure ; plot(shp);  %--------------
    %     [tri,P] =  boundaryFacets(shp);
    %
    %     if ~isempty(P)
    %         if ~isempty(setdiff(C,P,'rows'))
    %             P=[P;setdiff(C,P,'rows')]    ;
    %         end
    %     end
    %
    %     if rank(C(1:end-1,:))==2 && isempty(P_test)
    %         [~,  IndExtra]=ismember(extraP,P,'rows')  ;
    %         VecRelata=setdiff(1:size(P,1),IndExtra) ;
    %         [u,v]=find(tri==IndExtra)  ;
    %         tri(u,:)=[];
    %         for i2=1:size(tri,1)
    %             for j2=1:size(tri,2)
    %                 tri(  i2,j2)=  find(tri(  i2,j2)==VecRelata);
    %             end
    %         end
    %         P(IndExtra,:)=[];
    %     end
    %
    %     for jj=1:size(P,1)
    %         text(P(jj,1),P(jj,2),P(jj,3),  num2str(jj))
    %     end
    %
    %     edges3=zeros(size(edges2));
    %     for update3=1:size(edges2,1)
    %         PA=C(edges2(update3,1),:) ;
    %         PB=C(edges2(update3,2),:);
    %         [~,indA]=ismember(PA,P,'rows');
    %         [~,indB]=ismember(PB,P,'rows');
    %         if indA==0 || indB==0
    %             %             sdsfsf=234;
    %         end
    %         edges3(update3,:)=[indA,indB];
    %     end
    %
    %     faces2=cell(size(tri,1),2);
    %     for k=1:size(faces2,1)
    %         faces2{k,1}=3;
    %         faces2{k,2}=tri(k,:);
    %     end
    %
    %     coordinates=P;
    %     faces=faces2;
    %     edges=edges3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %     elseif isequal(allbtn.StartPoints ,src)  % use XYZ points as start
    
    readPoints_2(allbtn.StartPoints,[],allGraphics,ss_STEP) ;
    coordinates=  ss_STEP.UserData.UsePoints.pXYZ ;
    edges = ss_STEP.UserData.UsePoints. Edges ;
    f22=figure(22);  clf;   %later will be close
    scatter3(coordinates(:,1),coordinates(:,2),coordinates(:,3)) ; axis equal ;
    %           f23=figure(23);  clf;   %later will be close
    %        scatter3(coordinates(:,1),coordinates(:,2),coordinates(:,3)) ; axis equal ;
    [C,~,IC]=uniquetol(coordinates,'ByRows',true);
    if rank(C)==2
        %          extraP=2*norm(diff(coordinates))*rand(1,3) ;
        [coeff,score,latent] = pca(coordinates) ;
        extraP=mean(coordinates) + cross(coeff(:,1),coeff(:,2))';
        
        C=[C;extraP];
    else
        extraP=[];
    end
    shp = alphaShape(C(:,1),C(:,2),C(:,3))  ;
    %  figure ; plot(shp);
    [tri,P] =  boundaryFacets(shp);
    if ~isempty(P)
        if ~isempty(setdiff(C,P,'rows'))
            P=[P;setdiff(C,P,'rows')]    ;
        end
    end
    if rank(C(1:end-1,:))==2
        [~,  IndExtra]=ismember(extraP,P,'rows')  ;
        VecRelata=setdiff(1:size(P,1),IndExtra) ;
        [u,v]=find(tri==IndExtra)  ;
        tri(u,:)=[];
        for i2=1:size(tri,1)
            for j2=1:size(tri,2)
                tri(  i2,j2)=  find(tri(  i2,j2)==VecRelata);
            end
        end
        P(IndExtra,:)=[];
    end
    
    faces2=cell(size(tri,1),2);
    for k=1:size(faces2,1)
        faces2{k,1}=3;
        faces2{k,2}=tri(k,:);
    end
    faces=faces2;
    
end    % end of using STEP or sketch points
clearvars -except  coordinates edges faces f22 f23 allGraphics ss_STEP
% 04082019


% edge_length_PLY

edge_length_PLY = zeros(size(edges,1),1);
for edge_ID = 1:size(edges,1)
    edge_bgn = edges(edge_ID,1) ;
    edge_fin = edges(edge_ID,2);
    
    edge_length_PLY(edge_ID) = norm(coordinates(edge_bgn,:) - coordinates(edge_fin,:));
end
P= coordinates;
[min_edge_PLY, ~] = min(edge_length_PLY);

% %----------sorting lengths,   09142018
% [NewEdges0914,order1] = sort(edge_length_PLY) ;
% NewEdge = edges(order1,:) ;
% edges=  NewEdge ;
% %  [NewEdge(:) , edges(:)]
%-------------------sorting lengths,09142018

Explode_ratio=0.8;     % enlarge the assembly, need to be asked ealier
min_len_nt=[];

if isempty(min_len_nt)
    min_len_nt = 21;
end
scale = min_len_nt/min_edge_PLY*Explode_ratio;
%    scale=1;
Sca_coor=coordinates*scale ;
close(f22);


axes(allGraphics.allaxes.STEP); hold on;
%    fH=figure(35)  ;clf;hold on; axis equal;axis auto;
%    fH.Position=[680 49 900 900];
ax=gca;  ax.Position(3)=0.55; ax.Position(1)=0.05;
sH=scatter3(Sca_coor(:,1) ,Sca_coor(:,2),Sca_coor(:,3)) ;
sH.HitTest='off';

tH=cell(size(coordinates,1),1);
for tt=1:size(coordinates,1)
    %    tH{tt}=text(Sca_coor(tt,1) ,Sca_coor(tt,2),Sca_coor(tt,3),num2str(tt)   ) ;
    %    tH{tt}.HitTest='off';
end
%
edgeH=cell(size(edges,1),1);
for edgei=1:size(edges,1)
    A=Sca_coor(edges(edgei,1),:) ;
    B=Sca_coor(edges(edgei,2),:) ;
    C=[A;B];
    edgeH{edgei}= plot3(C(:,1),C(:,2),C(:,3));
    edgeH{edgei}.Color=[0.2,0.2,0.2];
end

edgeBelongFace=cell(size(edges,1),1);
for edgecheckface=1:size(edges,1)
    thisedge=edges(edgecheckface,:);
    for kface=1:size(faces,1)
        [A,~]=ismember(thisedge,faces{kface,2}) ;
        if sum(A)==2
            edgeBelongFace{edgecheckface}=union(edgeBelongFace{edgecheckface},kface);
        end
    end
end

box=[max(Sca_coor);min(Sca_coor) ];
DDbox=box(1,:)-box(2,:) ;

FaceNormal=zeros(size(faces,1),3);
for facei=1:size(faces,1)
    P1=Sca_coor(faces{facei,2}(1),:);
    P2=Sca_coor(faces{facei,2}(2),:);
    P3=Sca_coor(faces{facei,2}(3),:);
    NVec=cross( (P2-P1) ,(P3-P2)) ;
    NVec=NVec/norm(NVec) ;
    FaceNormal(facei,1:3)=NVec;
end

for scanFaceValid= 1:length(edgeBelongFace)  % 05152019
    QQ=edgeBelongFace{scanFaceValid} ;
    %     FaceNormal(QQ,1);
    IndNan=  isnan(FaceNormal(QQ,1) ) ;
    edgeBelongFace{scanFaceValid}=   QQ(~IndNan) ;
end

edge_ZYX_Vec=zeros(size(edges,1) ,9) ;
edge_center=zeros(size(edges,1) ,3) ;
for edgej=1:size(edges,1)
    thisedge=edges(edgej,:);
    ZVec=Sca_coor(edges(edgej,2),:)-Sca_coor(edges(edgej,1),:) ;
    ZVec=ZVec/norm(ZVec);
    
    %--  find corresponding two faces
    twofaces=edgeBelongFace{edgej};
    length(twofaces);
    if length(twofaces)==2
        
        N1=FaceNormal(twofaces(1) ,:) ;
        N2=FaceNormal(twofaces(2) ,:) ;
        if sum(N1+N2==[0,0,0])==3
            N2=[0,0,0] ;
        end
        YVec= cross(      ZVec, N1+N2) ;YVec =YVec/norm(YVec) ;
        XVec=cross(YVec,ZVec)  ;
        
        
        % hard code
        YVec = rand(1,3) ; 
        YVec=cross(ZVec , cross(YVec, ZVec) ) ;
        XVec=cross(YVec,ZVec)  ;
        %-------
    elseif  length(twofaces)==1    %------1
        N1=FaceNormal(twofaces(1) ,:)   ;
        YVec= cross(      ZVec, N1) ;YVec =YVec/norm(YVec) ;
        XVec=cross(YVec,ZVec)  ;
    else
        S=setdiff(edges , thisedge,'rows') ;
        nw=1 ;
        while 1   % incase select wrong reference point,   09142018
            if  ismember(thisedge(1),S) && rand<0.01
                EEE=thisedge;
                %                 rareHappen =1
            else
                EEE=flip(thisedge);
            end
            [aa,~]=find(S==EEE(1)) ;
            if isempty(aa)
                %                 RadVec=rand(1,3);RadVec=RadVec/norm(RadVec) ;
                if DDbox(1)==min(DDbox)
                    RedVec=[1,0,0];
                elseif DDbox(2)==min(DDbox)
                    RedVec=[0,1,0];
                else
                    RedVec=[0,0,1];
                end
                YVec=cross(    ZVec, RedVec) ;YVec =-YVec/norm(YVec) ;
                XVec=cross(YVec,ZVec)  ;
            else
                aa=aa(1);
                ZA=P(EEE(1),:)-P(EEE(2),:) ;
                RefPoint=setdiff(S(aa,:),EEE(1));
                ZB=P(EEE(1),:)-P(RefPoint,:) ;
                YVec= cross(    ZA, ZB) ;YVec =-YVec/norm(YVec) ;
                YVec=-cross(YVec,ZVec)  ;YVec =YVec/norm(YVec) ;
                XVec=cross(YVec,ZVec)  ;
            end
            if sum(~isnan(YVec))==3
                break;
            end
            
            if  nw==200
                RandV = rand(1,3) ;
                YVec=cross(ZVec,RandV)  ;YVec =YVec/norm(YVec) ;
                XVec=cross(YVec,ZVec)  ;
                break;
            end
            
            nw=nw+1 ;
        end  % incase select wrong reference point,   09142018
    end
    %     if isnan(YVec)
    %         sdfsf=3;
    %     end
    
    
    edge_ZYX_Vec(  edgej,:)=[ZVec, YVec,XVec];
    edge_center(edgej,:)=  mean([Sca_coor(thisedge(1),:) ; Sca_coor(thisedge(2),:)]) ;
    
    %          if ismember(edgej,[19,20])
    %              sdf=4
    %          end
    
end
%--visual edge info
local_triadL = 5 ;
LocalXYZ_h = cell(size(edges,1),3 ) ;
for ke=1:size(edges,1)
    G=edge_center(ke,:);
    V=local_triadL*edge_ZYX_Vec(ke,:);
    Zline=[G ; G+V(1:3)] ; Yline=[G ; G+V(4:6)] ;Xline=[G ; G+V(7:9)] ;
    LocalXYZ_h{ke,3}=plot3(Zline(:,1),Zline(:,2),Zline(:,3),'b'); LocalXYZ_h{ke,3}.LineWidth=4; %LocalXYZ_h{ke,3}.HitTest='off';
    LocalXYZ_h{ke,2}=plot3(Yline(:,1),Yline(:,2),Yline(:,3),'g');  LocalXYZ_h{ke,2}.LineWidth=4; % LocalXYZ_h{ke,2}.HitTest='off';
    LocalXYZ_h{ke,1}=plot3(Xline(:,1),Xline(:,2),Xline(:,3),'r');  LocalXYZ_h{ke,1}.LineWidth=4; % LocalXYZ_h{ke,1}.HitTest='off';
    
    LocalXYZ_h{ke,3}.UserData.Edgei_XYZ = [ke,3] ;    % label for distinguish .
    LocalXYZ_h{ke,2}.UserData.Edgei_XYZ = [ke,2] ;
    LocalXYZ_h{ke,1}.UserData.Edgei_XYZ = [ke,1] ;
    
end
% for ke=1:size(edges,1)
%     LocalXYZ_h{ke,1}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,LocalXYZ_h) ;
%     LocalXYZ_h{ke,2}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,LocalXYZ_h) ;
%     LocalXYZ_h{ke,3}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,LocalXYZ_h) ;
% end
% edge_ZYX_Vec(10:13,:)


hp1 = uipanel('Parent',ss_STEP,'Title','Insert Truss','FontSize',12,...
    'Position',[0.7  0.75   0.25 0.2]);
hp1.UserData.willbeD=1;

btnSelectEdge1= uicontrol('Style','togglebutton','Parent',hp1,'Units','normalized',...
    'Position',[0.1 0.7 0.5 0.2],  'String','Select Edge 1');
btnSelectEdge1.UserData.willbeD=1;
editH1 = uicontrol('Style', 'edit','Parent',hp1,'String', '5','Unit','normalized','Position', [0.7 0.7 0.3 0.2]);
editH1.UserData.willbeD=1;

btnSelectEdge2= uicontrol('Style','togglebutton','Parent',hp1,'Units','normalized',...
    'Position',[0.1 0.2 0.5 0.2],  'String','Select Edge 2');
btnSelectEdge2.UserData.willbeD=1;
editH2 = uicontrol('Style', 'edit','Parent',hp1,'String', '6','Unit','normalized','Position', [0.7 0.2 0.3 0.2]);
editH2.UserData.willbeD=1;


btnAdd= uicontrol('Style','pushbutton','Parent',ss_STEP,'Units','normalized',...
    'Position',[0.8 0.6 0.1 0.1],  'String','Add Truss');
btnAdd.UserData.willbeD=1;
btnAdd.Enable='off';
btnAdd.Callback=@(src,evn)AddTruss(src,evn,edgeH,edge_ZYX_Vec) ;
ax.UserData.btnAddH=btnAdd;
ax.UserData.LocalXYZ_h=LocalXYZ_h;

btnTogSelectorDelete= uicontrol('Style','togglebutton','Parent',ss_STEP,'Units','normalized',...
    'Position',[0.8 0.48 0.1 0.1],  'String','Add/Delete Edge');
btnTogSelectorDelete.Callback=@(src,evn)Select3(src,evn,edgeH,btnSelectEdge1,btnSelectEdge2 ) ;
btnTogSelectorDelete.UserData.willbeD=1;

btnSelectEdge1.Callback=@(src,evn)Select1(src,evn,edgeH,btnSelectEdge2, editH1,ax,btnTogSelectorDelete ) ;
btnSelectEdge2.Callback=@(src,evn)Select2(src,evn,edgeH,btnSelectEdge1, editH2,ax,btnTogSelectorDelete) ;

ax.UserData.extraEdge=[];

str1=strcat('Total Bundle = ',num2str(size(edges,1)) );
txt1 = uicontrol('Parent',ss_STEP,'Style','text','Units','normalized',...
    'Position',[0.7 0.3 0.3 0.15], 'String',str1);
ax.UserData.textH=txt1; txt1.UserData.willbeD=1;

btnExport= uicontrol('Style','pushbutton','Parent',ss_STEP,'Units','normalized',...
    'Position',[0.8 0.2 0.1 0.1],  'String','Finish');
btnExport.Callback=@(src,evn)exprotBundle(src,evn,edge_center,edge_ZYX_Vec,edgeH)  ;
btnExport.UserData.willbeD=1;
checkH = uicontrol('Parent',ss_STEP,'Style', 'checkbox','String', 'Axis auto/equal','Unit','normalized','Position', [0.7  0.7   0.2 0.05]);
checkH.Callback=@(src,evn)checkFcn(src,evn)  ;


%---------add legend as instructions, single row
ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
hLg= legend(ForLegend,'Click me for instructions','Location','northwest') ;

hLg.String={'\bf{Click me for instructions}'};
hLg.Interpreter='tex';  %latex
hLg.Orientation='horizontal';
ForLegend.Marker='.' ; ForLegend.Marker='none';
hLg.ButtonDownFcn=@(src,evn)LegendBoxing_SketchStep1( src,evn,ax );

hLg.Units='normalized'; hLg.AutoUpdate ='off';
hLg.Position=[0.0063 0.9728 0.1569 0.025];
%------------------

fH=gcf ; fH.UserData.MovieStep =2 ;
% fH.KeyPressFcn = @(src,evn)KeyPressEditSketch(src,evn) ;
axis equal;
ss_STEP.UserData.New_tarr=[];
exprotBundle(btnExport,[],edge_center,edge_ZYX_Vec,edgeH)   ;   % directly go to converting lines to bundles
end

% function KeyPressEditSketch(src,evn)
% switch evn.Key
%     case 'h'
%         implay('EditSketchSTEP.mp4');
% end
% end

function  ChangeLocalOrientation(src,evn,LocalXYZ_h)
% clc;
Edgei_XYZ= src.UserData.Edgei_XYZ  ; % [edgeInd, {1, 2,3 }=[xyz] ]
ax=gca;
DirectionMag = 5 ;
% Edgei_XYZ(1) = find(ax.UserData.RemainOriEdge==Edgei_XYZ(1) ) ;
% evn
if evn.Button==1   % left click -> flipping x and z direction
    XvecData = [LocalXYZ_h{Edgei_XYZ(1) ,1}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,1}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,1}.ZData;]  ;
    ZvecData = [LocalXYZ_h{Edgei_XYZ(1) ,3}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,3}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,3}.ZData;]  ;
    dirx = diff(XvecData')'; dirz =  diff(ZvecData')';  dirx=dirx/norm(dirx) ; dirz=dirz/norm(dirz) ;
    XvecData(:,2) =  XvecData(:,1)-  dirx*DirectionMag ;  
    ZvecData(:,2) =  ZvecData(:,1)-  dirz*DirectionMag ; 
    
%     XvecData(:,2) = 2* XvecData(:,1)  -  XvecData(:,2) ;
%     ZvecData(:,2) = 2* ZvecData(:,1)  -  ZvecData(:,2) ;
    LocalXYZ_h{Edgei_XYZ(1) ,1}.XData=XvecData(1,:) ;
    LocalXYZ_h{Edgei_XYZ(1) ,1}.YData=XvecData(2,:) ;
    LocalXYZ_h{Edgei_XYZ(1) ,1}.ZData=XvecData(3,:) ;
    
    LocalXYZ_h{Edgei_XYZ(1) ,3}.XData=ZvecData(1,:) ;
    LocalXYZ_h{Edgei_XYZ(1) ,3}.YData=ZvecData(2,:) ;
    LocalXYZ_h{Edgei_XYZ(1) ,3}.ZData=ZvecData(3,:) ;
    
elseif evn.Button==3   % right click -> rotating about z-axis 90 degree
    XvecData = [LocalXYZ_h{Edgei_XYZ(1) ,1}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,1}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,1}.ZData;]  ;
    YvecData = [LocalXYZ_h{Edgei_XYZ(1) ,2}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,2}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,2}.ZData;]  ;
%     XvecData(:,2) = 2* XvecData(:,1)  -  XvecData(:,2) ;
    dirx = diff(XvecData')'; diry =  diff(YvecData')';  dirx=dirx/norm(dirx) ; diry=diry/norm(diry) ;
    XvecData(:,2) =  XvecData(:,1)-  dirx*DirectionMag ; 
    YvecData(:,2) =  YvecData(:,1)+  diry*DirectionMag ; 
    
    LocalXYZ_h{Edgei_XYZ(1) ,1}.XData=YvecData(1,:) ;
    LocalXYZ_h{Edgei_XYZ(1) ,1}.YData=YvecData(2,:) ;
    LocalXYZ_h{Edgei_XYZ(1) ,1}.ZData=YvecData(3,:) ;
    
    LocalXYZ_h{Edgei_XYZ(1) ,2}.XData=XvecData(1,:) ;
    LocalXYZ_h{Edgei_XYZ(1) ,2}.YData=XvecData(2,:) ;
    LocalXYZ_h{Edgei_XYZ(1) ,2}.ZData=XvecData(3,:) ;
    
elseif    evn.Button==2   % middle click -> custom angle  % 03052020
    
    
    
    XvecData = [LocalXYZ_h{Edgei_XYZ(1) ,1}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,1}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,1}.ZData;]  ;
    YvecData = [LocalXYZ_h{Edgei_XYZ(1) ,2}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,2}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,2}.ZData;]  ;
    
    Vx = [LocalXYZ_h{Edgei_XYZ(1) ,1}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,1}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,1}.ZData;]   ;
    Vx = diff(Vx,1,2)   ; Vx=Vx/norm(Vx) ;
    Vy = [LocalXYZ_h{Edgei_XYZ(1) ,2}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,2}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,2}.ZData;]  ;
    Vy = diff(Vy,1,2)   ; Vy=Vy/norm(Vy) ;
    Vz = [LocalXYZ_h{Edgei_XYZ(1) ,3}.XData ; LocalXYZ_h{Edgei_XYZ(1) ,3}.YData ;LocalXYZ_h{Edgei_XYZ(1) ,3}.ZData;]  ;
    Vz = diff(Vz,1,2)   ; Vz=Vz/norm(Vz) ;
    
    Result = AskLocalSnapOrRotationInterval ;
    if strcmp(Result.Decision ,'Rotate')
        LocalAxis=Vz ;   theta= Result.RotationInterval ;
        if isempty(theta)
            fprintf('Invalid input!! \n')
            warndlg('Invalid input!!','Warning');
            return
        end
        
        RMat = RotationAxisToRotaionMatrix( LocalAxis,theta )      ;
        NewVxVy = [RMat*Vx  RMat*Vy  ]  ;
        
        X_3x2= [ XvecData(:,1) ,XvecData(:,1) + 2* NewVxVy(:,1)   ] ;
        Y_3x2= [ YvecData(:,1) ,YvecData(:,1) + 2* NewVxVy(:,2)  ] ;
        LocalXYZ_h{Edgei_XYZ(1) ,1}.XData=X_3x2(1,:) ;    LocalXYZ_h{Edgei_XYZ(1) ,1}.YData=X_3x2(2,:) ;    LocalXYZ_h{Edgei_XYZ(1) ,1}.ZData=X_3x2(3,:) ;
        LocalXYZ_h{Edgei_XYZ(1) ,2}.XData=Y_3x2(1,:) ;    LocalXYZ_h{Edgei_XYZ(1) ,2}.YData=Y_3x2(2,:) ;    LocalXYZ_h{Edgei_XYZ(1) ,2}.ZData=Y_3x2(3,:) ;
    elseif strcmp(Result.Decision ,'Snap')
        
        LocalAxis=Vz ;   thetas=-45:0.1:45 ;
        RMats = zeros(3,3,length(thetas)) ;
        NewVxs= zeros(3,length(thetas)) ;
        for k=1:length(thetas)
            RMats(:,:,k) = RotationAxisToRotaionMatrix( LocalAxis,thetas(k) )      ;
            NewVxs(:,k) =   RMats(:,:,k)*Vx      ;
        end
        %     NewVxVy = [RMat*Vx  RMat*Vy  ]  ;
        
        %         RMats
        NewVxs ;
        SnapXto=[ 1 ,0,0 ; -1 , 0, 0;   0 ,1,0 ;0 , -1, 0; 0 ,0,1 ;0 , 0, -1 ]';
        MaxVal=zeros(6,1) ;
        M1= repmat(SnapXto(:,1), 1,size(NewVxs,2));MaxVal(1) = max(sum(NewVxs.*M1 ) );
        M2= repmat(SnapXto(:,2), 1,size(NewVxs,2));MaxVal(2) = max(sum(NewVxs.*M2 ) );
        M3= repmat(SnapXto(:,3), 1,size(NewVxs,2));MaxVal(3) = max(sum(NewVxs.*M3 ) );
        M4= repmat(SnapXto(:,4), 1,size(NewVxs,2));MaxVal(4) = max(sum(NewVxs.*M4 ) );
        M5= repmat(SnapXto(:,5), 1,size(NewVxs,2));MaxVal(5) = max(sum(NewVxs.*M5 ) );
        M6= repmat(SnapXto(:,6), 1,size(NewVxs,2));MaxVal(6) = max(sum(NewVxs.*M6 ) );
        MaxVal;
        
        IndMaxCase =  find(MaxVal==max(MaxVal)); IndMaxCase=IndMaxCase(1) ;
        %         IndMaxCase=SnapXto(IndMaxCase,:) ;
        Mx_Max = repmat(SnapXto(:,IndMaxCase), 1,size(NewVxs,2));
        Vals= sum(NewVxs.*Mx_Max );
        IndSnap=find(Vals==max(Vals));IndSnap=IndSnap(1) ;
        RMatBest = RMats(:,:,IndSnap) ;
        NewVxVy = [RMatBest*Vx  RMatBest*Vy  ]  ;
        
        X_3x2= [ XvecData(:,1) ,XvecData(:,1) + 2* NewVxVy(:,1)   ] ;
        Y_3x2= [ YvecData(:,1) ,YvecData(:,1) + 2* NewVxVy(:,2)  ] ;
        LocalXYZ_h{Edgei_XYZ(1) ,1}.XData=X_3x2(1,:) ;    LocalXYZ_h{Edgei_XYZ(1) ,1}.YData=X_3x2(2,:) ;    LocalXYZ_h{Edgei_XYZ(1) ,1}.ZData=X_3x2(3,:) ;
        LocalXYZ_h{Edgei_XYZ(1) ,2}.XData=Y_3x2(1,:) ;    LocalXYZ_h{Edgei_XYZ(1) ,2}.YData=Y_3x2(2,:) ;    LocalXYZ_h{Edgei_XYZ(1) ,2}.ZData=Y_3x2(3,:) ;
        
        
    end
    
end




end


function checkFcn(src,~)
if src.Value==1
    axis normal;
    %     xlim auto; ylim auto; zlim auto;
else
    axis equal;
end
end



function exprotBundle(~,~,edge_center,edge_ZYX_Vec,edgeH)
ax=gca;

% fH=gcf;
fH= findobj(gcf,'tag','ss_STEP');
title(ax,'')
% ------scanning and ignore deleted edges
% Inds=[];
% for k=1:length(edgeH)
%     if strcmp(edgeH{k}.LineStyle,'-' )
%         Inds=union(Inds,k);
%     end
% end
% Inds
% edgeH=edgeH(Inds);

ax.UserData.edgeHhandle=edgeH;

%------
for ee1=1:length(edgeH)
    if   isfield(edgeH{ee1}.UserData,'scatterH')
        for ss1=1:length(edgeH{ee1}.UserData.scatterH)
            if isvalid(edgeH{ee1}.UserData.scatterH{ss1})
                edgeH{ee1}.UserData.scatterH{ss1}.Visible='off';
            end
        end
    end
end

ccf=findobj(fH);
for cci=1:length(ccf)
    if   isfield(ccf(cci).UserData,'willbeD')
        ccf(cci).Visible='off'  ;
    end
end

N=0;
RemainOriEdge = [] ;
%----Original Edge
for k=1:length(edgeH)
    if (sum(edgeH{k}.Color==[1,0,0])==3  ||  sum(edgeH{k}.Color==[0,1,0])==3 || sum(edgeH{k}.Color==[0.2,0.2,0.2])==3 ) && strcmp(edgeH{k}.LineStyle,'-')
        edgeH{k}.Color ;
        N=N+1;
        RemainOriEdge=union(RemainOriEdge,k) ;
    end
end
%------Extra Edge
k2=0;
if isfield(ax.UserData,'extraEdge')
    for k2=1:length( ax.UserData.extraEdge)
        cc=ax.UserData.extraEdge{k2}.plotH.Color;
        if sum(cc==[0.2,0.2,0.2])==3
            N=N+1;    RemainOriEdge=union(RemainOriEdge,k+k2) ;
        end
        ax.UserData.extraEdge{k2}.plotH.ButtonDownFcn='';
    end
end
%----------------------------
if isempty(k2);k2=0;end

for k=1:length(edgeH)+k2
    if ~ismember(k,RemainOriEdge)
        delete(ax.UserData.LocalXYZ_h{k,1}) ;
        delete(ax.UserData.LocalXYZ_h{k,2}) ;
        delete(ax.UserData.LocalXYZ_h{k,3}) ;
        
    end
end
ax.UserData.RemainOriEdge=RemainOriEdge ;
ax.UserData.LocalXYZ_h= ax.UserData.LocalXYZ_h(RemainOriEdge,:) ;




for ke=1:size(ax.UserData.LocalXYZ_h,1)
    ax.UserData.LocalXYZ_h{ke,3}.UserData.Edgei_XYZ = [ke,3] ;    % label for distinguish .
    ax.UserData.LocalXYZ_h{ke,2}.UserData.Edgei_XYZ = [ke,2] ;
    ax.UserData.LocalXYZ_h{ke,1}.UserData.Edgei_XYZ = [ke,1] ;
end
for ke=1:size(ax.UserData.LocalXYZ_h,1)
    ax.UserData.LocalXYZ_h{ke,1}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,ax.UserData.LocalXYZ_h) ;
    ax.UserData.LocalXYZ_h{ke,2}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,ax.UserData.LocalXYZ_h) ;
    ax.UserData.LocalXYZ_h{ke,3}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,ax.UserData.LocalXYZ_h) ;
end


PointsInfo=zeros(N,9);
VecInfo=zeros(N,9);  count=1;
EdgeRelativeTable=zeros(N,2);

ccOri=0;
for k=1:length(edgeH)
    if sum(edgeH{k}.Color==[1,0,0])==3  ||  sum(edgeH{k}.Color==[0,1,0])==3 || sum(edgeH{k}.Color==[0.2,0.2,0.2])==3
        PointsInfo(count,1:6)=[edgeH{k}.XData(1),edgeH{k}.YData(1),edgeH{k}.ZData(1),edgeH{k}.XData(2),edgeH{k}.YData(2),edgeH{k}.ZData(2)] ;
        PointsInfo(count,7:9)= edge_center(k,:);
        VecInfo(count,1:9)=edge_ZYX_Vec(k,:);
        EdgeRelativeTable(k,:)=[k,count];
        count=count+1;
        ccOri=ccOri+1;
        
        %----
        edgeH{k}.Color=[0,0,0];
        edgeH{k}.LineStyle='-';
        tH=text(edge_center(k,1),edge_center(k,2),edge_center(k,3), strcat('\leftarrow', num2str(count-1)),'HitTest','off' , 'clipping', 'on' );
        %       tH=annotation('textarrow',edge_center(k,1),edge_center(k,2),edge_center(k,3},'String','y = x ')
        
        %       tH.Position=tH.Position + 2*rand(1,3) ;
        tH.FontSize=18 ;
        %       sdfsf=234
    end
end
RemainExtraEdge=[];
if isfield(ax.UserData,'extraEdge')
    for k2=1:length( ax.UserData.extraEdge)
        cc=ax.UserData.extraEdge{k2}.plotH.Color;
        if sum(cc==[0.2,0.2,0.2])==3
            PointsInfo(count,1:6)=ax.UserData.extraEdge{k2}.XYZ ;
            PointsInfo(count,7:9)= ax.UserData.extraEdge{k2}.midP ;
            VecInfo(count,1:9)=ax.UserData.extraEdge{k2}.edge_ZYX_Vec ;
            EdgeRelativeTable(k+k2,:)=[k+k2,count];
            count=count+1;
            RemainExtraEdge=union(RemainExtraEdge,k2);
            %----
            ax.UserData.extraEdge{k2}.plotH.Color=[0,0,0];
            ax.UserData.extraEdge{k2}.plotH.LineStyle='-';
            text(ax.UserData.extraEdge{k2}.midP(1),ax.UserData.extraEdge{k2}.midP(2),ax.UserData.extraEdge{k2}.midP(3), strcat('\leftarrow',num2str(count-1)) ,'HitTest','off','FontSize',18, 'clipping', 'on'  );
        end
    end
end

ax.UserData.ExtarEdgeRemain=RemainExtraEdge;

%-------------
ax.UserData.Final.PointsInfo=PointsInfo;
ax.UserData.Final.VecInfo=VecInfo;
ax.UserData.RelateTable=EdgeRelativeTable;

DistBtwTwoEnds=zeros(N,1)  ;   %relative lengths unit
for dj=1:N
    DistBtwTwoEnds(dj)=norm(PointsInfo(dj,4:6)- PointsInfo(dj,1:3)) ;
end
PitchConst=0.34;

if ccOri~=count
    trusstr=num2str(min(DistBtwTwoEnds(ccOri+1:end)));
else
    trusstr='NA'    ;
end
% ArrEstimate=[min(DistBtwTwoEnds),min(DistBtwTwoEnds(1:ccOri)),
% return

% strsum=strcat('Minimum length of [All, Edge, truss]= [',num2str(min(DistBtwTwoEnds)),',',num2str(min(DistBtwTwoEnds(1:ccOri))), ',',...
%     trusstr, ']');
% strsum2=' Note: they are not final dimensions, please enter scales for blow-view and enlarge';
% strAl={strsum;strsum2};
d =ax.Parent;

% txt1 = uicontrol('Parent',ax.Parent,'Style','text','Units','normalized',...
%     'Position',[0.7 0.8  0.1 0.2],   'String',strAl,'Visible','off' );
strA2=strcat('After enlarge =>  [',num2str(min(DistBtwTwoEnds)), ']');
txt2 = uicontrol('Parent',ax.Parent,'Style','text','Units','normalized',...
    'Position',[0.82 0.8  0.1 0.2],   'String',strA2,'Visible','off');


sld1 = uicontrol('Style', 'slider','Parent',d,'Units','normalized',....
    'Min',0.5,'Max',4,'Value',1.09,'Position', [0.7 0.65 0.25 0.05]);    %  blow
sld2 = uicontrol('Style', 'slider','Parent',d,'Units','normalized',....
    'Min',0.05,'Max',4,'Value',1,'Position', [0.7 0.75 0.25 0.05]);     % scale

sld1.SliderStep=[0.002,0.05];
sld2.SliderStep=[0.002,0.05];

txts1 = uicontrol('Parent',d,'Style','text','Units','normalized',...
    'Position',[0.8 0.6 0.1 0.05],   'String','S_Blow = 1');
txts2 = uicontrol('Parent',d,'Style','text','Units','normalized',...
    'Position',[0.8 0.7 0.1 0.05],   'String','S_Enlarge = 1');

sld1.Callback= @(src,evn)sldscale1(src,evn,txts1) ;
btnExport2= uicontrol('Style','pushbutton','Parent',fH,'Units','normalized',...
    'Position',[0.82 0.02 0.1 0.13],  'String','Export Bundles');

InitTBP= sum(round(DistBtwTwoEnds/0.34))*4;
txtTotalBP = uicontrol('Parent',d,'Style','text','Units','normalized',...
    'Position',[0.92 0.05 0.05 0.1],   'String',strcat('Apx. Scaf = ',num2str(InitTBP))  ,'HorizontalAlignment','left' );

Mat=[(1:N)',round(DistBtwTwoEnds/PitchConst), 2*ones(N,1),2*ones(N,1)];   %default 2x2
Mat(end+1,:)=[NaN,NaN,2,2];
Mat(:,end+1:end+5)=zeros(size(Mat,1),5);   %add 4 more column for [kZ1x,kZ1y ,kZ2x,kZ2y]
%    Mat(:,end+1)=zeros(size(Mat,1),4);   %add 4 more column for [kZ1x,kZ1y ,kZ2x,kZ2y]


Datat=mat2cell(Mat, ones(1,N+1 ) , [1,1,1,1,1,1,1,1,1]  );
for rev=1:size(Datat,1)
    Datat{rev,1}= num2str(Datat{rev,1});
end

str = cellfun(@num2str,num2cell([1:10,12,14]),'uniformoutput',0);
strHC=str;
% HCChoice={'X_SQ', '6HB', '10HB' ,'14HB','18HB','4HB','X_HC','SavedPart'};

Ext={'xSQ1';'xSQ2';'xSQ3';'xSQ4';'xSQ5';'6HB'; '10HB' ;'14HB';'18HB';'4HB';'48HB';'xHB1';'xHB2';'xHB3';'xHB4';'xHB5' ;'SavedPart' };
for kk=1: length(Ext)
    strHC{end+1}=Ext{kk}   ;
end
% HCChoice={'6HB', '10HB' ,'14HB','18HB'};

% TooltipStr = ['AAA ','BBB']   ;
s= '   ';
TooltipStr = 'This table is used to convert lines into bundles with corresponding geometries';
TooltipStr = [TooltipStr newline 'Edge : edge index related to line models in left-hand side'] ;
TooltipStr = [TooltipStr newline 'Length : length of bundles in cylinder direction, unit: base-pair'] ;
TooltipStr = [TooltipStr newline 'Nx/HC : Square lattice or honeycomb lattice :'] ;
TooltipStr = [TooltipStr newline 'Numeric: number of layers in x-direction.  #HB: default common HC shape'] ;
TooltipStr = [TooltipStr newline 'X_SQ and XHB are for cutomized cross-sections. Use buttons above to specify'] ;
TooltipStr = [TooltipStr newline 'SavedPart: assign one pre-saved geometries to all edges with this selection'] ;
TooltipStr = [TooltipStr newline 'Ny: number of layers in y-direction, only valid if Nx/HC is numberic'] ;
TooltipStr = [TooltipStr newline 'dZ1/dNx dZ1/dNy dZ2/dNx dZ2/dNy: gradience at both ends in X and Y directions, unit: base-pairs per 2nm spacing'] ;
TooltipStr = [TooltipStr newline 'Shift: initial orientaion as the base'] ;
TooltipStr= 'Temp. disable....' ;
% str
t = uitable(fH,'Units','normalized','Position',[0.7 0.17 0.25 0.45],....
    'ColumnWidth','auto','ColumnFormat',({[] [] strHC str}),...
    'ColumnEditable', true,'Data',Datat,....
    'ColumnName',{'Edge'; 'Length(BP)' ;'Nx/HC  ';'Ny';'dZ1/dNx';'dZ1/dNy';'dZ2/dNx';'dZ2/dNy' ;'shift'   } ,'Tag','STEPtable','Tooltip', TooltipStr);



sld2.Callback= @(src,evn)sldscale2(src,evn,txts2,DistBtwTwoEnds,txt2,t,txtTotalBP) ;
t.UserData.OriTableC2=round(DistBtwTwoEnds/PitchConst);
btnExport2.Callback=@(src,evn)CheckFinal(src,evn,t,sld1,sld2)  ;
BaColor =  repmat(t.BackgroundColor, 2+floor(N/2) ,1) ;
BaColor=BaColor(1:N+1,:) ;
t.UserData.oriBaColor = BaColor ;
t.BackgroundColor=BaColor ;

for extrai = 1 :length(ax.UserData.extraEdge)
    ax.UserData.extraEdge{extrai}.plotH.HitTest='on';
    edgeH{end+1} = ax.UserData.extraEdge{extrai}.plotH ;
end

for k=1:length(edgeH)
    edgeH{k}.ButtonDownFcn=@(src,evn)MapEdgeWithTable(src,evn,t,edgeH ) ;
end

ax.UserData.Incase_edgeH=edgeH ;

t.CellEditCallback=  @(src,evn)TablECK(src,evn,txtTotalBP,edgeH);

% t.CellSelectionCallback=@(src,evn)TablSelect(src,evn,edgeH);  %

btnRoundL= uicontrol(fH,'Style','pushbutton','Units','normalized',...
    'Position',[0.7 0.93 0.05 0.05],  'String','Round 10.5X');
btnRoundL.Callback=@(src,evn)RoundBP(src,evn,t)  ;

btnAssignCustomHC= uicontrol(fH,'Style','pushbutton','Units','normalized',...
    'Position',[0.7 0.85 0.05 0.05],  'String','AssignCustomHC');
btnAssignCustomHC.Callback=@(src,evn)AssignCustomHC ;

btnAssignCustomSQ= uicontrol(fH,'Style','pushbutton','Units','normalized',...
    'Position',[0.7 0.82 0.05 0.05],  'String','AssignCustomSQ','Tag','btnCustomSQ');   % AssignCustomSQ
btnAssignCustomSQ.Callback=@(src,evn)AssignCustomSQ ;
% Locations= strcat(pwd,filesep,'images\') ;
% % screensize = get( groot, 'Screensize' );
% btnAssignCustomSQ.Units='pixels' ;
% [a,~]=imread(strcat(Locations,'CustomSQ.jpg')) ; a=imresize(a, fliplr(btnAssignCustomSQ.Position(1,3:4)));
% a(1:2,:,:)=0 ; a(end-1:end,:,:)=0 ; a(:,1:2,:)=0 ; a(:,end-1:end,:)=0 ;
% btnAssignCustomSQ.Units='normalized'; btnAssignCustomSQ.CData= a ; btnAssignCustomSQ.String='';



btnLoadTable= uicontrol(fH,'Style','pushbutton','Units','normalized',...
    'Position',[0.9 0.82 0.05 0.05],  'String','LoadTableData' ,'Visible' ,'off');
btnLoadTable.Callback=@(src,evn)LoadTData(src,evn,t) ;

btnSaveTable= uicontrol(fH,'Style','pushbutton','Units','normalized',...
    'Position',[0.87 0.81 0.08 0.17],  'String','SaveTableData');
btnSaveTable.Callback=@(src,evn)saveTable(src,evn,t ,sld1,ax) ;

btnCalCurveLocal= uicontrol(fH,'Style','pushbutton','Units','normalized',...
    'Position',[0.78 0.81 0.08 0.17],  'String','CalCurve');
btnCalCurveLocal.Callback=@(src,evn)CalCurveLocal(src,evn , ax.UserData.LocalXYZ_h ,ax,t ,edgeH) ;





AssignIcon( btnAssignCustomHC,'CustomHC.jpg' ) ;  btnAssignCustomHC.TooltipString='Create customized crosssection for honeycomb lattice' ;
AssignIcon( btnAssignCustomSQ,'CustomSQ.jpg' ) ;  btnAssignCustomSQ.TooltipString='Create customized crosssection for square lattice' ;
AssignIcon( btnRoundL,'Round10p5.jpg' );  btnRoundL.TooltipString='Round to 10.5 bp' ;
AssignIcon( btnLoadTable,'LoadTable.jpg' ); btnLoadTable.TooltipString='Load the pre-saved table to recreate the same cylinder model. ' ;
AssignIcon( btnExport2,'ExportToAssembly.jpg' );  btnExport2.TooltipString='Convert the line model into the cylinder model, according to the table. ' ;
% AssignIcon( btnSaveTable,'SaveTable.jpg' ); btnSaveTable.TooltipString='Save current table. Overwrite the file "STEP_table_Data.mat" under current directory, which can be moved to somethere.' ;

AssignIcon( btnCalCurveLocal,'Extrude.jpg' ) ;  btnCalCurveLocal.TooltipString='Automatic calculate bundle orientations and edge gradients for the extrude method.' ;
AssignIcon( btnSaveTable,'Sweep.jpg' ); btnSaveTable.TooltipString='Call the sweep GUI to sweep along the spline path. Assign the number of breaks later.' ;


align([btnAssignCustomSQ , btnAssignCustomHC ,btnRoundL],'Left' ,'distribute' ) ; 
% align([btnAssignCustomSQ btnAssignCustomHC btnRoundL],'distribute','bottom');
align([btnCalCurveLocal btnSaveTable ],'distribute','bottom');

%             hLg= legend(ForLegend,'Click me for instructions','Location','northwest','Tag','SketchLG' ) ;
hLg= findobj(gctab,'Type','legend' )  ;
hLg.ButtonDownFcn=@(src,evn)LegendBoxing_SketchStep2( src,evn,ax ) ;
fH=gcf ; fH.UserData.MovieStep =3 ;
end % end of export bundle

function  CalCurveLocal(src,evn , LocalXYZ_h ,ax,t  ,edgeH)
fprintf('Calculating local coordinates, gradient values, shifts. \n')
PointsInfo = ax.UserData.Final.PointsInfo  ; % used for calculating the line graph
AdjLines = zeros( size(PointsInfo ,1) ,size(PointsInfo ,1) ) ;

for  i = 1 : size(PointsInfo ,1)
    for j = i+1 : size(PointsInfo ,1)
        FourPoints = [ PointsInfo(i,1:3) ; PointsInfo(i,4:6) ;PointsInfo(j,1:3) ;PointsInfo(j,4:6) ] ;
        Uni = unique( FourPoints ,'rows') ;
        if size(Uni,1) <4
            AdjLines(i,j)= 1 ; AdjLines(j,i)= 1 ;
        end
    end
end
%-----hard code
% AdjLines(1,4) =0 ;AdjLines(4,1) =0 ;
% AdjLines(1,5) =0 ;AdjLines(5,1) =0 ;
%


gAdjM2=graph(AdjLines);

% figure(445);clf; plot(gAdjM2) ; % debug use

CCmp =  conncomp(gAdjM2) ; % Instead of spanning tree alg, here develop similar tree based on depth-first-search by groups
predM2_CM = zeros(size(CCmp)) ; RootLine=[];
for k =1: max(CCmp)
    SubLines =  find(CCmp==k) ;
    AllPointsinSub =  [PointsInfo(SubLines,1:3) , PointsInfo(SubLines,4:6)] ;
    AllPointsinSubRR = [PointsInfo(SubLines,1:3) ;PointsInfo(SubLines,4:6)] ;
    CCMx = zeros(size(AllPointsinSub,1) ,2) ;
    for ll= 1:size(AllPointsinSub,1)
        Pt_L =  AllPointsinSub(ll,1:3) ;
        LIA = ismembertol(AllPointsinSubRR,Pt_L,0.01,'ByRows',true) ;   CCMx(ll,1) = sum(LIA)  ;
        Pt_R =  AllPointsinSub(ll,4:6) ;
        LIB = ismembertol(AllPointsinSubRR,Pt_R,0.01,'ByRows',true) ;   CCMx(ll,2) = sum(LIB)  ;
    end
    [U,~] = find(CCMx==1)  ;
    PossibleRoot = SubLines(U) ;
    if ~isempty(PossibleRoot)
        Rooot=PossibleRoot(1) ;
    elseif  sum(sum(CCMx==2))==numel(CCMx) % circle
        Rooot= SubLines(1);
    elseif max(max(CCMx))>2
        fprintf('Detect branch. Consider changing the line model. \n');
        return
    end
    dfs = dfsearch(gAdjM2,Rooot) ;
    predM2_CM( dfs(2:end))=dfs(1:end-1) ;
    RootLine= union(RootLine,Rooot) ;
end
predM2 =predM2_CM ;
% [T2,predM2] = minspantree(gAdjM2,'Method','sparse','Type','forest' ) ; predM2 ;
% RootLine = find(predM2 ==0)

fprintf('Number of root = %i\n',sum(predM2==0));
RootLine
CountBranch =  predM2 ;
CountBranch(CountBranch==0) = [] ;
BinEdge = 0.5:1:max(CountBranch)+1 ;
[N,edges] = histcounts(CountBranch,BinEdge) ;
if max(N) > 1
    fprintf('Detect branch. Consider changing the line model. \n');
    return
end

BundleOrientation = getCurrentOrientation(ax) ;   % [z , y ,x ]
LL_local = norm( diff([LocalXYZ_h{1,1}.XData;LocalXYZ_h{1,1}.YData;LocalXYZ_h{1,1}.ZData]' ) ) ;
LL_local = 5;
%hard code, make annimation
% writerObj = VideoWriter('Extude_annimation','MPEG-4');
% writerObj.Quality=100;
% writerObj.FrameRate = 1.5;  % Default 30
% open(writerObj);
%--------------

for k =1 :  length(RootLine)     %%local coordinate
    LastOne_twist = false ;
    
    for runtwice_for_twice = 1 :2
        if runtwice_for_twice==2 && LastOne_twist==false
            continue
        end
        
        CurrentLine = RootLine(k) ;
        Next = find(predM2 ==CurrentLine) ; cc= 0;
        while 1
            if isempty(Next)  % prevent single line
                break
            end
            if length(Next)>1 ;  fprintf('Detect branch. Consider changing the line model. \n');  return;   end
            %         Nvec = cross( BundleOrientation(CurrentLine,1:3 ) , BundleOrientation(Next,1:3 )  ) ;   Nvec= Nvec/norm(Nvec) ;
            FourPoints = [ PointsInfo(CurrentLine,1:3) ; PointsInfo(CurrentLine,4:6) ;PointsInfo(Next,1:3) ;PointsInfo(Next,4:6) ] ;
            [IndV,~] =ismembertol(FourPoints(1:2,:) ,  FourPoints(3:4,:) ,0.01 ,'Byrows',true) ;
            VertexThis=  FourPoints(IndV ,: ) ;
            isV1_outward =   dot( FourPoints(~IndV ,: )-VertexThis , BundleOrientation(CurrentLine,1:3 )) >0 ;
            [IndV2,~] =ismembertol(FourPoints(3:4,:) ,  FourPoints(1:2,:) ,0.01 ,'Byrows',true) ;
            isV2_outward =  dot( FourPoints([false ;false;~IndV2],: )-VertexThis , BundleOrientation(Next,1:3 )) >0 ;
            if ~xor(isV1_outward, isV2_outward )  % flip Z and X on Next
                EV.Button=1 ;
                ChangeLocalOrientation(LocalXYZ_h{Next ,1},EV,LocalXYZ_h) ;
                BundleOrientation = getCurrentOrientation(ax) ;   % update data
            end
            
            V1=  BundleOrientation(CurrentLine,1:3) ; V1=V1/norm(V1) ;  % local Z of bundle i
            V2=  BundleOrientation(Next,1:3) ;  V2=V2/norm(V2) ;   % local Z of bundle i+1
            Nvec = cross( V1 , V2  ) ;   Nvec= Nvec/norm(Nvec) ;
            theta = acos(dot(V1,V2 ) );
%             if imag(theta)~=0
%                 sddsf=3
%             end
            if theta~=0 && imag(theta)==0  % in case of parallel
                RMat = RotationAxisToRotaionMatrix( Nvec, theta*180/pi )   ;
            else 
                RMat=eye(3)   ;
            end
            
            if LastOne_twist==true % loop-runtwice_for_twice, add individual twist
                RMat_twist = RotationAxisToRotaionMatrix(   BundleOrientation(Next, 1:3), -IndivTwist )   ;
                RMat = RMat_twist*RMat ;
            end
            %         V1*RMat'  = V2 ;
            New_V2_y =  BundleOrientation(CurrentLine, 4:6)*RMat'  ;
            New_V2_x =  BundleOrientation(CurrentLine, 7:9)*RMat' ;
            MidPoint = ax.UserData.Final.PointsInfo(Next,7:9) ;
            y_XYZ = [ MidPoint ; MidPoint + LL_local*New_V2_y ]' ;
            x_XYZ = [ MidPoint ; MidPoint + LL_local*New_V2_x ]' ;
            
            LocalXYZ_h{Next,1}.XData=x_XYZ(1,:) ;  LocalXYZ_h{Next,1}.YData=x_XYZ(2,:) ;  LocalXYZ_h{Next,1}.ZData=x_XYZ(3,:) ;
            LocalXYZ_h{Next,2}.XData=y_XYZ(1,:) ;  LocalXYZ_h{Next,2}.YData=y_XYZ(2,:) ;  LocalXYZ_h{Next,2}.ZData=y_XYZ(3,:) ;
            BundleOrientation = getCurrentOrientation(ax) ;   % update data
            %---------
            NextNext = find(predM2 ==Next) ;
            if isempty(NextNext) || cc>300
                if AdjLines(RootLine(k),Next  ) >0  % compensate twist angle if closed-loop, Feb 28, 2021
                    LastOne_twist= true ;
                    CurrentLine=Next ;Next = RootLine(k) ;
                    
                    V1=  BundleOrientation(CurrentLine,1:3) ; V1=V1/norm(V1) ;  % local Z of bundle i
                    V2=  BundleOrientation(Next,1:3) ;  V2=V2/norm(V2) ;   % local Z of bundle i+1
                    Nvec = cross( V1 , V2  ) ;   Nvec= Nvec/norm(Nvec) ;
                    theta = acos(dot(V1,V2 ) );
                    if theta~=0  % in case of parallel
                        RMat_lastTwist = RotationAxisToRotaionMatrix( Nvec, theta*180/pi )   ;
                    else
                        RMat_lastTwist=eye(3)   ;
                    end
                    Proj_LastToFirst_V2_x =  BundleOrientation(CurrentLine, 7:9)*RMat_lastTwist' ;
                    First_fixed_x =  BundleOrientation(RootLine(k), 7:9) ;
                    First_fixed_y =  BundleOrientation(RootLine(k), 4:6) ;
                    
                    tan2x_onfirst = dot(Proj_LastToFirst_V2_x ,  First_fixed_x )  ;
                    tan2y_onfirst = dot(Proj_LastToFirst_V2_x ,  First_fixed_y )  ;
                    
                    closedLoop_twist = atan2d(tan2y_onfirst,tan2x_onfirst)
                    IndivTwist =  closedLoop_twist/(cc+2 )
                    % will be used in the next for loop-runtwice_for_twice
                    if abs(IndivTwist )<1 % set a threshold to avoid run second time
                        LastOne_twist= false ;
                    end
                end
%                 if ~isempty(Next)
%                     for kkk= 1 : length(edgeH)
%                         set(edgeH{kkk},'Color' , [0 0 0 ]) ;
%                         set(edgeH{kkk},'LineWidth' , 0.5) ;
%                     end
%                     set(edgeH{Next},'Color' , [0.9,0.9,0.6]) ;
%                     set(edgeH{CurrentLine},'Color' , [0.9,0.5,0.6]) ;
%                     set(edgeH{CurrentLine},'LineWidth' , 5) ;
%                     set(edgeH{Next},'LineWidth' , 5) ;
%                 end
%                 thisimage= getframe(gcf) ;
%                 writeVideo(writerObj, thisimage);                
                
                break
            else
%                 if ~isempty(Next)
%                     for kkk= 1 : length(edgeH)
%                         set(edgeH{kkk},'Color' , [0 0 0 ]) ;
%                         set(edgeH{kkk},'LineWidth' , 0.5) ;
%                     end
%                     set(edgeH{Next},'Color' , [0.9,0.9,0.6]) ;
%                     set(edgeH{CurrentLine},'Color' , [0.9,0.5,0.6]) ;
%                     set(edgeH{CurrentLine},'LineWidth' , 5) ;
%                     set(edgeH{Next},'LineWidth' , 5) ;
%                 end
%                 thisimage= getframe(gcf) ;
%                 writeVideo(writerObj, thisimage);
                
                
                CurrentLine=Next ;Next = find(predM2 ==CurrentLine) ;
            end
%                 edgeH{IndLine2}.Color=[0.2,0.2,0.7] ;edgeH{IndLine2}.LineWidth=3 ;

 
            
            cc=cc+1 ;
        end
    end
end
cc;
% close(writerObj);


for k =1: size(t.Data ,1)
    t.Data{k,5 } =0 ;
    t.Data{k,6 } =0 ;
    t.Data{k,7 } =0 ;
    t.Data{k,8 } =0 ;
    t.Data{k,9 } =0 ;
end


Const_BPperTurn = 48 ;
BundleOrientation = getCurrentOrientation(ax) ;   % [z , y ,x ]
ratiospacing = 0.34 /2.5 ; % 2.5nm spacing
% ratiospacing = 0.34 /3.2 ; % 3.2nm spacing

for k =1 :  length(RootLine)     %%local coordinate
    if length(t.Data{k,3})<3  || strcmp(t.Data{k,3}(1:end-1) , 'xSQ')   % bundle is square
        PeriodRound = 32 ;
    else
        PeriodRound = 21;
    end
    
    CurrentLine = RootLine(k) ;
    Next = find(predM2 ==CurrentLine) ; cc= 0; ConsiderRoot = 0 ;
    ShiftBase = 0;
    while 1
        if isempty(Next)  % prevent single line
            break
        end
        
        V1=  BundleOrientation(CurrentLine,1:3) ; V1=V1/norm(V1) ;
        V2=  BundleOrientation(Next,1:3) ;  V2=V2/norm(V2) ;
        Nvec = cross( V1 , V2  ) ;   Nvec= Nvec/norm(Nvec) ;
        
        alpha = acosd(dot(V1,V2)) ;
        
        if alpha<150
            Mag = [2*tand(alpha/2)/ratiospacing  ; 0  ]    ; % 2.5 nm spacing
        else
            Mag = [2*tand(150/2)/ratiospacing  ; 0  ] ;    % 2.5 nm spacing
        end
        %         Mag = [2*tand(alpha/2)/0.17  ; 0  ] ;   % 2nm spacing
        
        %          Mag = [alpha/360*Const_BPperTurn  ; 0  ] ;   % 48 Magic number
        
        
        V1_x = BundleOrientation(CurrentLine,7:9) ; V1_x=V1_x/norm(V1_x) ;
        V1_y = BundleOrientation(CurrentLine,4:6) ; V1_y=V1_y/norm(V1_y) ;
        Beta = acosd(dot(V1_x,Nvec)) ;
        if alpha~=0  % in case of parallel
            RMat = RotationAxisToRotaionMatrix( Nvec, alpha/2 )   ;
        else
            RMat=eye(3)   ;
        end
        NormalPlaneVec_z =  BundleOrientation(CurrentLine, 1:3)*RMat' ;
        NormalPlaneVec_z=NormalPlaneVec_z/ norm(NormalPlaneVec_z) ;
        %------
        %-------
        NormalPlaneVec_q = cross(Nvec ,NormalPlaneVec_z ) ;
        
        LongElipDir = NormalPlaneVec_q*inv(RMat' )  ;
        IniPhase = acosd( dot(LongElipDir,V1_x ))  ;
        IniPhase2 = acosd( dot(LongElipDir,V1_y ))  ;
        if IniPhase2>90
            IniPhase = -IniPhase ;
        end
        
        FourPoints = [ PointsInfo(CurrentLine,1:3) ; PointsInfo(CurrentLine,4:6) ;PointsInfo(Next,1:3) ;PointsInfo(Next,4:6) ] ;
        [IndV,~] =ismembertol(FourPoints(1:2,:) ,  FourPoints(3:4,:) ,0.01 ,'Byrows',true) ;
        VertexThis=  FourPoints(IndV ,: ) ;
        
        isV1_outward =   dot( FourPoints(~IndV ,: )-VertexThis , BundleOrientation(CurrentLine,1:3 )) >0 ;
        if isnan(IniPhase)
            RtM=eye(2)   ;
        else
            RtM = Rt(-IniPhase) ;
        end
        %         IniPhase
        if isV1_outward  % Current Z1 site + Next Z2 site
            Grads= RtM*Mag   ;
            t.Data{CurrentLine,5 } = -(Grads(1)/2 ) ;
            t.Data{CurrentLine,6 } = -(Grads(2)/2 ) ;
            %             t.Data{CurrentLine,5 } = -round(Grads(1)/2 ) ;
            %             t.Data{CurrentLine,6 } = -round(Grads(2)/2 ) ;
            
            t.Data{Next,7 } = (Grads(1)/2 ) ;
            t.Data{Next,8 } = (Grads(2)/2 ) ;
            %             t.Data{Next,7 } = round(Grads(1)/2 ) ;
            %             t.Data{Next,8 } = round(Grads(2)/2 ) ;
            
            %       ShiftBase ;
            if ConsiderRoot==0
                ShiftInThis =  t.Data{CurrentLine,9 } ;
                L_this = t.Data{Next,2 } ;
                NewShiftBase= mod( ShiftInThis - L_this ,PeriodRound) ;
                t.Data{Next,9 } = NewShiftBase-1  ;
            end
            %       sdfs=3
            
        else
            Grads= RtM*Mag   ;
            t.Data{CurrentLine,7 } = -(Grads(1)/2 ) ;
            t.Data{CurrentLine,8 } = -(Grads(2)/2 ) ;
            %             t.Data{CurrentLine,7 } = -round(Grads(1)/2 ) ;
            %             t.Data{CurrentLine,8 } = -round(Grads(2)/2 ) ;
            t.Data{Next,5 } = (Grads(1)/2 ) ;
            t.Data{Next,6 } = (Grads(2)/2 ) ;
            %             t.Data{Next,5 } = round(Grads(1)/2 ) ;
            %             t.Data{Next,6 } = round(Grads(2)/2 ) ;
            
            if ConsiderRoot==0
                ShiftInThis =  t.Data{CurrentLine,9 } ;
                L_this = t.Data{CurrentLine,2 } ;
                NewShiftBase= mod( ShiftInThis + L_this ,PeriodRound) ;
                t.Data{Next,9 } = NewShiftBase+1  ;
            end
            
        end
        %---------
        if   ConsiderRoot==1
            break;
        end
        
        NextNext = find(predM2 ==Next) ;
        if isempty(NextNext) || cc>300
            if AdjLines(RootLine(k),Next  ) >0
                ConsiderRoot=1 ;
                CurrentLine=Next ;Next = RootLine(k) ;
            else
                break
            end
        else
            CurrentLine=Next ;Next = find(predM2 ==CurrentLine) ;
        end
        
        cc=cc+1 ;
    end
    %     sdfs=3
end

fprintf('Finish Calculating local coordinates, gradient values, shifts. \n')

end

function M = Rt(theta)

M = [cosd(theta) , sind(theta) ; -sind(theta) , cosd(theta) ] ;

end



function saveTable(~,~,t , sld1,ax)
fH=gcf;
% get Xsection
if length(t.Data{1,3})<3  || strcmp(t.Data{1,3}(1:end-1) , 'xSQ')   % bundle is square
    if length(t.Data{1,3})<3
        Nx=str2double(t.Data{1,3}) ; Ny=str2double(t.Data{1,4});  %Ny=str2double(t.Data{iBun,4});
        if isnan( Nx); Nx=t.Data{1,3};end
        if isnan( Ny); Ny=t.Data{1,4};end
        Cyl1=[2,2] ;
        InplaneXY=zeros(Nx*Ny,2);        k=1;
        for i=1:Nx
            for j=1:Ny
                InplaneXY(k,:)= Cyl1+ (i-1)*[2,0] + (j-1)*[0,2] ;k=k+1;
            end
        end
    else % custom SQ
        [X,Y] = meshgrid(0:2:40,0:2:40) ;
        SQlattice.SQcenter= [X(:) ,Y(:) ] ;
        load('CustumSQLattice.mat','CustumSQLattice') ;
        SelectCase =  str2num( t.Data{1,3}(end) ) ;
        InplaneXY = SQlattice.SQcenter(CustumSQLattice{SelectCase} ,:) ;
    end
    XsecLattice =  'SQ' ; period =32 ;
    %      AddBundle=  BundleCylinderSQ(1,[],Z1Arr,Z2Arr,InplaneXY) ;
else    %  ith bundle is HC        if strcmp(t.Data{iBun,3},'6HB')
    if strcmp(t.Data{1,3},'6HB')
        [ InplaneXY ] = GiveHCinPlaneP( 1 ) ;    nCyl= 6 ;
    elseif strcmp(t.Data{1,3},'10HB')
        [ InplaneXY ] = GiveHCinPlaneP( 2 )  ;   nCyl= 10 ;
    elseif strcmp(t.Data{1,3},'14HB')
        [ InplaneXY ] = GiveHCinPlaneP( 3 )  ;   nCyl= 14 ;
    elseif strcmp(t.Data{1,3},'18HB')
        [ InplaneXY ] = GiveHCinPlaneP( 4 )  ;   nCyl= 18 ;
    elseif strcmp(t.Data{1,3},'4HB')
        [ InplaneXY ] = GiveHCinPlaneP( 5 )  ;   nCyl= 4 ;
    elseif strcmp(t.Data{1,3},'48HB')
        [ InplaneXY ] = GiveHCinPlaneP( 7 )  ;   nCyl= 48 ;
    elseif   strcmp(t.Data{1,3}(1:end-1) , 'xHB')
        SelectCase =  str2num( t.Data{1,3}(end) ) ;
        [ InplaneXY ] = GiveHCinPlaneP( 6 ,SelectCase)  ;   nCyl= size(InplaneXY,1) ;   % create custome crosssection.
    end
    XsecLattice =  'HC' ;
    period =21 ;
end
%%%
PointsInfo = ax.UserData.Final.PointsInfo  ; % used for calculating the line graph
AdjLines = zeros( size(PointsInfo ,1) ,size(PointsInfo ,1) ) ;

for  i = 1 : size(PointsInfo ,1)
    for j = i+1 : size(PointsInfo ,1)
        FourPoints = [ PointsInfo(i,1:3) ; PointsInfo(i,4:6) ;PointsInfo(j,1:3) ;PointsInfo(j,4:6) ] ;
        Uni = unique( FourPoints ,'rows') ;
        if size(Uni,1) <4
            AdjLines(i,j)= 1 ; AdjLines(j,i)= 1 ;
        end
    end
end

%
% ax=findobj(fH,'Tag','SS','Type','Axes');
VecAll = getCurrentOrientation(ax) ;
InplaneXY ;axes(ax) ;
ss_STEP= gctab ;    FFC= ss_STEP.UserData.FFC ;     L_bundle= cell2mat(t.Data(1:end-1,2)) ;

if isfield(ss_STEP.UserData,'New_tarr' )
    if ~isempty(ss_STEP.UserData.New_tarr)
        Cont_FFC_data=GUI_ContinuousSlicing( FFC,InplaneXY,VecAll, L_bundle,XsecLattice ,[],ss_STEP.UserData.New_tarr ,AdjLines)    ;
    else
        Cont_FFC_data=GUI_ContinuousSlicing( FFC,InplaneXY,VecAll, L_bundle,XsecLattice ,[],[] ,AdjLines)    ;
    end
else
    Cont_FFC_data=GUI_ContinuousSlicing( FFC,InplaneXY,VecAll, L_bundle,XsecLattice ,[],[],AdjLines)    ;
end


% % % % Cont_FFC_data= FFC.VisualizeFFC_noGUI_wNormalVector(InplaneXY,VecAll,L_bundle,XsecLattice) ;
%-------Save previous tarr for reload and adjust
Cont_FFC_data.Cylinders ;
New_tarr=Cont_FFC_data.Bundle.New_tarr ;
ss_STEP.UserData.New_tarr=New_tarr;
% return
% get data from freeformcurve class--------
BundleData= Cont_FFC_data.Bundle ;
Cylinders_Data= Cont_FFC_data.Cylinders;

Ltable = round(Cylinders_Data.ds_Cyli_ForBundle/0.34)  ; % convert to base
Ltable(Ltable<1) = 1;
%  IniBase = 34*ones(Cylinders_Data.N ,1) ;
Z1cumAvg = 34;
Midcum = 60 ;

% Z1siteRef = 55*ones(Cylinders_Data.N ,1 ) ;
for buni = 1:Cylinders_Data.Nbundles
    
    %     IniBase = Z1cumAvg*ones(Cylinders_Data.N ,1)  ;
    %     Z1Arr= IniBase  ;
    %     Z2Arr = IniBase +  Ltable(:,buni)-1 ;
    %     round(mean(Ltable(:,buni)))
    
    Z1Z2 = round([Midcum-0.5*Ltable(:,buni) , Midcum+0.5*Ltable(:,buni) ]) ; % for trapzoid bundle
    
%     Z1Z2 = round([Z1siteRef , Z1siteRef+1*Ltable(:,buni) ]) ;
    % for initializing next bundle, non-eq-trapzoid bundle
    
    
    while min(Z1Z2(:,1))< 30 %prevent negative
        Z1Z2 = Z1Z2+ period*ones(size(Z1Z2))     ;
    end
    
    if buni~=1
        while  ~ismember( [mod( round(mean(Z1Z2(:,1))), period )- mod(round(mean(SavePrevZ1Z2(:,2))), period) ] ,[1 , 1-period ])
            Z1Z2 = Z1Z2+ 1*ones(size(Z1Z2))     ;
        end
    end
    
    Z1Arr=Z1Z2(:,1) ;  Z2Arr=Z1Z2(:,2) ;
    
    %-------
%     Z1siteRef = Z1Z2(:,2)+1 ; % for initializing next bundle, non-eq-trapzoid bundle
%     Z1siteRef = Z1siteRef - mean(Z1siteRef)*ones(Cylinders_Data.N ,1 )+  (mod(mean(Z1siteRef), period) + 2*period)*ones(Cylinders_Data.N ,1 )  
    
    if strcmp(XsecLattice,'SQ')
        AddBundle=  BundleCylinderSQ(1,[],Z1Arr',Z2Arr',InplaneXY) ;
    else
        AddBundle=  BundleCylinderHC(1,[],Z1Arr',Z2Arr',InplaneXY) ;
    end
    
    if buni==1
        InportHyperBundle=hyperbundle(1,AddBundle); % first bundle
    else
        InportHyperBundle=AddBundles(InportHyperBundle,AddBundle);
    end
    SavePrevZ1Z2 = Z1Z2 ;
    
end
%-----------
VecAll= [ BundleData.TanVec , cross( BundleData.TanVec ,BundleData.NVec),BundleData.NVec];
%             QQ2 =cross(PTN_vec_onSpline(:,4:6 ), Temp ) ;
mag = sqrt(sum( BundleData.TanVec.^2 ,2) ) ;  mag= 1./mag;
VecAll(:,1:3) = BundleData.TanVec.*repmat(mag,1,3);

QQ2 = cross( BundleData.TanVec ,BundleData.NVec);
mag = sqrt(sum(QQ2.^2 ,2) ) ;  mag= 1./mag;
VecAll(:,4:6) = QQ2.*repmat(mag,1,3);

mag = sqrt(sum( BundleData.NVec.^2 ,2) ) ;  mag= 1./mag;
VecAll(:,7:9) = BundleData.NVec.*repmat(mag,1,3);


FFPI= sld1.Value * BundleData.Center;  % explode
% OFFPI=ax.UserData.Final.PointsInfo;

for bundi=1: length(InportHyperBundle.containBundle )
    %     bundi=UsedEdge(k) ;
    Bundle=InportHyperBundle.containBundle{bundi};
    TransR=[VecAll(bundi, 7:9)',VecAll(bundi, 4:6)' ,VecAll(bundi, 1:3)'];
    
    Bundle.TransformMatrix2(1:3,1:3)=TransR;
    Q1XYZG= Bundle.CylinderXYZGlobal;
    Q1Center= [mean([Q1XYZG(:,1);Q1XYZG(:,4)]),mean([Q1XYZG(:,2);Q1XYZG(:,5)]),mean([Q1XYZG(:,3);Q1XYZG(:,6)])];
    TargetD=FFPI(bundi,1:3);
    %     (TargetD-Q1Center)
    
    Bundle.TransformMatrix2(1:3,4)=(TargetD-Q1Center);
    InportHyperBundle.containBundle{bundi}=Bundle;
end

% prefercase = AskPrePairPreference();
prefercase.pair =1 ;
prefercase.NConn= round(size(InplaneXY,1)/2 ) ;
trial=0;
while trial<=1000
    [Doable,CellPairList,CellUnpair]=checkpairableH(InportHyperBundle,prefercase.pair);
    if  (Doable==1   ||  trial>=800)
        break;
    end
    trial=trial+1;
end
if Doable==1
    for k=2: length(CellPairList)
        CellPairList{k} =CellPairList{1} ;
    end
end
%----------
PremPair.Doable=Doable;PremPair.CellPairList=CellPairList;PremPair.CellUnpair=CellUnpair;
AdjM=zeros(Cylinders_Data.Nbundles ,Cylinders_Data.Nbundles  ) ;
for k =1:Cylinders_Data.Nbundles-1
    AdjM(k,k+1)=prefercase.NConn;  AdjM(k+1,k)=prefercase.NConn;
end

axes(findobj(fH,'Tag','AssemblyMain')) ;
cltab;
InportHyperBundle=InportHyperBundle.CreateTFWindow(PremPair,AdjM) ;

%------
TargetTab=findobj(fH,'Tag','ss_Assembly') ;
s_Mech= findobj(fH,'Tag','s_Mech') ;
s_Mech.Children.SelectedTab=TargetTab ;
%------------
%----------

for k=1: length(InportHyperBundle.containBundle)
    %     XvecData = [LocalXYZ_h{k ,1}.XData ; LocalXYZ_h{k ,1}.YData ;LocalXYZ_h{k ,1}.ZData;]  ;
    %     YvecData = [LocalXYZ_h{k ,2}.XData ; LocalXYZ_h{k ,2}.YData ;LocalXYZ_h{k,2}.ZData;]  ;
    %     ZvecData = [LocalXYZ_h{k ,3}.XData ; LocalXYZ_h{k ,3}.YData ;LocalXYZ_h{k,3}.ZData;]  ;
    Center = BundleData.Center(k,:) ; LL = 1 ;
    dz=VecAll(k,1:3)  ; dz=LL*dz/norm(dz) ;
    dy=VecAll(k,4:6)  ; dy=LL*dy/norm(dy) ;
    dx=VecAll(k,7:9)  ; dx=LL*dx/norm(dx) ;
    MM = [Center ; Center+dx ;Center ; Center+dy ;Center ; Center+dz ];
    
    InportHyperBundle.containBundle{k}.LocalCoordinatFromLineModel =MM';
end

fprintf('Finish exporting FFC. \n')

end


function LoadTData(src,evn,t)
src.Visible= 'off';
% [FileName,PathName] = uigetfile('*.mat','Select STEP data file');
%
% load(strcat(PathName,FileName),'STEP_table_Data') ;
% if size(t.Data ,1)== size(STEP_table_Data,1)
%     t.Data=STEP_table_Data ;
% else
%     fprintf('Data is not matched. Please select another one. \n')
% end
% %      sdfs=3;
end

%
% function TablSelect(src,evn,edgeH)
% evn
% evn.Indices;
% if isempty(evn.Indices); return;end
% SelectRow = evn.Indices(1) ;
% IndLine=SelectRow ;
% if IndLine == size(src.Data,1)
%     return
% end
% for k=1:length(edgeH)
%     edgeH{k}.Color=zeros(1,3) ;    edgeH{k}.LineWidth=0.5 ;
% end
% src.BackgroundColor=src.UserData.oriBaColor ;
% src.BackgroundColor(IndLine,:) = [0.9,0.9,0.6];
% edgeH{IndLine}.Color=[0.2,0.2,0.7] ;edgeH{IndLine}.LineWidth=3 ;
% end

function MapEdgeWithTable(src,evn,t ,edgeH)
IndLine= 0; TrackInds =[];
for k=1:length(edgeH)
    if isequal(edgeH{k},src);IndLine=k;end
    edgeH{k}.Color=zeros(1,3) ;    edgeH{k}.LineWidth=0.5 ;
    if strcmp(edgeH{k}.LineStyle,'-' )
        TrackInds=union(TrackInds,k);
    end
    
end
IndLine2=find(IndLine==TrackInds);
if ~isempty(IndLine2)
    t.BackgroundColor=t.UserData.oriBaColor ;
    t.BackgroundColor(IndLine2,:) = [0.9,0.9,0.6];
    edgeH{IndLine}.Color=[0.9,0.9,0.6] ;edgeH{IndLine}.LineWidth=2 ;
end
end


function RoundBP(src,evn,t)

for k=1:size(t.Data,1)-1
    if length(t.Data{k,3})<3  || strcmp(t.Data{k,3}(1:end-1) , 'xSQ')   % bundle is square
        PeriodRound = 10.67 ;
    else
        PeriodRound = 10.5 ;
    end
    OValue=t.Data{k,2} ;
    Nval= round(round(OValue/PeriodRound)*PeriodRound)-1;
    t.Data{k,2}=Nval ;
end

end

function TablECK(src,evn,txtTotalBP,edgeH)
HCChoice={'xSQ1';'xSQ2';'xSQ3';'xSQ4';'xSQ5';'6HB'; '10HB' ;'14HB';'18HB';'4HB';'48HB';'xHB1';'xHB2';'xHB3';'xHB4';'xHB5' ;'SavedPart' };

% HCChoice={'X_SQ', '6HB', '10HB' ,'14HB','18HB','4HB','48HB','X_HC','SavedPart'};
SQYChoice=cellfun(@num2str,num2cell([1:10,12,14]),'uniformoutput',0) ;
clc;
test=11111;
% src.ColumnFormat{4}

if strcmp(evn.EventName,'CellEdit')
    if evn.Indices(2)==4
        if strcmp(src.Data{evn.Indices(1),evn.Indices(2)-1},'HoneyComb')
            src.ColumnFormat{4}=HCChoice;
        else
            src.ColumnFormat{4}=SQYChoice;
        end
    end
end


if strcmp(evn.EventName,'CellEdit')
    if evn.Indices(1)==size(src.Data,1) %  change last row
        switch  evn.Indices(2)
            case 3
                if strcmp(evn.EditData, 'HoneyComb')
                    for k=1:size(src.Data,1)
                        src.Data{k,3}= 'HoneyComb' ;
                        %                     src.Data{k,4}= HCChoice ;
                    end
                else
                    for k=1:size(src.Data,1)
                        src.Data{k,3}=  evn.EditData;
                    end
                end
                
            case 4
                for k=1:size(src.Data,1)-1
                    src.Data{k,4}=  evn.EditData;
                end
        end
    else    %change individual bundle
        switch  evn.Indices(2)
            case 3
                if strcmp(evn.EditData, 'HoneyComb')
                    src.Data{evn.Indices(1),3}= 'HoneyComb' ;
                else
                    src.Data{evn.Indices(1),3}=  evn.EditData;
                end
            case 4
                if strcmp(src.Data{evn.Indices(1),evn.Indices(2)-1},'HoneyComb')
                    dsfsfsf=324 ;
                else
                end
                %
        end
    end
end


cc=0;
for k=1:size(src.Data,1)-1
    if length(src.Data{k,3})<3
        inthisbundle=src.Data{k,2}*str2num_CM(src.Data{k,3})*str2num_CM(src.Data{k,4});
    elseif strcmp(src.Data{k,3}(1:3) , 'xSQ')
        load('CustumSQLattice.mat','CustumSQLattice') ;
        Ind = str2num(src.Data{k,3}(end) ) ;
        inthisbundle=src.Data{k,2}*length(CustumSQLattice{Ind});
    else
        STR=src.Data{k,3} ;   STR=STR(1:end-2); Ind = str2num(src.Data{k,3}(end) ) ;
        
        inthisbundle=src.Data{k,2}*str2num_CM(STR) ;
        if isempty(str2num_CM(STR))
            load('CustumHCLattice.mat','CustumHCLattice') ;
            inthisbundle=src.Data{k,2}*length(CustumHCLattice{Ind});          % nCyl as in custom crosssection.  Case to case.
        end
    end
    cc=cc+inthisbundle;
end
txtTotalBP.String=strcat('Ap.Scaf = ',num2str(cc));

if strcmp(evn.EventName,'CellEdit')
    evn.Indices;
    if isempty(evn.Indices); return;end
    SelectRow = evn.Indices(1) ;
    IndLine=SelectRow ;
    if IndLine == size(src.Data,1)
        return
    end
    TrackInds=[];
    for k=1:length(edgeH)
        edgeH{k}.Color=zeros(1,3) ;    edgeH{k}.LineWidth=0.5 ;
        if strcmp(edgeH{k}.LineStyle,'-' )
            TrackInds=union(TrackInds,k);
        end
    end
    IndLine2=TrackInds(IndLine);
    src.BackgroundColor=src.UserData.oriBaColor ;
    src.BackgroundColor(IndLine,:) = [0.9,0.9,0.6];
    edgeH{IndLine2}.Color=[0.2,0.2,0.7] ;edgeH{IndLine2}.LineWidth=3 ;
end

end

function CheckFinal(src,evn,t,sld1,sld2)
% fH=gcf;
% ax=gca;
ax=findobj(src.Parent.Children,'Type','Axes') ;
%--------check each bundle has even cylinders
N_Edge=size(t.Data,1)-1 ;
EffecXY=zeros(N_Edge,1);
for k=1:N_Edge
    if length(t.Data{k,3})==1
        if mod(t.Data{k,3},2)==0 || mod(t.Data{k,4},2)==0
            EffecXY(k)=1;
        end
    else
        EffecXY(k)=1;
    end
end

if nnz(EffecXY) ~=N_Edge
    warndlg('Some Edge has odd cylinder','!! Warning !!')
    return
end
fH=gcf;
fbar = waitbar(0,'Please wait...');
pause(.5)
%-----------------
RAll=sld1.Value*sld2.Value;
% figure(223); clf;hold on; axis equal;
% ax2=gca;
% scatter3( [FFPI(:,1);FFPI(:,4);FFPI(:,7)] ,[FFPI(:,2);FFPI(:,5);FFPI(:,8)],[FFPI(:,3);FFPI(:,6);FFPI(:,9)] );
% scatter3( [OFFPI(:,1);OFFPI(:,4);OFFPI(:,7)] ,[OFFPI(:,2);OFFPI(:,5);OFFPI(:,8)],[OFFPI(:,3);OFFPI(:,6);OFFPI(:,9)] );
% NX1=str2num_CM(t.Data{1,3});   NY1=str2num_CM(t.Data{1,4});
% kZ1x=round(t.Data{1,5}); kZ1y=round(t.Data{1,6}); kZ2x=round(t.Data{1,7}); kZ2y=round(t.Data{1,8});
% TotalBPLength=t.Data{1,2} ;

ChoseLoadSavedPart=0 ;
for k=1:size(t.Data,1)-1
    if  strcmp(t.Data{k,3},'SavedPart'  )
        ChoseLoadSavedPart = 1 ;
    end
end

if ChoseLoadSavedPart==1  % Ask
    PartFolder=[ pwd '\parts\'];
    [FileName,PathName] = uigetfile('*.mat','DialogTitle',PartFolder);
    SSS=open(strcat(PathName,FileName));
    LoadingPart=SSS.Psave.GUISavePart ;
    if isempty(LoadingPart)  % for Hub bundle
        LoadingPart=SSS.Psave ;
    end
end


waitbar(.2,fbar,'Loading your data');
InportHyperBundle=[];
N_ExtraEdge=length(ax.UserData.ExtarEdgeRemain);
first = 1 ; UsedEdge = [] ;
for iBun=1:N_Edge
    
    %     kZ1x=(t.Data{iBun,5}); kZ1y=(t.Data{iBun,6}); kZ2x=(t.Data{iBun,7}); kZ2y=(t.Data{iBun,8});
    kZ1x=round(t.Data{iBun,5}); kZ1y=round(t.Data{iBun,6}); kZ2x=round(t.Data{iBun,7}); kZ2y=round(t.Data{iBun,8});
    
    % slopes for any bundle, including SQ and HC
    TotalBPLength=t.Data{iBun,2} ;
    %     if  TotalBPLength<5
    %         continue;   % allow ignore certain edges/bundle shorter than 5
    %     end
    
    if length(t.Data{iBun,3})<3  || strcmp(t.Data{iBun,3}(1:end-1) , 'xSQ')   % bundle is square
        if length(t.Data{iBun,3})<3
            Nx=str2double(t.Data{iBun,3}) ; Ny=str2double(t.Data{iBun,4});%  Ny=str2double(t.Data{iBun,4});
            
            if isnan( Nx); Nx=t.Data{iBun,3};end
            if isnan( Ny); Ny=t.Data{iBun,4};end
            if iBun<=N_Edge-N_ExtraEdge
                shiftZ1=34 + round(t.Data{iBun,9})   ;
            else
                shiftZ1=39 + round(t.Data{iBun,9})  ;
            end
            Cyl1=[2,2] ;
            InplaneXY=zeros(Nx*Ny,2);
            k=1;
            for i=1:Nx
                for j=1:Ny
                    InplaneXY(k,:)= Cyl1+ (i-1)*[2,0] + (j-1)*[0,2] ;k=k+1;
                end
            end
            %       HardCode=[1,1,2]
            %               Z1ArrT=shiftZ1*ones(NX1,NY1)+kZ1x*[0:NX1-1]'*ones(1,NY1)+kZ1y*ones(NX1,1)  ;
            %               Z2ArrT=(shiftZ1+TotalBPLength)*ones(NX1,NY1)+kZ2x*[0:NX1-1]'*ones(1,NY1)+kZ2y*ones(NX1,1)  ;
            Z1ArrT=shiftZ1*ones(Nx,Ny)+kZ1x*[0:Nx-1]'*ones(1,Ny) + kZ1y*ones(Nx,1)*[0:Ny-1]  ;    % linear
            Z2ArrT=(shiftZ1+TotalBPLength)*ones(Nx,Ny)+kZ2x*[0:Nx-1]'*ones(1,Ny)+kZ2y*ones(Nx,1)*[0:Ny-1]  ;
            
%             if mod(iBun,2) ==0  % hard code
%                 Z2ArrT(3:4,:) = Z2ArrT(3:4,:)-30 ;
%             else
%                 Z1ArrT(3:4,:) = Z1ArrT(3:4,:)+30 ;
%             end
            
            
        else % custom SQ
            if iBun<=N_Edge-N_ExtraEdge
                shiftZ1=34 + round(t.Data{iBun,9})   ;
            else
                shiftZ1=39 + round(t.Data{iBun,9})  ;
            end
            Cyl1=[2,2] ;
            [X,Y] = meshgrid(0:2:40,0:2:40) ;
            SQlattice.SQcenter= [X(:) ,Y(:) ] ;
            load('CustumSQLattice.mat','CustumSQLattice') ;
            SelectCase =  str2num( t.Data{iBun,3}(end) ) ;
            %                 HCIndex{6}=CustumHCLattice ;
            InplaneXY = SQlattice.SQcenter(CustumSQLattice{SelectCase} ,:) ;
            %             InplaneXY= round( InplaneXY -  min(InplaneXY) + Cyl1 ) ;
            RefCyl=InplaneXY(1,:) ;
            %             Z1ArrT=shiftZ1*ones(size(InplaneXY,1),1) +  kZ1x*  [0:Nx-1]'*ones(1,Ny) + kZ1y*ones(Nx,1)*[0:Ny-1]  ;    % linear
            Z1ArrT=shiftZ1*ones(size(InplaneXY,1),1) +  kZ1x/2*(InplaneXY(:,1)- RefCyl(1)) + kZ1y/2*(InplaneXY(:,2)- RefCyl(2)) ;    % linear
            Z2ArrT=(shiftZ1+TotalBPLength)*ones(size(InplaneXY,1),1) +  kZ2x/2*(InplaneXY(:,1)- RefCyl(1))  + kZ2y/2*(InplaneXY(:,2)- RefCyl(2))  ;
            
        end
        %--------------------------------
        Z1ArrT= round(Z1ArrT- mean(mean(Z1ArrT))+ shiftZ1  ) ;     %09042018
        Z2ArrT= round(Z2ArrT- mean(mean(Z2ArrT))+ (shiftZ1+TotalBPLength) )  ;  %09042018
        
        minZ1= min(min(Z1ArrT)) ;
        ShiftToPositve = 32*ceil((shiftZ1-minZ1)/32) ;
        Z1ArrT=Z1ArrT+   ShiftToPositve    ;
        Z2ArrT=Z2ArrT+   ShiftToPositve    ;
        
        Z1Arr=reshape(Z1ArrT',[1,size(Z1ArrT,1)*size(Z1ArrT,2)]);
        Z2Arr=reshape(Z2ArrT',[1,size(Z2ArrT,1)*size(Z2ArrT,2)]);
        
        AddBundle=  BundleCylinderSQ(1,[],Z1Arr,Z2Arr,InplaneXY) ;
    elseif  strcmp(t.Data{iBun,3},'SavedPart')     % ChoseLoadSavedPart==1
        
        if  isa(LoadingPart,'BundleCylinderSQ')      %if SQ
            Part1=BundleCylinderSQ(1,[], LoadingPart.Zbase1  ,LoadingPart.Zbase2 ,LoadingPart.CylInplanePosition)   ;
        elseif  isa(LoadingPart ,'BundleCylinderHC')       %if HC
            Part1=BundleCylinderHC(1,[], LoadingPart.Zbase1  ,LoadingPart.Zbase2 ,LoadingPart.CylInplanePosition)   ;
        end
        
        AddBundle= Part1 ;
    else    %  ith bundle is HC
        %         hello=1
        
        if strcmp(t.Data{iBun,3},'6HB')
            [ XY ] = GiveHCinPlaneP( 1 ) ;    nCyl= 6 ;
        elseif strcmp(t.Data{iBun,3},'10HB')
            [ XY ] = GiveHCinPlaneP( 2 )  ;   nCyl= 10 ;
        elseif strcmp(t.Data{iBun,3},'14HB')
            [ XY ] = GiveHCinPlaneP( 3 )  ;   nCyl= 14 ;
        elseif strcmp(t.Data{iBun,3},'18HB')
            [ XY ] = GiveHCinPlaneP( 4 )  ;   nCyl= 18 ;
        elseif strcmp(t.Data{iBun,3},'4HB')
            [ XY ] = GiveHCinPlaneP( 5 )  ;   nCyl= 4 ;
        elseif strcmp(t.Data{iBun,3},'48HB')
            [ XY ] = GiveHCinPlaneP( 7 )  ;   nCyl= 48 ;
        elseif   strcmp(t.Data{iBun,3}(1:end-1) , 'xHB')
            SelectCase =  str2num( t.Data{iBun,3}(end) ) ;
            [ XY ] = GiveHCinPlaneP( 6 ,SelectCase)  ;   nCyl= size(XY,1) ;   % create custome crosssection.
        end
        %   figure; scatter(XY(:,1) ,XY(:,2) ); axis equal  ;  % if want to  see crossection
        
        shiftZ1=24 + round(t.Data{iBun,9})  ;
        Z1Arr=shiftZ1*ones(1, nCyl);   % initial values without gradient
        Z2Arr=(t.Data{iBun,2}+shiftZ1)*ones(1, nCyl); % initial values without gradient
        
        RefCyl=XY(1,:) ;
        GradZ1 =   round(kZ1x*( (XY(:,1)-RefCyl(1) )/ 2)') +   floor(kZ1y*((XY(:,2)-RefCyl(2) )/2 )' ) ;
        %         floor(3*(((XY(:,2)-RefCyl(2) )/2 ))')
        GradZ2 =   round(kZ2x*((XY(:,1)-RefCyl(1) )/ 2)') +   floor(kZ2y*((XY(:,2)-RefCyl(2) )/2 )' ) ;
        %         GradZ1 =   round(kZ1x*( (XY(:,1)-RefCyl(1) )/ sqrt(3))') +   floor(kZ1y*((XY(:,2)-RefCyl(2) )/2 )' ) ;
        %         GradZ2 =   round(kZ2x*((XY(:,1)-RefCyl(1) )/ sqrt(3))') +   floor(kZ2y*((XY(:,2)-RefCyl(2) )/2 )' ) ;
        
        
        Z1Arr=Z1Arr+ GradZ1 ; Z2Arr=Z2Arr+ GradZ2 ;  % considering  gradient
        
        Z1Arr= round(Z1Arr- mean(mean(Z1Arr))+ shiftZ1  ) ;   % shift for orientation    %09042018
        Z2Arr= round(Z2Arr- mean(mean(Z2Arr))+ (shiftZ1+TotalBPLength) )  ;  %09042018
        %        end
        minZ1= min(min(Z1Arr)) ;  ShiftToPositve = 21*ceil((shiftZ1-minZ1)/21)+21   ;  % shift positive
        Z1Arr=Z1Arr+   ShiftToPositve    ;
        Z2Arr=Z2Arr+   ShiftToPositve    ;
        
        %           if ismember(iBun , [9:12])   % hard code
        %                 Z1Arr(1:2)=Z1Arr(1:2) +21  ;
        %           elseif  ismember(iBun , [2:4])   % hard code
        %                Z1Arr(1:2)=Z1Arr(1:2) -11  ;
        %                Z2Arr(1:2)=Z2Arr(1:2) +11  ;
        %           end
        
        %         Z1Arr(1:2) = Z1Arr(1:2) +6 ;  % hard
        %         if ~ismember(iBun,[9:12])
        %         Z2Arr(1:2) = Z2Arr(1:2) -6 ;
        %         end
        AddBundle=  BundleCylinderHC(1,[],Z1Arr,Z2Arr,XY) ;
        %         AddBundle.Tol=5 ;  % hard
    end
    
    if first==1
        InportHyperBundle=hyperbundle(1,AddBundle); first=0 ; % first bundle
        UsedEdge=union(UsedEdge,iBun) ;
    else
        InportHyperBundle=AddBundles(InportHyperBundle,AddBundle);
        UsedEdge=union(UsedEdge,iBun) ;
    end
    
end
waitbar(.3,fbar,'Moving bundles');

% ax.UserData.Final.PointsInfo;
% VecAll=ax.UserData.Final.VecInfo;   % default orientation of bundles
BundleOrientation = getCurrentOrientation(ax) ;
VecAll = BundleOrientation ;
FFPI= RAll*(ax.UserData.Final.PointsInfo);
OFFPI=ax.UserData.Final.PointsInfo;

for k=1: length(InportHyperBundle.containBundle )
    bundi=UsedEdge(k) ;
    Bundle=InportHyperBundle.containBundle{k};
    TransR=[VecAll(bundi, 7:9)',VecAll(bundi, 4:6)' ,VecAll(bundi, 1:3)'];
    
    Bundle.TransformMatrix2(1:3,1:3)=TransR;
    Q1XYZG= Bundle.CylinderXYZGlobal;
    Q1Center= [mean([Q1XYZG(:,1);Q1XYZG(:,4)]),mean([Q1XYZG(:,2);Q1XYZG(:,5)]),mean([Q1XYZG(:,3);Q1XYZG(:,6)])];
    TargetD=FFPI(bundi,7:9);
    
    Bundle.TransformMatrix2(1:3,4)=(TargetD-Q1Center);
    InportHyperBundle.containBundle{k}=Bundle;
    %         if bundi==18
    %         sdfsf=3
    %     end
    
end


prefercase = AskPrePairPreference();


waitbar(.5,fbar,'checking pairing');

trial=0;
while trial<=1000
    [Doable,CellPairList,CellUnpair]=checkpairableH(InportHyperBundle,prefercase.pair);
    if  (Doable==1   ||  trial>=800)
        break;
    end
    trial=trial+1;
end

% fprintf('Using hard code for pairing !! \n')
% GroupPair=setdiff(2:15,4:6);% hard code
% for k=1:length(GroupPair)    % hard code
%    CellPairList{GroupPair(k)}=CellPairList{1} ;
% end
% CellPairList{2}=CellPairList{1} ;CellPairList{3}=CellPairList{1} ;  % hard code
% CellPairList{4}=CellPairList{1} ;
for k =1:length(CellPairList)   % hard  code
CellPairList{k} =CellPairList{1} ; 
end

PremPair.Doable=Doable;
PremPair.CellPairList=CellPairList;
PremPair.CellUnpair=CellUnpair;
% PremPair
OriNodes=unique([OFFPI(:,1:3);OFFPI(:,4:6)],'rows') ;
AdjM=zeros(N_Edge,N_Edge);
nConn=prefercase.NConn;    % number of connection between bundle---------------
% cm=0;
for k=1:size(OriNodes,1)   %  for adjM of orininal edge
    Node=OriNodes(k,:);
    AttacheTo1= find(ismember(OFFPI(:,1:3),Node,'rows'));  %Edge index
    AttacheTo2= find(ismember(OFFPI(:,4:6),Node,'rows'));   %Edge index
    
    OnEdgePoints=[];
    for oldedge=1:length(ax.UserData.edgeHhandle)
        if isfield(ax.UserData.edgeHhandle{oldedge}.UserData,'InsertNode')
            if ismember(Node,ax.UserData.edgeHhandle{oldedge}.UserData.InsertNode,'rows')
                OnEdgePoints=union( ax.UserData.RelateTable(oldedge,2),OnEdgePoints);
            end
        end
    end
    
    Pooled=union(AttacheTo1,AttacheTo2);
    Pooled=union(OnEdgePoints,Pooled);
    Pooled=setdiff(Pooled,0);
    length(Pooled);
    if length(Pooled)>1
        for kk2=1:length(Pooled)-1
            AdjM(Pooled(kk2),Pooled(kk2+1))=nConn;
            AdjM(Pooled(kk2+1),Pooled(kk2))=nConn;
        end
        AdjM(Pooled(1),Pooled(end))=nConn;   % circular links, optional
        AdjM(Pooled(end),Pooled(1))=nConn;
    elseif length(Pooled)==1
        for otheredge=1:size(OFFPI,1)
            if otheredge~=Pooled
                Mat=[OFFPI(otheredge,1:3);OFFPI(otheredge,4:6);Node];
                
                ss= dot(( Mat(3,1:3)-Mat(1,1:3) ) ,( Mat(2,1:3)-Mat(1,1:3) ))/norm(Mat(2,1:3)-Mat(1,1:3) )/norm(Mat(2,1:3)-Mat(1,1:3) ) ;  % see if interpolate or extrapolate
                area = norm(cross(Mat(2,:)-Mat(1,:),Mat(3,:)-Mat(1,:)));
                if area<=1e-6  && (ss>=0  && ss<=1 )
                    AdjM(Pooled,otheredge)=nConn;
                    AdjM(otheredge,Pooled)=nConn;
                    %                 cm=cm+1;
                    %                 if cm==5
                    %                     sdfsf=2
                    %                 end
                    %
                    %                 [otheredge,Pooled,cm]
                end
            end
        end
    end
end
AdjM=AdjM(UsedEdge, UsedEdge) ;
% [u,v]=find(AdjM); [u,v]

N_ExtraEdge=length(ax.UserData.ExtarEdgeRemain);
N_OriEdge=N_Edge-N_ExtraEdge ;


% end
% clc;

axes(findobj(fH,'Tag','AssemblyMain')) ;
cltab;

% %------------hard code
% try
%     PrevSave=SaveLoadTransMandConn(0  ) ;
%     if length(PrevSave.SaveTransM) ==length(InportHyperBundle.containBundle)
%         AdjM =PrevSave.SaveAdjM;
%
%         for k=1 :length(InportHyperBundle.containBundle)
%             InportHyperBundle.containBundle{k}.TransformMatrix2 = PrevSave.SaveTransM{k} ;
%         end
%         fprintf(' number of bundle matches Temp.mat, use saved data \n')
%     end
% end
%---------------------  for changing the way to solidify the step file into
% cylinder model, if one trial design has been finished about
% transformation matrix and adjacent matrix, Use
% "SaveLoadTransMandConn(1  )" to save parameter. Later re-cylinderization
% step will check number of bundle.
waitbar(0.8,fbar,'Creating Mechanism');


InportHyperBundle=InportHyperBundle.CreateTFWindow(PremPair,AdjM) ;

waitbar(1,fbar,'finishing');

TargetTab=findobj(fH,'Tag','ss_Assembly') ;
% tt=gctab;
s_Mech= findobj(fH,'Tag','s_Mech') ;
% tt=findobj(fH,'Tag','ss_Assembly') ;
s_Mech.Children.SelectedTab=TargetTab ;
close(fbar) ;
%----------
LocalXYZ_h=ax.UserData.LocalXYZ_h ;
for k=1: length(InportHyperBundle.containBundle)
    XvecData = [LocalXYZ_h{k ,1}.XData ; LocalXYZ_h{k ,1}.YData ;LocalXYZ_h{k ,1}.ZData;]  ;
    YvecData = [LocalXYZ_h{k ,2}.XData ; LocalXYZ_h{k ,2}.YData ;LocalXYZ_h{k,2}.ZData;]  ;
    ZvecData = [LocalXYZ_h{k ,3}.XData ; LocalXYZ_h{k ,3}.YData ;LocalXYZ_h{k,3}.ZData;]  ;
    
    
    InportHyperBundle.containBundle{k}.LocalCoordinatFromLineModel =[XvecData,YvecData,ZvecData];
end

end

function  BundleOrientation=getCurrentOrientation(CurretnAxes)
% BundleOrientation =  edge_ZYX_Vec , nrows for edges , dim=[n,9]
BundleOrientation=zeros(size(CurretnAxes.UserData.LocalXYZ_h,1)  , 9 ) ;

for Buni= 1:size(BundleOrientation,1)
    xVec = [diff(CurretnAxes.UserData.LocalXYZ_h{Buni,1}.XData ) ,  diff(CurretnAxes.UserData.LocalXYZ_h{Buni,1}.YData ),diff(CurretnAxes.UserData.LocalXYZ_h{Buni,1}.ZData )] ;
    yVec = [diff(CurretnAxes.UserData.LocalXYZ_h{Buni,2}.XData ) ,  diff(CurretnAxes.UserData.LocalXYZ_h{Buni,2}.YData ),diff(CurretnAxes.UserData.LocalXYZ_h{Buni,2}.ZData )] ;
    zVec = [diff(CurretnAxes.UserData.LocalXYZ_h{Buni,3}.XData ) ,  diff(CurretnAxes.UserData.LocalXYZ_h{Buni,3}.YData ),diff(CurretnAxes.UserData.LocalXYZ_h{Buni,3}.ZData )] ;
    xVec=xVec/norm(xVec) ;    yVec=yVec/norm(yVec) ;    zVec=zVec/norm(zVec) ;
    BundleOrientation(Buni,:) =[ zVec,yVec , xVec] ;
    
    
end



% BundleOrientation=0
end

function sldscale1(src,evn,txts1)
txts1.String=strcat('S_Blow = ' ,num2str(src.Value));
end

function sldscale2(src,evn,txts2,DistBtwTwoEnds,txt2,t,txtTotalBP)
txts2.String=strcat('S_Enlarge = ' ,num2str(src.Value));
mmnanometer= round(min(src.Value*DistBtwTwoEnds)*10)/10 ;
mmbase= round(mmnanometer/0.34) ;

txt2.String=strcat('After enlarge =>  [', num2str(mmnanometer), '(nm)',{' '}, num2str(mmbase),' (bp) ]');
TD=t.Data;
t.UserData.OriTableC2;
for k=1:size(TD,1)-1
    TD{k,2}=round(t.UserData.OriTableC2(k)*src.Value);
end
t.Data=TD;
% sdf=123
TablECK(t,evn,txtTotalBP);
end




function AddTruss(src,evn,edgeH,edge_ZYX_Vec)
ax=gca;
nExtraEdge= length(ax.UserData.extraEdge);
XYZ=[ax.UserData.singleGreenH.XData, ax.UserData.singleGreenH.YData,ax.UserData.singleGreenH.ZData,...
    ax.UserData.singleRedH.XData, ax.UserData.singleRedH.YData, ax.UserData.singleRedH.ZData];

pH=plot3([XYZ(:,1);XYZ(:,4)]' ,[XYZ(:,2);XYZ(:,5)]',[XYZ(:,3);XYZ(:,6)]' );
pH.HitTest='off';pH.Color=[0.2,0.2,0.2];
EE=[];
for edgi=1:length(edgeH)
    %    if sum(edgeH{edgi}.Color~=[0.2 ,0.2,0.2])==3
    if sum(edgeH{edgi}.Color==[1,0,0])==3
        EE=union(EE,edgi);
        if isfield(edgeH{edgi}.UserData,'InsertNode')
            edgeH{edgi}.UserData.InsertNode=[edgeH{edgi}.UserData.InsertNode; XYZ(4:6)];
        else
            edgeH{edgi}.UserData.InsertNode=    XYZ(4:6);
        end
        
    elseif sum(edgeH{edgi}.Color==[0,1,0])==3
        EE=union(EE,edgi);
        
        if isfield(edgeH{edgi}.UserData,'InsertNode')
            edgeH{edgi}.UserData.InsertNode=[edgeH{edgi}.UserData.InsertNode; XYZ(1:3)];
        else
            edgeH{edgi}.UserData.InsertNode=    XYZ(1:3);
        end
    end
    
    
    
end

ZVec=XYZ(4:6)-XYZ(1:3) ; ZVec=ZVec/norm(ZVec);
Za=edge_ZYX_Vec(EE(1),1:3);
Zb=edge_ZYX_Vec(EE(2),1:3);
if norm( cross(Za,Zb))==0  %parallel
    Ytemp=Za;
else
    Ytemp=0.5*(Za+Zb);
end
YVec=cross(cross(  Ytemp,  ZVec),ZVec);YVec=YVec/norm(YVec);
XVec=cross(YVec,ZVec);XVec=XVec/norm(XVec);

G=[ mean([XYZ(1),XYZ(4)]),mean([XYZ(2),XYZ(5)]),mean([XYZ(3),XYZ(6)])];
V=[ZVec,YVec,XVec];
Zline=[G ; G+V(1:3)] ; Yline=[G ; G+V(4:6)] ;Xline=[G ; G+V(7:9)] ;
pZ=plot3(Zline(:,1),Zline(:,2),Zline(:,3),'b','LineWidth',4);% pZ.HitTest='off';
pY=plot3(Yline(:,1),Yline(:,2),Yline(:,3),'g','LineWidth',4); %pY.HitTest='off';
pX=plot3(Xline(:,1),Xline(:,2),Xline(:,3),'r','LineWidth',4); %pX.HitTest='off';

% LocalXYZ_h= ax.UserData.LocalXYZ_h  ;  % dynamic orientation data
ax.UserData.LocalXYZ_h{end+1,1}=pX ;
ax.UserData.LocalXYZ_h{end,2}=pY ;
ax.UserData.LocalXYZ_h{end,3}=pZ ;
nEdge =size(ax.UserData.LocalXYZ_h ,1) ;
ax.UserData.LocalXYZ_h{nEdge,1}.UserData.Edgei_XYZ = [nEdge,1] ;
ax.UserData.LocalXYZ_h{nEdge,2}.UserData.Edgei_XYZ = [nEdge,2] ;
ax.UserData.LocalXYZ_h{nEdge,3}.UserData.Edgei_XYZ = [nEdge,3] ;

for k=1:nEdge
    ax.UserData.LocalXYZ_h{k,1}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,ax.UserData.LocalXYZ_h) ;
    ax.UserData.LocalXYZ_h{k,2}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,ax.UserData.LocalXYZ_h) ;
    ax.UserData.LocalXYZ_h{k,3}.ButtonDownFcn =@(src,evn)ChangeLocalOrientation(src,evn,ax.UserData.LocalXYZ_h) ;
    
end

%     LocalXYZ_h{ke,1}=plot3(Xline(:,1),Xline(:,2),Xline(:,3),'r');  LocalXYZ_h{ke,1}.LineWidth=4; % LocalXYZ_h{ke,1}.HitTest='off';


ax.UserData.extraEdge{end+1}.connect2Edge=EE;
ax.UserData.extraEdge{end}.XYZ=XYZ;
ax.UserData.extraEdge{end}.plotH=pH;
ax.UserData.extraEdge{end}.midP=[ mean([XYZ(1),XYZ(4)]),mean([XYZ(2),XYZ(5)]),mean([XYZ(3),XYZ(6)])];
ax.UserData.extraEdge{end}.edge_ZYX_Vec=[ZVec,YVec,XVec];
ax.UserData.extraEdge{end}.plotHx=pX;
ax.UserData.extraEdge{end}.plotHy=pY;
ax.UserData.extraEdge{end}.plotHz=pZ;

% ax.UserData.extraEdge{end}.scaCenter=scatter3(mean([XYZ(1),XYZ(4)]),mean([XYZ(2),XYZ(5)]),mean([XYZ(3),XYZ(6)]);



delete(ax.UserData.singleGreenH);
delete(ax.UserData.singleRedH);
src.Enable='off';


updateNEdge(edgeH);

end

function Select3(src,evn,edgeH,btnSelectEdge1,btnSelectEdge2)
ax=gca;
if  src.Value ==1
    
    btnSelectEdge2.Value=0;
    btnSelectEdge1.Value=0;
    src2.Value=2;
    Select1(src2,[],edgeH,[],[],ax,[])
    Select2(src2,[],edgeH,[],[],ax,[])
    
    
    for k=1:length(edgeH)
        edgeH{k}.HitTest='on';
        edgeH{k}.ButtonDownFcn=@(src,evn)hitedge3(src,evn,edgeH);
    end
    
    for k2=1:length(ax.UserData.extraEdge)
        ax.UserData.extraEdge{k2}.plotH.HitTest='on';
        ax.UserData.extraEdge{k2}.plotH.ButtonDownFcn=@(src,evn)hitedge3(src,evn,edgeH);
    end
    
    title(ax,'Please recover/cancel Edge')
else
    
    for k=1:length(edgeH)
        edgeH{k}.ButtonDownFcn='';
        edgeH{k}.HitTest='off';
    end
    for k2=1:length(ax.UserData.extraEdge)
        ax.UserData.extraEdge{k2}.plotH.HitTest='off';
        ax.UserData.extraEdge{k2}.plotH.ButtonDownFcn='';
    end
    
    
    title(ax,'')
end

end

function hitedge3(src,evn,edgeH)

if evn.Button==1   %left click-> Add
    src.Color=[0.2 ,0.2, 0.2];
    src.LineStyle='-';
elseif evn.Button==3
    src.Color=[0.8 ,0.8, 0.8];
    src.LineStyle=':';
end
updateNEdge(edgeH) ;

end


function Select1(src,evn,edgeH,btnSelectEdge2, editH1,ax,btnTogSelectorDelete)  %src = btnSelectEdge1
if  src.Value ==1
    btnSelectEdge2.Value=0;
    src2.Value=2;
    Select2(src2,[],edgeH,[],[],ax,[]);
    btnTogSelectorDelete.Value=0;
    Select3(src2,[],edgeH,[],[]);
    
    
    for k=1:length(edgeH)
        if ~(sum(edgeH{k}.Color==[0,1,0])==3 || sum(edgeH{k}.Color==[1,0,0])==3 || sum(edgeH{k}.Color==[0.8,0.8,0.8])==3 )
            edgeH{k}.Color=[0.2,0.2,0.2];
            edgeH{k}.LineStyle= '--'; edgeH{k}.UserData.edgeindex=k;
            edgeH{k}.ButtonDownFcn=@(src,evn)hitedge1(src,evn,edgeH,editH1);
            edgeH{k}.HitTest='on';
        end
        
        if sum(edgeH{k}.Color==[1,0,0])==3
            edgeH{k}.ButtonDownFcn=@(src,evn)hitedge1(src,evn,edgeH,editH1);
            edgeH{k}.HitTest='on';
        end
    end
    title(ax,'Please select Edge 1(red)')
else
    for k=1:length(edgeH)
        edgeH{k}.ButtonDownFcn='';
        edgeH{k}.HitTest='off';
    end
    title(ax,'')
end

end

function Select2(src,evn,edgeH,btnSelectEdge1, editH2,ax,btnTogSelectorDelete)  %src = btnSelectEdge2
if  src.Value ==1
    btnSelectEdge1.Value=0;src2.Value=2;
    Select2(src2,[],edgeH,[],[],ax,[]);
    btnTogSelectorDelete.Value=0;
    Select3(src2,[],edgeH,[],[]);
    
    for k=1:length(edgeH)
        if  ~(sum(edgeH{k}.Color==[1,0,0])==3 || sum(edgeH{k}.Color==[0,1,0])==3  || sum(edgeH{k}.Color==[0.8,0.8,0.8])==3)
            edgeH{k}.Color=[0.2,0.2,0.2];
            edgeH{k}.LineStyle= '--'; edgeH{k}.UserData.edgeindex=k;
            edgeH{k}.ButtonDownFcn=@(src,evn)hitedge2(src,evn,edgeH,editH2);
            edgeH{k}.HitTest='on';
        end
        if sum(edgeH{k}.Color==[0,1,0])==3
            edgeH{k}.ButtonDownFcn=@(src,evn)hitedge1(src,evn,edgeH,editH1);
            edgeH{k}.HitTest='on';
        end
    end
    title(ax,'Please select Edge 2(green)')
else
    for k=1:length(edgeH)
        edgeH{k}.ButtonDownFcn='';
        edgeH{k}.HitTest='off';
    end
    title(ax,'')
    
end

end


function hitedge1(src,evn,edgeH,editH1)  %src =one selected edge
selEdge=src.UserData.edgeindex ;
NumberofEqDidPoint= str2double( editH1.String);
for k=1:length(edgeH)
    if k==selEdge
        edgeH{selEdge}.Color=[1,0,0];
        edgeH{selEdge}.LineStyle= '-';
        
        if isfield(edgeH{k}.UserData, 'scatterH')
            for ss=1:length(edgeH{k}.UserData.scatterH)
                delete(edgeH{k}.UserData.scatterH{ss});
            end
        end
        
        xarr=linspace(edgeH{selEdge}.XData(1),edgeH{selEdge}.XData(2),NumberofEqDidPoint+2) ;
        yarr=linspace(edgeH{selEdge}.YData(1),edgeH{selEdge}.YData(2),NumberofEqDidPoint+2) ;
        zarr=linspace(edgeH{selEdge}.ZData(1),edgeH{selEdge}.ZData(2),NumberofEqDidPoint+2) ;
        
        %     edgeH{selEdge}.UserData.scatterH=scatter3(xarr',yarr',zarr','r');
        %     edgeH{selEdge}.UserData.scatterH.ButtonDownFcn=@(src,env)selectScatterRed(src,evn);
        for i=1:length(xarr)
            edgeH{selEdge}.UserData.scatterH{i}=scatter3(xarr(i),yarr(i),zarr(i),'r');
            edgeH{selEdge}.UserData.scatterH{i}.ButtonDownFcn=@(src,env)selectScatterRed(src,evn);
        end
    elseif sum(edgeH{k}.Color==[0,1,0])==3 || sum(edgeH{k}.Color==[0.8,0.8,0.8])==3
        
    else
        edgeH{k}.Color=[0.2,0.2,0.2];
        edgeH{k}.LineStyle= '-';
        if isfield(edgeH{k}.UserData, 'scatterH')
            for ss=1:length(edgeH{k}.UserData.scatterH)
                delete(edgeH{k}.UserData.scatterH{ss});
            end
        end
    end
end
end

function TF1 = contains(Oneline,CriticalWords)   %work for V2015b, check 3/17
TF1=false;
for jj=1:length(CriticalWords)
    k = strfind(Oneline,CriticalWords{jj}) ;
    if ~isempty(k)
        TF1=true;
    end
end
end


function hitedge2(src,evn,edgeH,editH2)  %src =one selected edge
selEdge=src.UserData.edgeindex ;
NumberofEqDidPoint= str2double( editH2.String);

for k=1:length(edgeH)
    if k==selEdge
        edgeH{selEdge}.Color=[0,1,0];
        edgeH{selEdge}.LineStyle= '-';
        
        if isfield(edgeH{k}.UserData, 'scatterH')
            for ss=1:length(edgeH{k}.UserData.scatterH)
                delete(edgeH{k}.UserData.scatterH{ss});
            end
        end
        xarr=linspace(edgeH{selEdge}.XData(1),edgeH{selEdge}.XData(2),NumberofEqDidPoint+2) ;
        yarr=linspace(edgeH{selEdge}.YData(1),edgeH{selEdge}.YData(2),NumberofEqDidPoint+2) ;
        zarr=linspace(edgeH{selEdge}.ZData(1),edgeH{selEdge}.ZData(2),NumberofEqDidPoint+2) ;
        for i=1:length(xarr)
            edgeH{selEdge}.UserData.scatterH{i}=scatter3(xarr(i),yarr(i),zarr(i),'g');
            edgeH{selEdge}.UserData.scatterH{i}.ButtonDownFcn=@(src,env)selectScatterGreen(src,evn);
        end
    elseif sum(edgeH{k}.Color==[1,0,0])==3 || sum(edgeH{k}.Color==[0.8,0.8,0.8])==3
        
    else
        edgeH{k}.Color=[0.2,0.2,0.2];
        edgeH{k}.LineStyle= '-';
        
        if isfield(edgeH{k}.UserData, 'scatterH')
            for ss=1:length(edgeH{k}.UserData.scatterH)
                delete(edgeH{k}.UserData.scatterH{ss});
            end
        end
    end
end
end

function selectScatterRed(src,evn)
ax=gca;
%---finding which node
Target=evn.IntersectionPoint ;
XYZ=[src.XData', src.YData' ,src.ZData' ]- ones(length(src.XData),1)*Target ;
d= sqrt(XYZ(:,1).^2 + XYZ(:,2).^2 + XYZ(:,3).^2);
ind=find(d==min(d));
%--
if isfield(ax.UserData,'singleRedH')
    delete(ax.UserData.singleRedH) ;
end

ax.UserData.singleRedH=scatter3(src.XData(ind),src.YData(ind),src.ZData(ind),'r','filled');
ax.UserData.singleRedH.SizeData=52;

if isfield(ax.UserData,'singleGreenH')
    if isvalid(ax.UserData.singleGreenH)
        ax.UserData.btnAddH.Enable='on';
    end
end

end


function selectScatterGreen(src,evn)
ax=gca;
%---finding which node
Target=evn.IntersectionPoint ;
XYZ=[src.XData', src.YData' ,src.ZData' ]- ones(length(src.XData),1)*Target ;
d= sqrt(XYZ(:,1).^2 + XYZ(:,2).^2 + XYZ(:,3).^2);
ind=find(d==min(d));
%--
if isfield(ax.UserData,'singleGreenH')
    delete(ax.UserData.singleGreenH);
end

ax.UserData.singleGreenH=scatter3(src.XData(ind),src.YData(ind),src.ZData(ind),'g','filled');
ax.UserData.singleGreenH.SizeData=52;

if isfield(ax.UserData,'singleRedH')
    if isvalid(ax.UserData.singleRedH)
        ax.UserData.btnAddH.Enable='on';
    end
end
end

function updateNEdge(edgeH)
ax=gca;
sdfsf=324;
N=0;

%----Original Edge
for k=1:length(edgeH)
    if sum(edgeH{k}.Color==[1,0,0])==3  ||  sum(edgeH{k}.Color==[0,1,0])==3 || sum(edgeH{k}.Color==[0.2,0.2,0.2])==3
        N=N+1;
    end
end
%------Extra Edge
if isfield(ax.UserData,'extraEdge')
    for k2=1:length( ax.UserData.extraEdge)
        cc=ax.UserData.extraEdge{k2}.plotH.Color;
        if sum(cc==[0.2,0.2,0.2])==3
            N=N+1;
        end
    end
end

ax.UserData.textH.String=strcat('Total Bundle =',num2str(N));
end

function num=str2num_CM(str)
try num=str2num(str);
catch num=str;
end
end

