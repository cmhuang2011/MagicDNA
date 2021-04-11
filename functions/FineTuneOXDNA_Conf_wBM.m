function FineTuneOXDNA_Conf_wBM

PathName1=pwd;
Topology_filename='prova.top' ;
Conf_filename='prova3322.conf';

% [Topology_filename,PathName1,FilterIndex] = uigetfile({'*.dat;*.top','Topogyformat '},'Select the topology file');
Spwd=pwd;
cd(PathName1);
% [Conf_filename,PathName2,FilterIndex]= uigetfile({'*.dat;*.conf','Configuration' },'Select the configuration file');
% [JSON_filename,PathName3,FilterIndex]= uigetfile({'*.json','JSONfile' },'Select the json file');
cd(Spwd);

PathName2=PathName1;
delimiterIn = ' ';
A= importdata(strcat(PathName1,filesep,Topology_filename),delimiterIn,1) ;

FirstRow=strsplit(A.textdata{1,1});
NBase=str2double(FirstRow{1});
% NStrand=str2double(FirstRow{2});
TSeq= A.textdata((2:NBase+1)',2);

QQ=A.textdata((2:NBase+1)',1);
TStrand=zeros(size(QQ));
for k=1:length(QQ)
    TStrand(k)=str2double(QQ{k});
end

B= importdata(strcat(PathName2,filesep,Conf_filename),delimiterIn,3) ;
ConfBB=B;
T= B.data;
%-----------------
%  figure(11);clf;hold on; axis equal;
size(T) ;
pB=T(:,1:3);
Bvec=T(:,4:6);
Axvec=T(:,7:9);

% colorArr=zeros(size(pB,1),3);
% for i=1:size(colorArr,1)
%     if strcmp( TSeq{i},'A')
%         colorArr(i,:)=[1,0,0];
%     elseif strcmp( TSeq{i},'T')
%         colorArr(i,:)=[0,1,0];
%     elseif strcmp( TSeq{i},'C')
%         colorArr(i,:)=[0,0,1];
%     elseif strcmp( TSeq{i},'G')
%         colorArr(i,:)=[0.6,0.6,0];
%     end
% end
% ind=1:size(pB,1);

fH=figure;
clf;
hold on; axis equal;

Ts=unique(TStrand);
Coeff=-0.4;

% Coeff=1.35; line model


Save=zeros(NBase,1);
cc=1;
% mmMM=zeros(3,2);
mmMM=[ min(T(:,1:3)); max(T(:,1:3))]' ;
pHss=cell(1, length(Ts)) ;


[a,b]=hist(TStrand,unique(TStrand)) ;a0=a  ;b0=b;
scafInitNotFirst=0;


% if 1~= find(a==max(a))
%     sacfInd=find(a==max(a));
%     QQT=T;
%     QQTStrand=TStrand;
%     PrevStappleInd = sum(a(1:sacfInd-1))  ;
%     ScafIndAll=PrevStappleInd+1:a(sacfInd)+PrevStappleInd  ;
%     T(ScafIndAll,:)=[];
%     TStrand(ScafIndAll,:)=[];
%     T=[QQT(ScafIndAll,:) ;T];
%     TStrand=[QQTStrand(ScafIndAll,:) ;TStrand];
%     QQQ=TStrand;
%     [a1,~,~]=unique(TStrand,'stable') ;
%     %     IndOrder= unique(TStrand,'stable') ;
%     for subs=1:length(a1)
%         TStrand(QQQ==a1(subs))=subs;
%     end
%     scafInitNotFirst=1;
%     
% end
[a,b]=hist(TStrand,unique(TStrand)) ;
StramdIndTable=[a;a0];


GloIndex=1;


for strandi=1: length(Ts)
    included=TStrand==strandi;
    
    BVechere=T(included,4:6);
    
    xpp=T(included,1) +Coeff*BVechere(:,1)   ;   %Backbone position
    ypp=T(included,2) +Coeff*BVechere(:,2) ;
    zpp=T(included,3) +Coeff*BVechere(:,3) ;
    
    PartXYZ=[xpp,ypp,zpp];
%     
%     if strandi== find(a==max(a))
% %         xc = xpp';
% %         yc = ypp';
% %         zc = zpp';
% %         col = (1:length(ypp))*1000;  % This is the color
% %         pHss{strandi}=surface([xc;xc],[yc;yc],[zc;zc],[col;col],...
% %             'facecol','no', 'edgecol','interp', 'linew',2);
%         %         pHss{strandi}.Color=[0,0,1];
%         pHss{strandi} =plot3(xpp,ypp,zpp,'.-');     %backbone
%         pHss{strandi}.UserData=GloIndex:length(xpp)  ;
%         pHss{strandi}.LineWidth=0.1;
%     else
        pHss{strandi} =plot3(xpp,ypp,zpp,'.-');     %backbone
        pHss{strandi}.LineWidth=0.1;
        
        
        pHss{strandi}.Color=[1,0,0];
        pHss{strandi}.UserData=GloIndex:GloIndex+length(xpp)-1  ;
        %            pHss{strandi}.HitTest='off'  ;
        pHss{strandi}.MarkerFaceColor=pHss{strandi}.Color;
%     end
    
    GloIndex=GloIndex+length(xpp);
    
    for k=1:size(PartXYZ,1)-1
        BackPA=PartXYZ(k,:);
        BackPB=PartXYZ(k+1,:) ;
        d=norm(BackPA-BackPB) ;
        Save(cc)=norm(BackPA-BackPB) ;
        cc=cc+1;
    end
    
    
end
clearanceA=3;
xlim(mmMM(1,:)+[-clearanceA ,clearanceA]) ; ylim(mmMM(2,:)+[-clearanceA ,clearanceA]);zlim(mmMM(3,:)+[-clearanceA ,clearanceA]);
grid on;

%%--------------------------------------------
%%

%%
%
cd(PathName1);
BelongTransM=[];
try
    load(strcat(PathName1 ,filesep, 'BM3.mat'))  ;  % read BelongTransM from previous export
    BelongTransM=BM3 ;
    
catch
    load(strcat(PathName1, filesep ,'BM.mat'))  ;  % read BelongTransM from previous export
    %  BelongTransM=BM;
end
% hard
CylBaseBynumIndex=[];
cd(Spwd);

% RBT_Conf
%   [~,b]=hist(ColorCode(:,3),unique(ColorCode(:,3)));
%
nBundle=max(BelongTransM);
%
aH=gca;
aH.Position=[0.05,0.05, 0.7, 0.9];
%
Str={'--'};
for k=1:nBundle
    Str{k}=strcat('Bundle  ',num2str(k));
end
Str{end+1}='All';
%         %---------
popupH = uicontrol('Style', 'popup',...
    'String', Str,'Unit','normalized','Position', [0.8 0.87 0.1 0.08]);
checkH = uicontrol('Style', 'checkbox','String', 'Trans/Rotate','Unit','normalized','Position', [0.8 0.65 0.1 0.05]);
editH = uicontrol('Style', 'edit','String', '1','Unit','normalized','Position', [0.83 0.85 0.05 0.05]);
txtH2 = uicontrol('Style','text','Unit','normalized','FontSize',14,'Position', [0.89 0.85 0.05 0.05],....
    'String','Unit: nm');
checkH.Callback=@(src,evn)checkFcn(src,evn,editH,txtH2)  ;





set(fH,'KeyPressFcn',{@(src,evn)keyMove(src,evn,fH,pHss,BelongTransM,TStrand,T ,popupH ,checkH ,editH)  }) ;

popupH.Callback=@(src,evn) popupFcn(src,evn,pHss, fH ,BelongTransM,TStrand,T   );

xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
btn1 = uicontrol('Style', 'pushbutton', 'String', 'Export','Unit','normalized', 'Position', [0.8 0.1 0.1 0.1] ,...
    'Callback', {@(src,evn)ExportNewConf(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable ,PathName2  )});

for k=1:length(pHss)
    pHss{k}.ButtonDownFcn=@(src,evn)lineselect(src,evn,popupH, pHss , fH ,BelongTransM,TStrand,T   );
end
checkH_View = uicontrol('Style', 'checkbox','String', 'Axis auto/equal','Unit','normalized','Position', [0.8 0.55 0.05 0.05]);
checkH_View.Callback=@(src,evn)checkFcn2(src,evn,aH)  ;

sld1 = uicontrol('Style', 'slider','Parent',fH,'Units','normalized',....
    'Min',-3,'Max',3,'Value',0,'Position', [0.8 0.72 0.15 0.02]);
sld2 = uicontrol('Style', 'slider','Parent',fH,'Units','normalized',....
    'Min',0.005,'Max',5,'Value',1,'Position', [0.8 0.75 0.15 0.02]);

sld1.Callback= @(src,evn)sldscale1(src,evn,pHss,BelongTransM) ;

sld2.Callback= @(src,evn)sldscale2(src,evn,pHss,BelongTransM) ;


sld1.SliderStep=[0.002,0.05];
end

%
%--------------------------
function sldscale1(src,~,pHss,~)
%change scaffold color
arrX= linspace(0.01,1, size(pHss{1}.XData,2)) ;
n=2^(src.Value) ;
arrY=arrX.^n;
pHss{1}.CData = [arrY; arrY] ;
if src.Value<-2.5
    
    pHss{1}.CData=ones(size(  pHss{1}.CData)) ;
    pHss{1}.FaceColor=[0,0,1];
end


end

function sldscale2(src,~,pHss,~)
%change staple linewidth
%     sdff=4;
pHss{1}.Marker='.';
%     pHss{1}.Marker='none';


for ph2=1:length(pHss)
    pHss{ph2}.MarkerSize=8;
    pHss{ph2}.LineWidth= src.Value;
    if src.Value<=0.02
        pHss{ph2}.LineWidth= src.Value;
        pHss{ph2}.Visible='off';
    else
        pHss{ph2}.Visible='on';
    end
    
end

end


function lineselect(src,evn,popupH,pHss , fH ,BelongTransM,TStrand,T )

xy=evn.IntersectionPoint(1:3);
XYZAll=[src.XData(1,:) ;src.YData(1,:)  ;   src.ZData(1,:)]' ;
%              XYZAll=[pHss{1}.XData( ;pHss{1}.YData  ;   pHss{1}.ZData]' ;

dXY= XYZAll-  ones(size(XYZAll,1),1)*xy ;
ds=dXY(:,1).^2 +dXY(:,2).^2  +dXY(:,3).^2;
Ind= find(ds==min(ds)) ; Ind=Ind(1);
BelongTransM(Ind);
BelongTransM(src.UserData(Ind));
popupH.Value=  BelongTransM(src.UserData(Ind)) ;
popupFcn(popupH,[],pHss, fH ,BelongTransM,TStrand,T  )
end

%
% function  T=convertEulerToRotation(Vars_6)
%
% RMatx=[1 0 0; 0 cos(Vars_6(1)) sin(Vars_6(1)); 0 -sin(Vars_6(1)) cos(Vars_6(1))];
% RMaty=[cos(Vars_6(2)) 0 -sin(Vars_6(2)); 0 1 0; sin(Vars_6(2)) 0 cos(Vars_6(2))];
% RMatz=[cos(Vars_6(3)) sin(Vars_6(3)) 0; -sin(Vars_6(3))  cos(Vars_6(3))  0; 0 0 1];
% M=RMatx*RMaty*RMatz;
% T=[M, [Vars_6(4);Vars_6(5);Vars_6(6)] ; 0 0 0 1];
% end

%---------------------------------------------------------------------------------------------------------

%------------
function ExportNewConf(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand ,StramdIndTable ,PathName2 )
Coeff2=-0.4;


% gathering data into matrix(n by 3) form
BackBone=zeros(size(T,1),3) ;

OldBackBone=[T(:,1)+Coeff2*T(:,4) ,  T(:,2)+Coeff2*T(:,5) ,T(:,3)+Coeff2*T(:,6) ];
NewT=zeros(size(T));

ScafInd=find(StramdIndTable(2,:)==max(StramdIndTable(1,:)));
ScafL=max(StramdIndTable(1,:));

for strandi2=1:length(pHss)
    BackBone( pHss{strandi2}.UserData ,1:3)=[ pHss{strandi2}.XData(1,:)' ,  pHss{strandi2}.YData(1,:)' , pHss{strandi2}.ZData(1,:)' ];  %Backbone position
end



for iTransf= 1:max(BelongTransM)
    IndexAll=BelongTransM==iTransf;
    A3byN=  transpose(OldBackBone(IndexAll,1:3)) ;
    B3byN=  transpose(BackBone(IndexAll,1:3) );
    OldBVec= transpose(T(IndexAll,4:6)) ;
    OldNVec= transpose(T(IndexAll,7:9)) ;
    
    
    %                    OP3=PlotXYZV2Int;
    [regParams,~,w]=absor(A3byN,B3byN);
%     w.errmax
    PosV=    regParams.R*A3byN + regParams.t*ones(1,size(B3byN,2)) ;
    BVecNew=    regParams.R*OldBVec ;
    NVecNew=    regParams.R*OldNVec ;
    
    NewT( IndexAll,1:9)=   [PosV'-Coeff2* BVecNew',   BVecNew',   NVecNew'];
end

% if ScafInd~=1  %scaf strand is not the first in topology file------found bug in V3
%     MoveAheadind= sum(StramdIndTable(2,1:ScafInd-1));
%     
%     NewT=[ NewT(ScafL+1:ScafL+MoveAheadind ,: ) ; NewT(1:ScafL,:); NewT(ScafL+MoveAheadind+1:end,:)];
% end


mmNewT=min(NewT(:,1:3)) ;
MMNewT=max(NewT(:,1:3)) ;

% boxsize= MMNewT-mmNewT +[30,30,30];
boxsize=1.8*ceil(abs(mmNewT-MMNewT)/50)*50   ;

NewT(:,1:3)= NewT(:,1:3) - ones(size(NewT,1),1)*( 0.5*(mmNewT+MMNewT)) ;

fprintf('Exporting prova33.conf \n') ;
file2_name='prova33.conf' ;
fileID = fopen([PathName2 filesep file2_name],'w');
% fileID = fopen(file2_name,'w');
fprintf(fileID,'t = 0\n');
E0=0;
fprintf(fileID,'b = %3.0f %3.0f %3.0f\n',max(boxsize),max(boxsize),max(boxsize));
fprintf(fileID,'E = %8.6f %8.6f %8.6f\n',E0,E0,E0);
for printout=1: size(NewT,1)
    PP=NewT(printout,1:3);
    BB=NewT(printout,4:6);
    NN=NewT(printout,7:9);
    
    fprintf(fileID,'%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f\n',[PP,BB,NN ,zeros(1,3),zeros(1,3)  ] );
    
end
fclose(fileID);

end


function keyMove(src,evn,fH,pHss,BelongTransM,TStrand,T ,popupH ,checkH ,editH)
UseAll=0;
GindexNeedTomove= find(BelongTransM ==   popupH.Value);
if isempty(GindexNeedTomove)
    GindexNeedTomove=1:size(BelongTransM,1) ; UseAll=1;
end
%             sdfsf=3
%             GindexNeedTomove= find(BelongTransM~=0) ;   % special for all bases;

if  checkH.Value==0   %translation
    TransVec=[];
    LinearIncrement= str2double(editH.String);
    if isnan(LinearIncrement)
        return
    end
    
    switch evn.Character
        case 'a'
            TransVec= -[LinearIncrement,0,0] ;
        case 'q'
            TransVec=[LinearIncrement,0,0]   ;
        case 's'
            TransVec=[0,-LinearIncrement,0]      ;
        case 'w'
            TransVec=[0,LinearIncrement,0,0]    ;
        case 'd'
            TransVec=[0,0,-LinearIncrement]  ;
        case 'e'
            TransVec=[0,0,+LinearIncrement];
    end
    if ~isempty(TransVec)
        
        for kstp=2: length(pHss)
            if ~isempty(intersect( pHss{kstp}.UserData,GindexNeedTomove))
                if UseAll==0
                    IndMove= popupH.Value==BelongTransM(pHss{kstp}.UserData) ;
                else
                    IndMove= 1:length(pHss{kstp}.XData );
                end
                
                %                              IndMove=find(ones(size(IndMove))) ; % special for allw
                pHss{kstp}.XData(IndMove)= pHss{kstp}.XData(IndMove)+ TransVec(1);
                pHss{kstp}.YData(IndMove)= pHss{kstp}.YData(IndMove)+ TransVec(2);
                pHss{kstp}.ZData(IndMove)= pHss{kstp}.ZData(IndMove)+ TransVec(3);
            end
        end
        %----scaf
        Inscaf=intersect(pHss{1}.UserData,GindexNeedTomove) ;
        pHss{1}.XData(:,Inscaf)= pHss{1}.XData(:,Inscaf)+ TransVec(1);
        pHss{1}.YData(:,Inscaf)= pHss{1}.YData(:,Inscaf)+ TransVec(2);
        pHss{1}.ZData(:,Inscaf)= pHss{1}.ZData(:,Inscaf)+ TransVec(3);
        
        
    end
else   %rotation
    RotaionIncrement= str2double(editH.String);
    if isnan(RotaionIncrement)
        return
    end
    
    
    RotaionIncrement=RotaionIncrement*pi/180;
    if isempty(RotaionIncrement)
        RotaionIncrement=0.05;
    end
    theta=RotaionIncrement;
    RMat=[];
    switch evn.Character
        case 'a'
            RMat=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
        case 'q'
            RMat=[1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        case 's'
            RMat=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
        case 'w'
            RMat=[cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
        case 'd'
            RMat=[cos(theta) -sin(theta) 0; sin(theta)  cos(theta)  0; 0 0 1];
        case 'e'
            RMat=[cos(theta) sin(theta) 0; -sin(theta)  cos(theta)  0; 0 0 1];
    end
    if ~isempty(RMat)
        
        Inscaf=intersect(pHss{1}.UserData,GindexNeedTomove) ;
        if ~isempty(Inscaf)
            Gcenter=[mean(pHss{1}.XData(1,Inscaf) );mean(pHss{1}.YData(1,Inscaf)) ;mean(pHss{1}.ZData(1,Inscaf)) ];
        else
            G_XYZ=[-1,-1,-1];
            for kstp=2: length(pHss)
                SelectStp= popupH.Value==BelongTransM(pHss{kstp}.UserData) ;
                if sum(SelectStp)>0
                    G_XYZ=union(G_XYZ,[ pHss{kstp}.XData(SelectStp)', pHss{kstp}.YData(SelectStp)' pHss{kstp}.ZData(SelectStp)'] ,'rows') ;
                end
                
            end
            G_XYZ=setdiff(G_XYZ,[-1,-1,-1],'rows');
            Gcenter=mean(G_XYZ)';
        end
        for kstp=2: length(pHss)
            if ~isempty(intersect( pHss{kstp}.UserData,GindexNeedTomove))
                if UseAll==0
                    IndMove= popupH.Value==BelongTransM(pHss{kstp}.UserData) ;
                else
                    IndMove=1:length( pHss{kstp}.XData) ;
                end
                %                          IndMove=find(ones(size(IndMove))) ; % special for allw
                
                
                OriXYZArray=[pHss{kstp}.XData(IndMove) ;pHss{kstp}.YData(IndMove) ;pHss{kstp}.ZData(IndMove)];
                NewXYZ=  RMat*(OriXYZArray-Gcenter*ones(1,size(OriXYZArray,2)))+ Gcenter*ones(1,size(OriXYZArray,2));
                
                pHss{kstp}.XData(IndMove)=NewXYZ(1,:);
                pHss{kstp}.YData(IndMove)=NewXYZ(2,:);
                pHss{kstp}.ZData(IndMove)=NewXYZ(3,:);
            end
        end
        
        SacfXYZ=[pHss{1}.XData(1,Inscaf);pHss{1}.YData(1,Inscaf);pHss{1}.ZData(1,Inscaf)] ;
        NewScaf=RMat*(SacfXYZ-Gcenter*ones(1,size(SacfXYZ,2)))+ Gcenter*ones(1,size(SacfXYZ,2));
        pHss{1}.XData(:,Inscaf)=[NewScaf(1,:)] ;
        pHss{1}.YData(:,Inscaf)=[NewScaf(2,:)] ;
        pHss{1}.ZData(:,Inscaf)=[NewScaf(3,:)] ;
        
    end
end
switch evn.Character
    case 'r'
        aH.XLimMode='manual';
        XXX= aH.XLim +2;    xlim(XXX);
    case 'f'
        aH.XLimMode='manual';
        XXX= aH.XLim -2;     xlim(XXX)
    case 't'
        aH.YLimMode='manual';
        YYY= aH.YLim +2;    ylim(YYY);
    case 'g'
        aH.YLimMode='manual';
        YYY= aH.YLim -2;    ylim(YYY);
    case 'y'
        aH.ZLimMode='manual';
        ZZZ= aH.ZLim +2;    zlim(ZZZ);
    case 'h'
        aH.ZLimMode='manual';
        ZZZ= aH.ZLim -2;    zlim(ZZZ);
    otherwise
        axis equal; xlim auto; ylim auto; zlim auto ;
end

end


function checkFcn(src,~,editH,txtH2)
if  src.Value==0
    txtH2.String = 'Unit : nm';
    editH.String = '5' ;
else
    txtH2.String = 'Unit : deg';
    editH.String = '30' ;
end
end

function popupFcn(src,evn,pHss, fH ,BelongTransM,TStrand,T  )
GindexNeedTomove= find(BelongTransM ==   src.Value);
for kstp=2: length(pHss)
    if ~isempty(intersect( pHss{kstp}.UserData,GindexNeedTomove))
        pHss{kstp}.MarkerSize=6;
        pHss{kstp}.Color=[ 0,0.5,0.6] ;
    else
        pHss{kstp}.MarkerSize=10;
        pHss{kstp}.Color=[ 1,0,0] ;
    end
end
end




