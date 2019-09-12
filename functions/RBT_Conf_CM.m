% Use to load oxDNA topology and configuration(input prova, or last conf.)
% square lattice
% 

% author: Chao-Min Huang, OSU , 03/29/2018


% addpath('jsonlab') ;
% addpath('absor') ;



% %  
[Topology_filename,PathName1,FilterIndex] = uigetfile({'*.dat;*.top','Topogyformat '},'Select the topology file');
Spwd=pwd; 
cd(PathName1);
[Conf_filename,PathName2,FilterIndex]= uigetfile({'*.dat;*.conf','Configuration' },'Select the configuration file');
[JSON_filename,PathName3,FilterIndex]= uigetfile({'*.json','JSONfile' },'Select the json file');
cd(Spwd);

%--------
  dat=loadjson(strcat(PathName3,JSON_filename));
ColorCode=[-1 ,-1 ,-1]; 
EffecCyl=[];
  for k=1:length(dat.vstrands)
      if ~isempty(dat.vstrands{k}.stap_colors)
       ColorCode=[ColorCode;   [k*ones(size(dat.vstrands{k}.stap_colors,1),1), dat.vstrands{k}.stap_colors   ]];
       EffecCyl=union(EffecCyl,k);
      end
  end
 ColorCode=setdiff(ColorCode,[-1,-1,-1],'rows') ;
  [~,b]=hist(ColorCode(:,3),unique(ColorCode(:,3)));   colorRGB=b ;  
  nBundle=length(b);
 %------------ 




delimiterIn = ' ';
A= importdata(strcat(PathName1,Topology_filename),delimiterIn,1) ;

FirstRow=strsplit(A.textdata{1,1});
NBase=str2double(FirstRow{1});
% NStrand=str2double(FirstRow{2});
TSeq= A.textdata((2:NBase+1)',2);

QQ=A.textdata((2:NBase+1)',1);
TStrand=zeros(size(QQ));
for k=1:length(QQ)
TStrand(k)=str2double(QQ{k});
end

B= importdata(strcat(PathName2,Conf_filename),delimiterIn,3) ;
ConfBB=B;
T= B.data;
T0=T ;
 %-----------------
%  figure(11);clf;hold on; axis equal;

pB=T(:,1:3);
Bvec=T(:,4:6);
Axvec=T(:,7:9);

ind=1:size(pB,1);

fH=figure;
clf;
hold on; axis equal;

Ts=unique(TStrand);
Coeff=-0.4;
Save=zeros(NBase,1);
cc=1;
mmMM=[ min(T(:,1:3)); max(T(:,1:3))]' ;
[a,b]=hist(TStrand,unique(TStrand)) ;a0=a  ;b0=b;
scafInitNotFirst=0;

 if 1~= find(a==max(a))
     sacfInd=find(a==max(a));
     QQT=T;
     QQTStrand=TStrand;
     PrevStappleInd = sum(a(1:sacfInd-1))  ;
     ScafIndAll=PrevStappleInd+1:a(sacfInd)+PrevStappleInd  ;
     T(ScafIndAll,:)=[];
     TStrand(ScafIndAll,:)=[];
     T=[QQT(ScafIndAll,:) ;T];
     TStrand=[QQTStrand(ScafIndAll,:) ;TStrand];
     QQQ=TStrand;
      [a1,~,~]=unique(TStrand,'stable') ;
%     IndOrder= unique(TStrand,'stable') ;
     for subs=1:length(a1)
     TStrand(QQQ==a1(subs))=subs;
     end
     scafInitNotFirst=1;
 end
 
[a,b]=hist(TStrand,unique(TStrand)) ;
a0;
StramdIndTable=[a;a0];
  
  
GloIndex=1; Coeff=-0.4;

pHss=cell(1, length(Ts)) ;

for strandi=1: length(Ts)
    included=TStrand==strandi;
    BVechere=T(included,4:6);
    xpp=T(included,1) +Coeff*BVechere(:,1)   ;   %Backbone position
    ypp=T(included,2) +Coeff*BVechere(:,2) ;
    zpp=T(included,3) +Coeff*BVechere(:,3) ;

    PartXYZ=[xpp,ypp,zpp];
%         pHss{strandi}.LineWidth=0.8;
      
      if strandi== find(a==max(a))
        xc = xpp';
        yc = ypp';
        zc = zpp';
        col = (1:length(ypp))*1000;  % This is the color
        pHss{strandi}=surface([xc;xc],[yc;yc],[zc;zc],[col;col],...
        'facecol','no', 'edgecol','interp', 'linew',2);
        %         pHss{strandi}.Color=[0,0,1];
        pHss{strandi}.UserData=GloIndex:length(xpp)  ;
         
      else
        pHss{strandi} =plot3(xpp,ypp,zpp,'.-');     %backbone
        pHss{strandi}.LineWidth=0.1;
        pHss{strandi}.Color=[1,0,0];  
        pHss{strandi}.UserData=GloIndex:GloIndex+length(xpp)-1  ;
        %            pHss{strandi}.HitTest='off'  ;
        pHss{strandi}.MarkerFaceColor=pHss{strandi}.Color;
      end
    
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
% xlim(mmMM(1,:)+[-clearanceA ,clearanceA]) ; ylim(mmMM(2,:)+[-clearanceA ,clearanceA]);zlim(mmMM(3,:)+[-clearanceA ,clearanceA]);

set(gca,'xlim', mmMM(1,:)+[-clearanceA ,clearanceA])
set(gca,'ylim', mmMM(2,:)+[-clearanceA ,clearanceA])
set(gca,'zlim',mmMM(3,:)+[-clearanceA ,clearanceA])

grid on;


%%

%%
drawnow;
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFig = get(handle(fH), 'JavaFrame'); 
jFig.setMaximized(true);
drawnow;
      
aH=gca;  
aH.Position=[0.05,0.05, 0.7, 0.9]; 

Str={'--'};
for k=1:nBundle
    Str{k}=strcat('Bundle  ',num2str(k));    
end
        %---------
popupH = uicontrol('Style', 'popup',...
'String', Str,'Unit','normalized','Position', [0.8 0.87 0.1 0.08]); 
checkH = uicontrol('Style', 'checkbox','String', 'Translation/Rotate','Unit','normalized','Position', [0.8 0.6 0.1 0.05]); 
editH = uicontrol('Style', 'edit','String', '1','Unit','normalized','Position', [0.83 0.85 0.05 0.05]); 
txtH2 = uicontrol('Style','text','Unit','normalized','FontSize',14,'Position', [0.89 0.85 0.05 0.05],....   
'String',' Unit: nm');   
checkH.Callback=@(src,evn)checkFcn(src,evn,editH,txtH2)  ; 

dxy=2.6        ;
NumList=[-1,-1,-1,-1];    %[k, num, col , row, GcX, GcY]
Eff=[];
for k=1:length(dat.vstrands)
   NumList=[NumList;   [k, dat.vstrands{k}.num,  dat.vstrands{k}.col,  dat.vstrands{k}.row ]    ];   
end
NumList=setdiff(NumList,[-1,-1,-1,-1],'rows');
NumList(:,end+1:end+2)=0;
for kk2=1:size(NumList,1)
    NumList(kk2,5)= dxy * NumList(kk2,3) ;
    NumList(kk2,6)= dxy * NumList(kk2,4) ;       
end
ValidAll=NumList ;
  
coeffZ=0.3898 ;
  
  
ConvT_CylBase=zeros(NBase,2); % [Cyy(row1 index) , Base]   %Appendix of T
 
for bbi=1 : NBase
    DD= ValidAll(:,5:6)- ones(size(ValidAll,1) ,1)*T(bbi,1:2)   ;
    d=DD(:,1).^2  +DD(:,2).^2 ;
    ConvT_CylBase(  bbi,1)= ValidAll(  (d==min(d)) ,1)  ;

    ConvT_CylBase(  bbi,2)=   round(T(bbi,3)/coeffZ);
end
 
ValidAll ;
BelongTransM=zeros(NBase,1);
cc=0;
for stapSSi= 1:    size(ColorCode,1)    % fill ds part
    CylInd=ColorCode(stapSSi,1) ;
    BaseStart=ColorCode(stapSSi,2) ;
    ColorIndex= find(colorRGB==ColorCode(stapSSi,3)) ;  % also mean which transformation matrix
    Next= dat.vstrands{CylInd}.stap(BaseStart+1,3:4)  ; 

    stapAndSacf= intersect(find(ConvT_CylBase(:,1)==CylInd),find(ConvT_CylBase(:,2)==BaseStart));
    BelongTransM(stapAndSacf)=ColorIndex;

    cc=cc+length(stapAndSacf);
    nloop=0;
    while sum(Next==[-1,-1])~=2
        CylIndNext= ValidAll(ValidAll(:,2)==Next(1),1) ;
        BaseNext= Next(2);
        stapAndSacf= intersect(find(ConvT_CylBase(:,1)==CylIndNext),find(ConvT_CylBase(:,2)==BaseNext)) ;

        BelongTransM(stapAndSacf)=ColorIndex;
        Next= dat.vstrands{CylIndNext}.stap(BaseNext+1,3:4)  ; 
        nloop=nloop+1;

        cc=cc+length(stapAndSacf);
    end
    nloop;
end
    sum(BelongTransM~=0)

    SSBelongTransM=BelongTransM ;
%-------------  stap
RemainNotAssign=find(BelongTransM==0);  LRNA1=length(RemainNotAssign)

StapNoAss=  RemainNotAssign(RemainNotAssign>sum(TStrand==1)) ;
nn=1;
while ~isempty(StapNoAss)
    for k=1:length(StapNoAss)
        Cand1= StapNoAss(k)+1  ;
        Cand2= StapNoAss(k)-1; 

        ThisStapBaseBC=  ConvT_CylBase( StapNoAss(k),:);
        Cand1BaseBC=  ConvT_CylBase( StapNoAss(k)+1,:);
        Cand2BaseBC=  ConvT_CylBase( StapNoAss(k)-1,:);
        %--------new added June 21
        normD1=norm(Cand1BaseBC-ThisStapBaseBC);
        normD2=norm(Cand2BaseBC-ThisStapBaseBC);
        if normD1<normD2
        BelongTransM(StapNoAss(k))=  BelongTransM(Cand1);
        else
        BelongTransM(StapNoAss(k))=  BelongTransM(Cand2);  
        end
        %------------
    end
    RemainNotAssign=find(BelongTransM==0);
    StapNoAss=  RemainNotAssign(RemainNotAssign>sum(TStrand==1)) ;
    length(StapNoAss);
    nn=nn+1;
    if nn>80
    break
    end

end
  
%--------------- scaffold
ConvT_CylBase;  % allnt, show by [GlobalShowingIndex, Bas] 
CylBaseBynumIndex=zeros(size(ConvT_CylBase));
CylBaseBynumIndex(:,1)=  ValidAll(ConvT_CylBase(:,1) ,2);
CylBaseBynumIndex(:,2)=ConvT_CylBase(:,2) ;

ValidAll;
RemainNotAssign=find(BelongTransM==0);   LRNA2=length(RemainNotAssign)
ScafNoAss=  RemainNotAssign(RemainNotAssign<=sum(TStrand==1)) ;
UnL=length(ScafNoAss) ;
nwh=1;
while ~isempty(ScafNoAss)     
       for kk=1:length(ScafNoAss)
             CylBase=   ConvT_CylBase(ScafNoAss(kk),:) ;
             CADNANOINDEX=  ValidAll( ValidAll(:,1)==CylBase(1),2);
             TwoSideConnec = dat.vstrands{CylBase(1)}.scaf(CylBase(2)+1,:)  ;%read cadnano scaf 

             [~,indLeft]= ismember( TwoSideConnec(1:2),CylBaseBynumIndex,'rows') ;
             if TwoSideConnec(1)~=CADNANOINDEX  || abs(CylBase(2)-TwoSideConnec(2))~=1
                 indLeft=0;
             end
             
             [~,indRight]= ismember( TwoSideConnec(3:4),CylBaseBynumIndex,'rows');
             if TwoSideConnec(3)~=CADNANOINDEX || abs(CylBase(2)-TwoSideConnec(4))~=1
                 indRight=0;
             end
             
             if indLeft==0 
                 if indRight~=0
                 LRHasAss=[ 0  , BelongTransM(indRight)~=0  ];
                 else
                  LRHasAss=[ 0  , 0 ];    
                 end
             elseif indRight==0
             LRHasAss=[ BelongTransM(indLeft)~=0  ,0  ];
             else 
             LRHasAss=[ BelongTransM(indLeft)~=0  , BelongTransM(indRight)~=0  ];                              
             end
            if LRHasAss(1)==1
                Trust=TwoSideConnec(1:2);
            elseif  LRHasAss(2)==1
              Trust=TwoSideConnec(3:4);
            else
               Trust=[];
            end         
            
            if ~isempty(Trust)  && sum(Trust==[-1,-1])~=2
                    CylFirstRep= ValidAll(ValidAll(:,2)==Trust(1),1) ;

                    Gindex= intersect(find(ConvT_CylBase(:,1)==CylFirstRep),find(ConvT_CylBase(:,2)==Trust(2)));
                    Gindex(Gindex>sum(TStrand==1))=[];

                    if  BelongTransM(Gindex)~=0

                      BelongTransM( ScafNoAss(kk))= BelongTransM(Gindex)   ;
                    else
                        nothinghappen=1 ;
                    end
            end
       end
    [nwh length(ScafNoAss) ]
    RemainNotAssign=find(BelongTransM==0);
    ScafNoAss=  RemainNotAssign(RemainNotAssign<=sum(TStrand==1)) ;
    nwh=nwh+1 ;

    if nwh>80
    break 
    end
       
end
%   BelongTransM(ScafNoAss)=3;
  for kkss=1:200
    RemainNotAssign=find(BelongTransM==0);
    ScafNoAss=  find(BelongTransM==0); 
    if isempty(ScafNoAss)
        break;
    end
    
    for SSNi=1:length(ScafNoAss)
        C2Cylinder=CylBaseBynumIndex(ScafNoAss(SSNi),1);
        N1Cyl=CylBaseBynumIndex(ScafNoAss(SSNi)-1,1);
        P1Cyl=CylBaseBynumIndex(ScafNoAss(SSNi)+1,1);
        C2CColRow=NumList(NumList(:,2)==C2Cylinder ,3:4);
        N1CColRow=NumList(NumList(:,2)==N1Cyl ,3:4);
        P1CColRow=NumList(NumList(:,2)==P1Cyl ,3:4);
        S1=C2Cylinder==N1Cyl;
        if norm(C2CColRow-N1CColRow)==1
        S1=true;
        end
        S2=C2Cylinder==P1Cyl;
        if norm(C2CColRow-P1CColRow)==1
        S2=true;
        end
        
        if BelongTransM( ScafNoAss(SSNi)-1)~=0 && S1
        BelongTransM(ScafNoAss(SSNi)) =  BelongTransM(ScafNoAss(SSNi)-1 );
        elseif  BelongTransM( ScafNoAss(SSNi)+1)~=0 && S2
        BelongTransM(ScafNoAss(SSNi)) =  BelongTransM(ScafNoAss(SSNi)+1 ) ;
        end
    end
  end
 
if  nwh>80  %still have unassigned
      for ssss=1:2
           for kk=1:length(ScafNoAss)
                targetXYZ= T(ScafNoAss(kk),1:3) ;

                d= (T(:,1:3)- ones(size(T,1) ,1)*targetXYZ ).^2 ;
                DD=d(:,1).^2+ d(:,2).^2+ d(:,3).^2  ;    
                DD(DD==0)=50000;
                DD( BelongTransM==0)=500000;

                [Sortd,~] = sort(DD);
                closest= Sortd(1:2); 
                closest=closest( randi([1,length(closest)]));
                Choose= find(DD==closest);
               if  BelongTransM(Choose)~=0
                  BelongTransM( ScafNoAss(kk))= BelongTransM(Choose)   ; 

               end
           end
            RemainNotAssign=find(BelongTransM==0);
            ScafNoAss=  RemainNotAssign(RemainNotAssign<=sum(TStrand==1)) ;

            sum(BelongTransM==0);
         if  sum(BelongTransM==0)==0
             break;
         end
      end
end
  
RemainNotAssign=find(BelongTransM==0) ;

[aa,bb]=hist(BelongTransM,unique(BelongTransM)) 

BelongTransM;
AsswhichTM=1 ;
    %% finish searching, shouldn't be any base along
    
for iPh=2: length(pHss)
    BundleInd=  unique(BelongTransM(pHss{iPh}.UserData)) ;
    rgbDec=colorRGB( BundleInd) ;
    RGBHex=dec2hex(rgbDec,6 ) ;
    Rc=hex2dec(RGBHex(1:2));
    Gc=hex2dec(RGBHex(3:4)) ; 
    Bc=hex2dec(RGBHex(5:6))  ;

    pHss{iPh}.Color=[Rc, Gc ,Bc ]/256;   
end
   


set(fH,'WindowKeyPressFcn',{@(src,evn)keyMove(src,evn,fH,pHss,BelongTransM,TStrand,T ,popupH ,checkH ,editH)  }) ;
% set(fH,'WindowKeyReleaseFcn',{@(src,evn)keyMove(src,evn,fH,pHss,BelongTransM,TStrand,T ,popupH ,checkH ,editH)  }) ;

popupH.Callback=@(src,evn) popupFcn(src,evn,pHss, fH ,BelongTransM,TStrand,T ,colorRGB  );
xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
btn1 = uicontrol('Style', 'pushbutton', 'String', 'Export','Unit','normalized', 'Position', [0.8 0.1 0.1 0.1] ,...
'Callback', {@(src,evn)ExportNewConf(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable ,PathName2  )});       

for k=1:length(pHss)
        pHss{k}.ButtonDownFcn=@(src,evn)lineselect(src,evn,popupH, pHss , fH ,BelongTransM,TStrand,T ,colorRGB  );
end
checkH.Position = [ 0.8000    0.8000    0.1000    0.0500 ] ;

sld1 = uicontrol('Style', 'slider','Parent',fH,'Units','normalized',....
'Min',1,'Max',3,'Value',2,'Position', [0.8 0.4 0.15 0.04]);   

sld2 = uicontrol('Style', 'slider','Parent',fH,'Units','normalized',....
'Min',0.005,'Max',5,'Value',1,'Position', [0.8 0.6 0.15 0.04]);    

sld1.Callback= @(src,evn)sldscale1(src,evn,pHss,BelongTransM) ;
sld2.Callback= @(src,evn)sldscale2(src,evn,pHss,BelongTransM) ;

sld1.SliderStep=[0.002,0.05]; sld1.BackgroundColor=[0.6,0.6,1] ;
txtSld1 = uicontrol('Style','text' ,'Units','normalized','Position',[0.8 0.35 0.15 0.05], 'String','Scaffold linewidth');
txtSld2 = uicontrol('Style','text' ,'Units','normalized','Position',[0.8 0.53 0.15 0.05], 'String','Staple linewidth');
sld2.BackgroundColor=[1,0.6,0.6] ;

fileTransInd_name=strcat(PathName2, 'BM.mat') ;
save(fileTransInd_name, 'BelongTransM') ;

txtSld1.HorizontalAlignment='left';  txtSld1.FontSize=16 ;
txtSld2.HorizontalAlignment='left';  txtSld2.FontSize=16 ;
                                       txtH2.FontSize=16 ;
                                      checkH.FontSize=16 ;
                                       editH. FontSize=12 ;
                                       btn1. FontSize=12 ;
                                       popupH.FontSize=16 ;
txtInstruction = uicontrol('Style','text' ,'Units','normalized','Position',[0.8 0.83 0.15 0.15], 'String','Staple linewidth');

chr = 'Rigid Body Transform Code V1';
chr = [chr newline 'mouse/popup to select bundle'];
chr = [chr newline 'Use keyboard to transform'];
chr = [chr newline 'q(+X)   w(+Y)   e(+Z)'];
chr = [chr newline 'a( -X)    s( -Y)    d( -Z)'];

txtInstruction.String =chr; txtInstruction. FontSize=14 ;
txtInstruction.HorizontalAlignment='left';                       
 
txtPop = uicontrol('Style','text' ,'Units','normalized','Position',[popupH.Position(1)+popupH.Position(3)+0.1 0 0.1 0.05]);
chr = 'Select';
chr = [chr newline 'Bundle'];
txtPop.String=chr ; txtPop.FontSize=16 ;  txtPop.HorizontalAlignment='left';
                          
                                       
chr = 'Increment';chr = [chr newline 'Unit : .85nm'];txtH2.String =chr; 

txtH2.Position(2)=0 ; txtH2.Position(3:4)= [0.08,  0.1] ; 
txtH2.HorizontalAlignment='left';


align([sld1 txtSld1 sld2,txtSld2,popupH , editH,checkH,btn1 ],'Left','fixed',5); 
align([sld1 txtSld1 sld2,txtSld2,popupH , editH,checkH,btn1 ,txtInstruction],'Left','distribute'); 

align([txtH2, editH],'fixed',5,'top');
align([txtPop, popupH],'fixed',5,'top');

align([sld1, txtSld1],'Left','fixed',3);
align([sld2, txtSld2],'Left','fixed',3);

% panel = uipanel(fH);
bg = uibuttongroup('Visible','on','Position',[0.79 0 .21 1]); 
uistack(bg,'bottom') ;bg.Title='control panel' ;


set(gca,'Fontsize',18) ;

    return
%      %-------new in V5
     CC=0.6; tol=0.0008 ; 
     BaseCentOnCly= [T0(:,1)+CC*T0(:,4),T0(:,2)+CC*T0(:,5),T0(:,3)+CC*T0(:,6)] ;
    UU_BCOC=uniquetol(BaseCentOnCly,tol,'ByRows',true) ;
    size(UU_BCOC);
    NofDsSsDomain= size(UU_BCOC,1) ;
    ScafStapCorrelate= zeros(NofDsSsDomain,2) ;
     for k=1: size(ScafStapCorrelate,1)
     [~,q2]=ismembertol(BaseCentOnCly,UU_BCOC(k,:),tol,'ByRows',true)  ;
     if sum(q2)==2
         ind=find(q2) ;
     ScafStapCorrelate(k,1)=ind(1) ;
     ScafStapCorrelate(k,2)=ind(2) ;
     end         
     end

    ScafStapCorrelate2=ScafStapCorrelate ;
    ScafStapCorrelate2(ScafStapCorrelate2(:,1)==0,:)=[] ;
    ScafStapCorrelate2=ScafStapCorrelate2 -ones(size(ScafStapCorrelate2)) ;   % index starts at 0, python/Oxdna index
%     
%     
        SfSpCorrATCG=cell(size(ScafStapCorrelate2,1),2) ;
     for k=1: size(ScafStapCorrelate2,1)
        SfSpCorrATCG{k,1}=  TSeq{ScafStapCorrelate2(k,1)+1} ;
        SfSpCorrATCG{k,2}=  TSeq{ScafStapCorrelate2(k,2)+1} ;
     end
     
     RateCorrect=zeros(size(SfSpCorrATCG,1),1) ;
     for k=1: size(ScafStapCorrelate2,1)
            if SfSpCorrATCG{k,1}=='C' && SfSpCorrATCG{k,2}=='G'
                RateCorrect(k)=1 ;
            elseif SfSpCorrATCG{k,1}=='G' && SfSpCorrATCG{k,2}=='C'
                RateCorrect(k)=1 ;
            elseif SfSpCorrATCG{k,1}=='A' && SfSpCorrATCG{k,2}=='T'
                RateCorrect(k)=1 ;
            elseif SfSpCorrATCG{k,1}=='T' && SfSpCorrATCG{k,2}=='A'
                RateCorrect(k)=1 ;
            end
     end
    [size(RateCorrect) ,sum(RateCorrect)]
    NofMutualTrap = size(ScafStapCorrelate2,1) 
    cd(PathName1);

     file3_name='dSRemain.conf'  
     fileID2 = fopen(file3_name,'w');
     
     for iF=1:size(ScafStapCorrelate2,1)
            fprintf(fileID2,'{\n' );
            fprintf(fileID2,'type = mutual_trap\n' );
            fprintf(fileID2,'particle = %u\n' ,ScafStapCorrelate2(iF,1));
            fprintf(fileID2,'ref_particle  = %u\n' ,ScafStapCorrelate2(iF,2));
            fprintf(fileID2,'stiff = %u \n' ,100 );
            fprintf(fileID2,'r0 = 1.2 \n'  );
            fprintf(fileID2,'}\n' );  
            fprintf(fileID2,'{\n' );
            fprintf(fileID2,'type = mutual_trap\n' );
            fprintf(fileID2,'particle = %u\n' ,ScafStapCorrelate2(iF,2));
            fprintf(fileID2,'ref_particle  = %u\n' ,ScafStapCorrelate2(iF,1));
            fprintf(fileID2,'stiff = %u \n' ,100 );
            fprintf(fileID2,'r0 = 1.2 \n' );    
            fprintf(fileID2,'}\n' );  
     end

     fclose(fileID2);
     cd(Spwd);
% %------end of new in V5
    

fprintf('(see if = 7213) the longest strand lengths (scaffold) = %i \n',max(a0)) ;
 
%
%--------------------------
    function sldscale1(src,evn,pHss,BelongTransM)
    %change scaffold color
%     sdf=3;
% %      arrX= linspace(0.01,1, size(pHss{1}.XData,2)) ;
%       arrX= linspace(0,0, size(pHss{1}.XData,2)) ;
%      n=2^(src.Value) ;
%      arrY=arrX.^n;
%      pHss{1}.CData = [arrY; arrY] ;
    pHss{1}.LineWidth=src.Value ;
    
    end

    function sldscale2(src,evn,pHss,~)
    %change staple linewidth
        for ph2=2:length(pHss)            
            pHss{ph2}.LineWidth= src.Value;         
            if src.Value<=0.02
                pHss{ph2}.LineWidth= src.Value;
                pHss{ph2}.Visible='off';
            else
                pHss{ph2}.Visible='on'; 
            end
        end
    end


    function lineselect(src,evn,popupH,pHss , fH ,BelongTransM,TStrand,T,colorRGB )
        xy=evn.IntersectionPoint(1:3);
        XYZAll=[src.XData(1,:) ;src.YData(1,:)  ;   src.ZData(1,:)]' ;   
        dXY= XYZAll-  ones(size(XYZAll,1),1)*xy ;
        ds=dXY(:,1).^2 +dXY(:,2).^2  +dXY(:,3).^2;
        Ind= find(ds==min(ds)) ; Ind=Ind(1);
        BelongTransM(Ind);
        BelongTransM(src.UserData(Ind));
        popupH.Value=  BelongTransM(src.UserData(Ind)) ;
        popupFcn(popupH,[],pHss, fH ,BelongTransM,TStrand,T,colorRGB  )
    end


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

        [regParams,q,w]=absor(A3byN,B3byN);
        w.errmax;
        PosV=    regParams.R*A3byN + regParams.t*ones(1,size(B3byN,2)) ;
        BVecNew=    regParams.R*OldBVec ;
        NVecNew=    regParams.R*OldNVec ;
        NewT( IndexAll,1:9)=   [PosV'-Coeff2*BVecNew',   BVecNew',   NVecNew'];
    end

    if ScafInd~=1  %scaf strand is not the first in topology file------found bug in V3
        MoveAheadind= sum(StramdIndTable(2,1:ScafInd-1));
        NewT=[ NewT(ScafL+1:ScafL+MoveAheadind ,: ) ; NewT(1:ScafL,:); NewT(ScafL+MoveAheadind+1:end,:)];
    end

    mmNewT=min(NewT(:,1:3)) ;
    MMNewT=max(NewT(:,1:3)) ;
    boxsize=floor(2*abs(mmNewT-MMNewT)/50)*50  ; 

    NewT(:,1:3)= NewT(:,1:3)-0.5*(MMNewT+mmNewT) +0.5*[max(boxsize),max(boxsize),max(boxsize)];  %shift center to box center

    
    
    
%-----------------------------------------------------------------------    
    file2_name='prova33.conf'
    fileID = fopen([PathName2 file2_name],'w');
    % fileID = fopen(file2_name,'w');
    fprintf(fileID,'t = 0\n');
    boxsize=floor(3*abs(mmNewT-MMNewT)/50)*50   ;
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
            GindexNeedTomove= find(BelongTransM ==   popupH.Value);
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
                            IndMove= popupH.Value==BelongTransM(pHss{kstp}.UserData) ;
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
                 theta= -RotaionIncrement;
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
                        IndMove= popupH.Value==BelongTransM(pHss{kstp}.UserData) ;
                        
                        OriXYZArray=[pHss{kstp}.XData(IndMove) ;pHss{kstp}.YData(IndMove) ;pHss{kstp}.ZData(IndMove)];
                        NewXYZ=  RMat*(OriXYZArray-Gcenter*ones(1,size(OriXYZArray,2)))+ Gcenter*ones(1,size(OriXYZArray,2));

                        pHss{kstp}.XData(IndMove)=NewXYZ(1,:);
                        pHss{kstp}.YData(IndMove)=NewXYZ(2,:);
                        pHss{kstp}.ZData(IndMove)=NewXYZ(3,:);
                    end
                 end 

                 SacfXYZ=[pHss{1}.XData(1,Inscaf);pHss{1}.YData(1,Inscaf);pHss{1}.ZData(1,Inscaf)] ;
                 NewScaf=RMat*(SacfXYZ-Gcenter*ones(1,size(SacfXYZ,2)))+ Gcenter*ones(1,size(SacfXYZ,2));
                        pHss{1}.XData(:,Inscaf)=[NewScaf(1,:);NewScaf(1,:)] ;
                        pHss{1}.YData(:,Inscaf)=[NewScaf(2,:);NewScaf(2,:)] ;
                        pHss{1}.ZData(:,Inscaf)=[NewScaf(3,:);NewScaf(3,:)] ;

                end
             end
                aH=gca;
              switch evn.Character
               
                  case 'r'
                  aH.XLimMode='manual'; 
                   XXX= aH.XLim +2;    xlim(XXX);
                  case 'f'
                   aH.XLimMode='manual';
                  XXX= aH.XLim -2;     xlim(XXX)
                  case 't'
                  aH.YLimMode='manual';
                   YYY= aH.YLim +2;    %ylim(YYY);
                    set(gca,'ylim',YYY)
                  case 'g'
                  aH.YLimMode='manual';
                  YYY= aH.YLim -2;   
%                   ylim(YYY);
                  set(gca,'ylim',YYY)
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
   chr = 'Increment';
            if  src.Value==0
                chr = [chr newline 'Unit : .85nm'] ;
               editH.String = '5' ;  %default translation
            else
                chr = [chr newline 'Unit : deg'] ;
               editH.String = '30' ;  %default rotation
            end
  txtH2.String =chr;         
   end
  
    function popupFcn(src,evn,pHss, fH ,BelongTransM,TStrand,T,colorRGB  )
         GindexNeedTomove= find(BelongTransM ==   src.Value);
                 for kstp=2: length(pHss)
                    BundleInd=  unique(BelongTransM(pHss{kstp}.UserData)) ;
                    rgbDec=colorRGB( BundleInd) ;
                    RGBHex=dec2hex(rgbDec,6 ) ;
                    Rc=hex2dec(RGBHex(1:2));
                    Gc=hex2dec(RGBHex(3:4)) ; 
                    Bc=hex2dec(RGBHex(5:6))  ;
                    pHss{kstp}.Color=[Rc, Gc ,Bc ]/256;           
                    
                    if ~isempty(intersect( pHss{kstp}.UserData,GindexNeedTomove))

                        pHss{kstp}.LineStyle=':' ;
                       pHss{kstp}.Color=[ 1,0,0] ;

                    else
                         pHss{kstp}.LineStyle='-' ; 
                    end
                 end         
    end
  

    
    
  