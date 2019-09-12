% Use to load oxDNA topology and configuration(input prova, or last conf.)
% honeycomb lattice 11282017
%  ArrangeJSON_V5_Euler

% addpath(genpath(pwd));
% new in V5: 1. optional export force file which uses mutual trap between a scaf
% base and a stap base.
%   2. cancel function to script

clear;  %close all;

% %  
[Topology_filename,PathName1,FilterIndex] = uigetfile({'*.dat;*.top','Topogyformat '},'Select the topology file');
Spwd=pwd; 
cd(PathName1);
[Conf_filename,PathName2,FilterIndex]= uigetfile({'*.dat;*.conf','Configuration' },'Select the configuration file');

[JSON_filename,PathName3,FilterIndex]= uigetfile({'*.json','JSONfile' },'Select the json file');
cd(Spwd);



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

colorArr=zeros(size(pB,1),3);
for i=1:size(colorArr,1)
    if strcmp( TSeq{i},'A')
        colorArr(i,:)=[1,0,0];
    elseif strcmp( TSeq{i},'T')
       colorArr(i,:)=[0,1,0];
    elseif strcmp( TSeq{i},'C')
       colorArr(i,:)=[0,0,1];
    elseif strcmp( TSeq{i},'G')
       colorArr(i,:)=[0.6,0.6,0];
    end
end
ind=1:size(pB,1);

fH=figure;
clf;
hold on; axis equal;

Ts=unique(TStrand);
Coeff=-0.4;
Save=zeros(NBase,1);
cc=1;
% mmMM=zeros(3,2);
mmMM=[ min(T(:,1:3)); max(T(:,1:3))]' ;
pHss=cell(1, length(Ts)) ;


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
  
  
GloIndex=1;
  
  
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

% Ascat=scatter3(xpp(1),ypp(1),zpp(1),'CData',[1,0,0]);Ascat.SizeData=1;
% Tscat=scatter3(xpp(1),ypp(1),zpp(1),'CData',[0,1,0]);Tscat.SizeData=1;
% Cscat=scatter3(xpp(1),ypp(1),zpp(1),'CData',[0,0,1]);Cscat.SizeData=1;
% Gscat=scatter3(xpp(1),ypp(1),zpp(1),'CData',[0.6,0.6,0]);Gscat.SizeData=1;
% legend([Ascat Tscat Cscat Gscat ],'A','T','C','G')

% for kdd=2:size(T,1)
%     p0=T(kdd,1:3);
%     BVec=T(kdd,4:6);
%     p1=p0 +Coeff*BVec ;
%     
%     pH=plot3([p0(1),p1(1)],[p0(2),p1(2)],[p0(3),p1(3)],'k');
%     sdf=324;
%     
%     if TStrand(kdd)==1
%         pH.Color=[0,0,0];
%     else
%        pH.Color=[0.5,0,0];
% %        pH.Visible='off';
%     end
% 
% end
% Save(Save==0)=[];
% figure(51)
% % figure
% scatter(1:length(Save),Save')
%%--------------------------------------------
%%

%%


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
  
  [~,b]=hist(ColorCode(:,3),unique(ColorCode(:,3)));
  
  nBundle=length(b);
  
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
        checkH = uicontrol('Style', 'checkbox','String', 'Trans/Rotate','Unit','normalized','Position', [0.8 0.65 0.1 0.05]); 
        editH = uicontrol('Style', 'edit','String', '1','Unit','normalized','Position', [0.83 0.85 0.05 0.05]); 
        txtH2 = uicontrol('Style','text','Unit','normalized','FontSize',14,'Position', [0.89 0.85 0.05 0.05],....   
        'String','Unit: nm');   
        checkH.Callback=@(src,evn)checkFcn(src,evn,editH,txtH2)  ; 

  NumList=[-1,-1,-1,-1];    %[k, num, col , row, GcX, GcY]
  Eff=[];
   for k=1:length(dat.vstrands)
   NumList=[NumList;   [k, dat.vstrands{k}.num,  dat.vstrands{k}.col,  dat.vstrands{k}.row ]    ];   
  end
 NumList=setdiff(NumList,[-1,-1,-1,-1],'rows');
 NumList(:,end+1:end+2)=0;
 
 
 dxy=1.277      ;

 %-------------------HC
  for kk2=1:size(NumList,1)
      if mod(NumList(kk2,4)+NumList(kk2,3),2 )==0
      NumList(kk2,6)= dxy * NumList(kk2,4)*3 ;

      else
        NumList(kk2,6)= dxy * NumList(kk2,4)*3+ dxy;
                   
      end
      NumList(kk2,5)= dxy * NumList(kk2,3)*sqrt(3) ;
%       NumList(kk2,5)= 100*dxy * NumList(kk2,3) ;
% %       NumList(kk2,6)= -3.84 * NumList(kk2,4)+65.79 ;       
%        NumList(kk2,6)= 1100*dxy * NumList(kk2,4) ;       
  end
    ValidAll=NumList ;
  
    
%-----------------    
    
  coeffZ=0.3898 ;
  
%   scatter(ValidAll(:,5),ValidAll(:,6),'ok' ) ;
  
  grid minor
  
 ConvT_CylBase=zeros(NBase,2); % [Cyy(row1 index) , Base]   %Appendix of T
 
 for bbi=1 : NBase
   DD= ValidAll(:,5:6)- ones(size(ValidAll,1) ,1)*T(bbi,1:2)   ;
   d=DD(:,1).^2  +DD(:,2).^2 ;
   ConvT_CylBase(  bbi,1)= ValidAll(  (d==min(d)) ,1)  ;
     
%    if isempty(find(d==min(d)))
%        sdfsf=234652;
%    end
   
   
    ConvT_CylBase(  bbi,2)=   round(T(bbi,3)/coeffZ);
 end
 
 ValidAll ;
BelongTransM=zeros(NBase,1);
cc=0;
  for stapSSi= 1:    size(ColorCode,1)    % fill ds part
      CylInd=ColorCode(stapSSi,1) ;
      BaseStart=ColorCode(stapSSi,2) ;
      ColorIndex= find(b==ColorCode(stapSSi,3)) ;  % also mean which transformation matrix
       
       Next= dat.vstrands{CylInd}.stap(BaseStart+1,3:4)  ; 
    
      
     stapAndSacf= intersect(find(ConvT_CylBase(:,1)==CylInd),find(ConvT_CylBase(:,2)==BaseStart));
     BelongTransM(stapAndSacf)=ColorIndex;
     
    if  length(stapAndSacf)~=2
        sdf55=234;
    end
     cc=cc+length(stapAndSacf);
     nloop=0;
     while sum(Next==[-1,-1])~=2
     
        CylIndNext= ValidAll(ValidAll(:,2)==Next(1),1) ;
        BaseNext= Next(2);
       if isempty( CylIndNext)
           sdfsf=45;
       end
%         ConvT_CylBase
        stapAndSacf= intersect(find(ConvT_CylBase(:,1)==CylIndNext),find(ConvT_CylBase(:,2)==BaseNext)) ;

%         if  sum(BelongTransM(stapAndSacf)~=0)>1
%             sdff=234;
%         end
         BelongTransM(stapAndSacf)=ColorIndex;
         Next= dat.vstrands{CylIndNext}.stap(BaseNext+1,3:4)  ; 
         nloop=nloop+1;
         
         cc=cc+length(stapAndSacf);
           if  length(stapAndSacf)~=2
          end
     end
     nloop;
  end

    
     sum(BelongTransM~=0)
    %-------  stap
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
%         rdp=rand() ;
%        if    TStrand(Cand1)== TStrand(StapNoAss(k)) && rdp>0.2
%            BelongTransM(StapNoAss(k))=  BelongTransM(Cand1);
%        elseif  TStrand(Cand2)== TStrand(StapNoAss(k))
%             BelongTransM(StapNoAss(k))=  BelongTransM(Cand2);
%        else
% %            sdfsfsdf=23;
%        end
       
       %-----------
       
%          BelongTransM(StapNoAss(k));
    end
        RemainNotAssign=find(BelongTransM==0);

     StapNoAss=  RemainNotAssign(RemainNotAssign>sum(TStrand==1)) ;
     length(StapNoAss);
      nn=nn+1;
      if nn>80
          break
      end
      
  end
  
%--------------Line saperate

% xLineCut=20;
% ForceAssInd= and(BelongTransM==2, T(:,1)<xLineCut);
% 
%  BelongTransM(  ForceAssInd)=1;

  
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
      if nwh==150
          sff=2;
      end
      
      
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
             
%              LRHasAss=zeros(1,2);
%              LRHasAss=[ BelongTransM(indLeft)~=0  , BelongTransM(indRight)~=0  ];
             
             
%              rdp=rand() ; %rdp=1;
            if LRHasAss(1)==1
                Trust=TwoSideConnec(1:2);
            elseif  LRHasAss(2)==1
              Trust=TwoSideConnec(3:4);
            else
               Trust=[];
            end

%------------------------
%           if nwh>80
%             rdp=rand() ;
%             if abs(CylBase(2)-TwoSideConnec(2))==1  && rdp>0.5
%                 Trust=TwoSideConnec(1:2);
%             elseif  abs(CylBase(2)-TwoSideConnec(4))==1
%               Trust=TwoSideConnec(3:4);
%             else
%                Trust=[];
%             end    
%           end
% ------------------            
            
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
  for kk=1:200
       RemainNotAssign=find(BelongTransM==0);
       ScafNoAss=  RemainNotAssign(RemainNotAssign<=sum(TStrand==1)) ;      
       if isempty(ScafNoAss)
           break;
       end
       
       
       for SSNi=1:length(ScafNoAss)
           C2Cylinder=CylBaseBynumIndex(ScafNoAss(SSNi),1);
           sdfs=3;
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
      for ssss=1:30
          
%           if ssss==99
%               fsdf=43;
%           end
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
                  
%                 CylBase=   ConvT_CylBase(ScafNoAss(kk),:) ; 
% 
%                 SameCyl= find(ConvT_CylBase(:,1)==CylBase(1) ) ; %subs later
% 
%                 d= (ConvT_CylBase(SameCyl,2)- CylBase(2)).^2 ;
%                 d(d==0)=50000; [Sortd,~] = sort(d);
%                 
%                 
%                closest= Sortd(1:2); 
%                
%                closest=closest( randi([1,length(closest)]));
% 
%                ConvT_CylBase( SameCyl(closest) ,:) ;

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
%     TransVec=[10,0,0]
%      Gcenter=[mean(pHss{1}.XData(Inscaf)),mean(pHss{1}.YData(Inscaf)) ,mean(pHss{1}.ZData(Inscaf)) ]
%     MovePlot( fH,pHss,BelongTransM,TStrand,T,AsswhichTM,TransVec )
  set(fH,'KeyPressFcn',{@(src,evn)keyMove(src,evn,fH,pHss,BelongTransM,TStrand,T ,popupH ,checkH ,editH)  }) ;
 
   popupH.Callback=@(src,evn) popupFcn(src,evn,pHss, fH ,BelongTransM,TStrand,T   );
               
   xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
  btn1 = uicontrol('Style', 'pushbutton', 'String', 'Export','Unit','normalized', 'Position', [0.8 0.1 0.1 0.1] ,...
            'Callback', {@(src,evn)ExportNewConf(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable ,PathName2  )});       
  btn2 = uicontrol('Style', 'pushbutton', 'String', 'Assembly','Unit','normalized', 'Position', [0.8 0.2 0.1 0.1] ,...
            'Callback', {@(src,evn)DoAssemble(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable   )});       
        
  btn3 = uicontrol('Style', 'pushbutton', 'String', 'AssemblyByFM','Unit','normalized', 'Position', [0.8 0.3 0.1 0.1] ,...
            'Callback', {@(src,evn)FMAssembly(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable   )});       
  btn4 = uicontrol('Style', 'pushbutton', 'String', 'Absor','Unit','normalized', 'Position', [0.8 0.4 0.1 0.1] ,...
            'Callback', {@(src,evn)AbsorAssembly(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable,CylBaseBynumIndex   )});       
  btn5 = uicontrol('Style', 'pushbutton', 'String', 'Fmin','Unit','normalized', 'Position', [0.91 0.2 0.1 0.1] ,...
            'Callback', {@(src,evn)FminconAssembly(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable,CylBaseBynumIndex   )});       
  btn6 = uicontrol('Style', 'pushbutton', 'String', 'Orthogonal','Unit','normalized', 'Position', [0.91 0.3 0.1 0.1] ,...
            'Callback', {@(src,evn)Orthog(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable,CylBaseBynumIndex,popupH   )});       
for k=1:length(pHss)
        pHss{k}.ButtonDownFcn=@(src,evn)lineselect(src,evn,popupH, pHss , fH ,BelongTransM,TStrand,T   );
end
  checkH_View = uicontrol('Style', 'checkbox','String', 'Axis auto/equal','Unit','normalized','Position', [0.8 0.55 0.05 0.05]); 
  checkH_View.Callback=@(src,evn)checkFcn2(src,evn,aH)  ; 
        
%      fileTransInd_name='BelongTransM.dat' ;
%     save(fileTransInd_name, 'BelongTransM' ) ;


sld1 = uicontrol('Style', 'slider','Parent',fH,'Units','normalized',....
'Min',-3,'Max',3,'Value',0,'Position', [0.8 0.72 0.15 0.02]);   
sld2 = uicontrol('Style', 'slider','Parent',fH,'Units','normalized',....
'Min',0.005,'Max',5,'Value',1,'Position', [0.8 0.75 0.15 0.02]);    

sld1.Callback= @(src,evn)sldscale1(src,evn,pHss,BelongTransM) ;

  sld2.Callback= @(src,evn)sldscale2(src,evn,pHss,BelongTransM) ;


sld1.SliderStep=[0.002,0.05];


    fileTransInd_name=strcat(PathName2, 'BM.mat') ;
     save(fileTransInd_name, 'BelongTransM') ;
%     fileTransInd_name=strcat(PathName2, 'CBB.mat') ;  %CylBaseBynumIndex
%      save(fileTransInd_name, 'CylBaseBynumIndex') ;
%    
     
%      scatter3(T(:,1),T(:,2),T(:,3))
     
%      figure
%      CC=0.6;
%      scatter3(T(:,1)+CC*T(:,4),T(:,2)+CC*T(:,5),T(:,3)+CC*T(:,6),'r')
%     axis equal   
    

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
    %----------
    % if scaffold strand ==1
    stapleStrand= unique(TStrand(ScafStapCorrelate2(:,2)+1))' ;
    ClosingStrand= setdiff(setdiff(1:max(TStrand) , stapleStrand) ,1 ) ;
    %---------
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
     fprintf('number of dsRemain= %i \n' , iF)
     %-------------
          file3_name='dSRemainNoClosing.conf'  
     fileID2 = fopen(file3_name,'w');
     cc3=0;
     for iF=1:size(ScafStapCorrelate2,1)
         if  ~ismember(TStrand(ScafStapCorrelate2(iF,1)+1),ClosingStrand) && ~ismember(TStrand(ScafStapCorrelate2(iF,2)+1),ClosingStrand)
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
            cc3=cc3+1 ;
         end
     end

     fclose(fileID2);
          fprintf('number of dsRemainNoClosing= %i \n' , cc3)

     
     
     
     
     cd(Spwd);
% %------end of new in V5
    
%     figure(222)
%     clf;
%     M3plotInd= BelongTransM==3;
%     M7plotInd= BelongTransM==7;
%     plot3(T(M3plotInd,1),T(M3plotInd,2),T(M3plotInd,3),'r');
%     hold on;
%      plot3(T(M7plotInd,1),T(M7plotInd,2),T(M7plotInd,3),'b')
%        sdgdg=45

%
%--------------------------
    function sldscale1(src,evn,pHss,BelongTransM)
    %change scaffold color
    sdf=3;
%      arrX= linspace(0.01,1, size(pHss{1}.XData,2)) ;
      arrX= linspace(0,0, size(pHss{1}.XData,2)) ;
     n=2^(src.Value) ;
     arrY=arrX.^n;
     pHss{1}.CData = [arrY; arrY] ;
    
    
    end

    function sldscale2(src,evn,pHss,~)
    %change staple linewidth
        for ph2=2:length(pHss)
            
         pHss{ph2}.LineWidth= src.Value;
         if src.Value<=0.02
          pHss{ph2}.LineWidth= src.Value;
            sdsf=34;
            pHss{ph2}.Visible='off';
         else
             pHss{ph2}.Visible='on'; 
         end
         
         sdfs=3;
         if  mean(pHss{ph2}.ZData)<221
         pHss{ph2}.Color=[0,0,0];
         else
           pHss{ph2}.Color=[1,0,0];   
         end
        end
        pHss{1}.FaceColor='black' ;
        
%          pHss{1}.Visible='off' ;
% %          ppp=plot3( pHss{1}.XData(1,:),pHss{1}.YData(1,:),pHss{1}.ZData(1,:),'k' )
%         
%         sdfs=3 
%         h = findobj(gca,'Type','line') ;
%         for k=1:length(h)
%             if sum( h(1).Color==[0,0,0])==3
%            h(k).LineWidth= 0.02;
%            
%             end
%         end
%         
        
    end


    function lineselect(src,evn,popupH,pHss , fH ,BelongTransM,TStrand,T )
  
        xy=evn.IntersectionPoint(1:3);
        dsfsf=3;
        XYZAll=[src.XData(1,:) ;src.YData(1,:)  ;   src.ZData(1,:)]' ;
%              XYZAll=[pHss{1}.XData( ;pHss{1}.YData  ;   pHss{1}.ZData]' ;
   
        dXY= XYZAll-  ones(size(XYZAll,1),1)*xy ;
        ds=dXY(:,1).^2 +dXY(:,2).^2  +dXY(:,3).^2;
        Ind= find(ds==min(ds)) ; Ind=Ind(1);
        BelongTransM(Ind);
         BelongTransM(src.UserData(Ind));
        popupH.Value=  BelongTransM(src.UserData(Ind)) ;
%         popupFcn(popupH,[],pH, staplePlotAx2, fH ,BelongTransM,ax2  );
        popupFcn(popupH,[],pHss, fH ,BelongTransM,TStrand,T  )
%         popupFcn(src,evn,pHss, fH ,BelongTransM,TStrand,T  )
    end


    function Orthog(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand,StramdIndTable,CylBaseBynumIndex ,popupH  )

       
        RefBundle=popupH.Value;
        nBundle2=max(BelongTransM);
        LastScaf=sum(TStrand==1) ;
        scaffolBelongM=BelongTransM(1:LastScaf);
        GlobalInd=scaffolBelongM==RefBundle ;
        
        CurrentXYZInScaf= [ pHss{1}.XData ;  pHss{1}.YData ; pHss{1}.ZData ]';  %before moving anything
        CurrentXYZScafCurBundle= [ pHss{1}.XData(GlobalInd) ;  pHss{1}.YData(GlobalInd) ; pHss{1}.ZData(GlobalInd) ]';  %before moving anything
        Ori=T(GlobalInd,1:3) ;
        
%         CurrentXYZ=
        
        
        [regParamsA,~,Acc2]=absor(CurrentXYZScafCurBundle',Ori')     ;  
        theta=pi/2;
        RMat=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];

        OrthM= RMat*regParamsA.R ;
        
        CurrentXYZInScaf=(OrthM*CurrentXYZInScaf')' ;
        pHss{1}.XData= CurrentXYZInScaf(:,1);
        pHss{1}.YData= CurrentXYZInScaf(:,2);
        pHss{1}.ZData= CurrentXYZInScaf(:,3);
        
        for k222=2:length(pHss)
                         origXYZ=[ pHss{k222}.XData ;  pHss{k222}.YData ; pHss{k222}.ZData ]' ;
                          N1XYZ=  (OrthM *origXYZ'  )' ;
                           pHss{k222}.XData=N1XYZ(:,1); 
                           pHss{k222}.YData=N1XYZ(:,2); 
                           pHss{k222}.ZData=N1XYZ(:,3); 
              
        end
        
         sdfsf=234;
        
        
    end

    function FminconAssembly(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand ,StramdIndTable ,CylBaseBynumIndex )
         ScaffoldChain=find(StramdIndTable(1,:)==max( StramdIndTable(1,:)));
%         sdfsf22=234
        nBundle2=max(BelongTransM);
        LastScaf=sum(TStrand==1) ;
        scaffolBelongM=BelongTransM(1:LastScaf);
        IndForceConnect= find(BelongTransM(1:LastScaf-1)~= BelongTransM(2:LastScaf)) ;
        TableList=[BelongTransM(IndForceConnect) ,BelongTransM(IndForceConnect+1) ,IndForceConnect];
        
        BundleGraphEdg= unique(TableList(:,1:2),  'rows');
        IndFlip=BundleGraphEdg(:,1)>BundleGraphEdg(:,2) ;
        BundleGraphEdg(IndFlip,:) = flip(  BundleGraphEdg(IndFlip,:) ,2);
        BundleGraphEdg=unique(BundleGraphEdg,'rows') ;  
        
        nEdge=size(BundleGraphEdg,1) ;
        TableList2=[TableList , TableList(:,3)+1] ;
        IndFlip2=TableList2(:,1)>TableList2(:,2) ;
        TableList2(IndFlip2,1:2) =[ TableList2(IndFlip2,2) , TableList2(IndFlip2,1) ];
        TableList2(IndFlip2,3:4) =[ TableList2(IndFlip2,4) , TableList2(IndFlip2,3) ];
        
        EffectTable=ones(size(TableList2,1),1) ;
        for iT=1:length(EffectTable)
            LeftP =[ TableList2(iT,3)];
            RightP =[ TableList2(iT,4)];
            LP_CB=CylBaseBynumIndex(LeftP,:) ;
            RP_CB= CylBaseBynumIndex(RightP,:) ;
            point=0;
            if  length(intersect( find(CylBaseBynumIndex(:,1)==LP_CB(1)) ,find(CylBaseBynumIndex(:,2)==LP_CB(2)+2 )))==2 ; point=point+1;end
            if  length(intersect( find(CylBaseBynumIndex(:,1)==LP_CB(1)) ,find(CylBaseBynumIndex(:,2)==LP_CB(2)-2 )))==2 ; point=point+1;end
            if  length(intersect( find(CylBaseBynumIndex(:,1)==RP_CB(1)) ,find(CylBaseBynumIndex(:,2)==RP_CB(2)+2 )))==2 ; point=point+1;end
            if  length(intersect( find(CylBaseBynumIndex(:,1)==RP_CB(1)) ,find(CylBaseBynumIndex(:,2)==RP_CB(2)-2 )))==2 ; point=point+1;end
            if point==0 ;  EffectTable(iT)=0; end
        end
        TableList2= TableList2(EffectTable==1   ,:) ;
         CurrentXYZInScaf= [ pHss{1}.XData ;  pHss{1}.YData ; pHss{1}.ZData ]';  %before moving anything
%         x0=zeros(1, 12*nBundle2 ); 
%         x0=reshape([eye(3), [0 ;0;0]] ,1,12);
        x0=[0,0,0,0,0,0];
        x0= repmat(x0,[1,nBundle2]);
        
%         EachM=reshape(x0,[3,4,nBundle2]) ;
        
        Aine = [];
        bine = [];
        Aeq = [];
        beq = [];
        fun=@(x)CalculateForceEnergy(x, TableList2, CurrentXYZInScaf) ;
%         nonlcon = @RotationMCons;
%         vv=reshape(-1*ones(1,12*nBundle2) ,
        Bmm=200;
        lb =-pi*ones(1,6*nBundle2) ; lb(4:6:end)=-Bmm ; lb(5:6:end)=-Bmm ; lb(6:6:end)=-Bmm ;
        ub = pi*ones(1,6*nBundle2) ; ub(4:6:end)=Bmm ; ub(5:6:end)=Bmm ; ub(6:6:end)=Bmm ;
        
%        ub = [0.5,0.8];
       options = optimoptions('fmincon','Display','iter','Algorithm','interior-point','MaxFunctionEvaluations',15000  );
        [x,fval]= fmincon(fun,x0,Aine,bine,Aeq,beq,lb,ub,[],options) ;
        fval
%         EachM=reshape(x,[3,4,nBundle2]) ;
         EachV= reshape(x,[1,6,nBundle2]) ;
         EachM=zeros(4,4,nBundle2);
         for kt2=1:nBundle2
         EachM(:,:,kt2)=convertEulerToRotation(EachV(:,:,kt2))  ;
         end
        
        for kB2=1:nBundle2
        indXX= scaffolBelongM==kB2 ;
        CurrentXYZInScaf(indXX,:) =  (EachM(1:3,1:3 ,kB2) *CurrentXYZInScaf(indXX,:)' +  EachM(1:3,4 ,kB2) )' ;
        end
%         [c,ceq] = RotationMCons(x)
        tol=1e-4;
%         if  abs(min(ceq))<tol &&  abs(max(ceq))<tol
        
         pHss{1}.XData= CurrentXYZInScaf(:,1) ;
         pHss{1}.YData= CurrentXYZInScaf(:,2) ;
         pHss{1}.ZData= CurrentXYZInScaf(:,3) ;
              for bundi=1: nBundle2
                  for k1=2:length(pHss)
                      stapGlobalInd=pHss{k1}.UserData ;
                      IndEff=  find(BelongTransM(stapGlobalInd)==bundi );
                      if ~isempty(IndEff)
                          origXYZ=[ pHss{k1}.XData(IndEff) ;  pHss{k1}.YData(IndEff) ; pHss{k1}.ZData(IndEff) ]' ;
                          N1XYZ=  (EachM(1:3,1:3 ,bundi) *origXYZ' +  EachM(1:3,4 ,bundi) )' ;
                           pHss{k1}.XData(IndEff)=N1XYZ(:,1); 
                           pHss{k1}.YData(IndEff)=N1XYZ(:,2); 
                           pHss{k1}.ZData(IndEff)=N1XYZ(:,3); 
                      end
                  end         
              end
%         else
%             notsatified=1
%         end
        
        sdfsf=234;
    end

    function  T=convertEulerToRotation(Vars_6)  
        
           RMatx=[1 0 0; 0 cos(Vars_6(1)) sin(Vars_6(1)); 0 -sin(Vars_6(1)) cos(Vars_6(1))];
           RMaty=[cos(Vars_6(2)) 0 -sin(Vars_6(2)); 0 1 0; sin(Vars_6(2)) 0 cos(Vars_6(2))];
           RMatz=[cos(Vars_6(3)) sin(Vars_6(3)) 0; -sin(Vars_6(3))  cos(Vars_6(3))  0; 0 0 1];     
        M=RMatx*RMaty*RMatz;
        T=[M, [Vars_6(4);Vars_6(5);Vars_6(6)] ; 0 0 0 1];
    end
%     function [c,ceq] = RotationMCons(x)
%       AllM = reshape(x, [3,4,  length(x)/12]) ;  
%        ceq=zeros( length(x)/2  ,1);
%        k_con=1;
%        for buni=1:length(x)/12
%            Mi=AllM(:,:,buni);
%          ceq(k_con)=  norm( [Mi(1,1) ,Mi(2,1),Mi(3,1)]) -1 ; k_con=k_con+1 ;
%          ceq(k_con)=  norm( [Mi(1,2) ,Mi(2,2),Mi(3,2)]) -1 ; k_con=k_con+1 ;
%          ceq(k_con)= dot(Mi(:,1), Mi(:,2))  ; k_con=k_con+1 ;
%          CrossUV=cross(Mi(:,1), Mi(:,2));
%          ceq(k_con)= CrossUV(1)-Mi(1,3)   ; k_con=k_con+1 ;
%          ceq(k_con)= CrossUV(2)-Mi(2,3)   ; k_con=k_con+1 ;
%          ceq(k_con)= CrossUV(3)-Mi(3,3)   ; k_con=k_con+1 ;
%        end
%     
%     c = [];
%     end

    function E=CalculateForceEnergy(x, TableList2, CurrentXYZInScaf)
        
%         sdf=24;
         E=0;   nBundle3=length(x)/6  ;
%          EachM=reshape(x,[3,4,nBundle3]) ;
         EachV= reshape(x,[1,6,nBundle3]) ;
         EachM=zeros(4,4,nBundle3);
         for kt=1:nBundle3
         EachM(:,:,kt)=convertEulerToRotation(EachV(:,:,kt))  ;
         end
         
        for k3=1:size(TableList2,1)
        TransL=EachM(:,:,TableList2(k3,1) ) ; P_L= CurrentXYZInScaf(TableList2(k3,3),:) ;
        TransR=EachM(:,:,TableList2(k3,2) ) ; P_R= CurrentXYZInScaf(TableList2(k3,4),:) ;
        New_PL=TransL(1:3,1:3)*P_L' + TransL(1:3,4) ;
        New_PR=TransR(1:3,1:3)*P_R' + TransR(1:3,4) ;
        
        E_k3= norm(New_PL-New_PR) ;
        E_k3= (E_k3 - 0.6)^2 ;  
        E= E + E_k3 ;
        end
    end

    function AbsorAssembly(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand ,StramdIndTable ,CylBaseBynumIndex )
        for ccc=1:50
            pause(0.1);
            ccc
         ScaffoldChain=find(StramdIndTable(1,:)==max( StramdIndTable(1,:)));
%         sdfsf22=234
        ScafoldBundleCenter=zeros(max(BelongTransM),3) ;
         alpha=-0.1;Amp=0.1;  alpha2=0.01;Amp2=5;
        nBundle2=max(BelongTransM);
        LastScaf=sum(TStrand==1) ;
        scaffolBelongM=BelongTransM(1:LastScaf);
        IndForceConnect= find(BelongTransM(1:LastScaf-1)~= BelongTransM(2:LastScaf)) ;
        TableList=[BelongTransM(IndForceConnect) ,BelongTransM(IndForceConnect+1) ,IndForceConnect];
        
        BundleGraphEdg= unique(TableList(:,1:2),  'rows');
        IndFlip=BundleGraphEdg(:,1)>BundleGraphEdg(:,2) ;
        BundleGraphEdg(IndFlip,:) = flip(  BundleGraphEdg(IndFlip,:) ,2);
        BundleGraphEdg=unique(BundleGraphEdg,'rows') ;  
        
        nEdge=size(BundleGraphEdg,1) ;
        TableList2=[TableList , TableList(:,3)+1] ;
        IndFlip2=TableList2(:,1)>TableList2(:,2) ;
        TableList2(IndFlip2,1:2) =[ TableList2(IndFlip2,2) , TableList2(IndFlip2,1) ];
        TableList2(IndFlip2,3:4) =[ TableList2(IndFlip2,4) , TableList2(IndFlip2,3) ];
        
        EffectTable=ones(size(TableList2,1),1) ;
        for iT=1:length(EffectTable)
            LeftP =[ TableList2(iT,3)];
            RightP =[ TableList2(iT,4)];
            LP_CB=CylBaseBynumIndex(LeftP,:) ;
            RP_CB= CylBaseBynumIndex(RightP,:) ;
            point=0;
            if  length(intersect( find(CylBaseBynumIndex(:,1)==LP_CB(1)) ,find(CylBaseBynumIndex(:,2)==LP_CB(2)+2 )))==2 ; point=point+1;end
            if  length(intersect( find(CylBaseBynumIndex(:,1)==LP_CB(1)) ,find(CylBaseBynumIndex(:,2)==LP_CB(2)-2 )))==2 ; point=point+1;end
            if  length(intersect( find(CylBaseBynumIndex(:,1)==RP_CB(1)) ,find(CylBaseBynumIndex(:,2)==RP_CB(2)+2 )))==2 ; point=point+1;end
            if  length(intersect( find(CylBaseBynumIndex(:,1)==RP_CB(1)) ,find(CylBaseBynumIndex(:,2)==RP_CB(2)-2 )))==2 ; point=point+1;end
            if point==0 ;  EffectTable(iT)=0; end
            
        sdff=424;
        end
        TableList2= TableList2(EffectTable==1   ,:) ;
        
        
        
         CurrentXYZInScaf= [ pHss{1}.XData ;  pHss{1}.YData ; pHss{1}.ZData ]';  %before moving anything
        %----------
%         
%         sumE_edge=zeros(size(BundleGraphEdg,1),1) ;
%         for cal_Eng=1:size(TableList2,1)
%        [~, EdgeInd]= ismember(TableList2(cal_Eng,1:2) ,BundleGraphEdg ,'rows') ;
%        PA= CurrentXYZInScaf( TableList2(cal_Eng,3),:)  ;
%        PB= CurrentXYZInScaf( TableList2(cal_Eng,4),:)  ;
%        Ei= norm(PA-PB) ;
% %      
%          sumE_edge(  EdgeInd)  =sumE_edge(  EdgeInd)   +Ei ;
%         end
%         %----------
%         AllEnergy= sum(sumE_edge)
%                 randEdge= BundleGraphEdg( 1  ,:)   ;
% randEdge=[3,4];
%   sdfsf=234;
%         randEdge= BundleGraphEdg( randi(nEdge)  ,:)  ; 
%         randEdge=BundleGraphEdg(sumE_edge==max(sumE_edge),:);
       randEdge= BundleGraphEdg( mod(ccc,nEdge)+1  ,:)  ; 
        
%        if sum(randEdge==[3,4])==2
%            sdfsggg=3
%        end
       
        TableInd= ismember(TableList2(:,1:2), randEdge ,'rows') ;
        
        GraphEdge=graph(BundleGraphEdg(:,1) ,BundleGraphEdg(:,2));
        Cen = centrality(GraphEdge,'closeness') ;
        
        LeftXYZ=  CurrentXYZInScaf(TableList2(TableInd,3),:) ;
        RightXYZ=  CurrentXYZInScaf(TableList2(TableInd,4),:) ;
  
        if Cen(randEdge(1))< Cen(randEdge(2))
            ccCase=1;
        elseif  Cen(randEdge(1))> Cen(randEdge(2))
           ccCase=2;
        else
           ccCase=1;
        end
        
        switch ccCase
            case 1
                [regParams1,~,Acc2]=absor(LeftXYZ',RightXYZ')     ;   
%                  regParams1.R = sqrtm(regParams1.R) ;  %----------
                 
                LeftBundle= unique(TableList2(TableInd,1));
%                 if LeftBundle==1
%                     continue;
%                 end
%                 
                
                RightBundle=setdiff(randEdge,LeftBundle);
                IndScafLeft=scaffolBelongM==LeftBundle;
                MovingSacf=CurrentXYZInScaf(IndScafLeft,:) ;
                MovingSacf=regParams1.R*MovingSacf' +  regParams1.t ;
                CurrentXYZInScaf(IndScafLeft,:)=MovingSacf';   %use absor fix to hinge axis, update        
                 xyz=RightXYZ;
                 r0=mean(xyz);
                 xyz=bsxfun(@minus,xyz,r0);
                 [~,~,V]=svd(xyz,0);
                 IndLeft=union(find(TableList2(:,1)==LeftBundle),find(TableList2(:,2)==LeftBundle)) ;
                 LeftAtLeft=TableList2(IndLeft,1)==LeftBundle;
                 chainL0= [CurrentXYZInScaf( TableList2(IndLeft(find(LeftAtLeft==1)) ,3)  ,:) ; CurrentXYZInScaf( TableList2(IndLeft(find(LeftAtLeft==0)) ,4)  ,:)] ;
                 chainR0= [CurrentXYZInScaf( TableList2(IndLeft(find(LeftAtLeft==1)) ,4)  ,:) ; CurrentXYZInScaf( TableList2(IndLeft(find(LeftAtLeft==0)) ,3)  ,:)] ;
                 %----------          
                 theta = 0.1:5:360 ;
                 Eng=zeros(size(theta));
                 for k_theta=1:length(theta)
                      M_Trans_theta=AxelRot(theta(k_theta), V(:,1)',r0) ;
                      k_ChainL= M_Trans_theta(1:3,1:3)*chainL0' + M_Trans_theta(1:3,4) ;
                      chainL_k=k_ChainL'  ;
                      DXYZ_k=chainL_k-chainR0;
                       Eng(k_theta)= sumsqr(DXYZ_k);
                 end
                 BestAng= theta((Eng==min(Eng)));
                 %---------------
                 M_Trans_Best=AxelRot(BestAng(1), V(:,1)',r0) ;
                 QQ2= M_Trans_Best(1:3,1:3)*CurrentXYZInScaf(IndScafLeft,:)' + M_Trans_Best(1:3,4) ;
                 CurrentXYZInScaf(IndScafLeft,:)=QQ2';

                  pHss{1}.XData=CurrentXYZInScaf(:,1);
                  pHss{1}.YData=CurrentXYZInScaf(:,2);
                  pHss{1}.ZData=CurrentXYZInScaf(:,3);
                  for k1=2:length(pHss)
                      stapGlobalInd=pHss{k1}.UserData ;
                      IndEff=  find(BelongTransM(stapGlobalInd)==LeftBundle );
                      if ~isempty(IndEff)
                          origXYZ=[ pHss{k1}.XData(IndEff) ;  pHss{k1}.YData(IndEff) ; pHss{k1}.ZData(IndEff) ]' ;
                          N1XYZ=  regParams1.R* origXYZ' +  regParams1.t ;
                          N2XYZ= M_Trans_Best(1:3,1:3)*N1XYZ + M_Trans_Best(1:3,4) ;
                           pHss{k1}.XData(IndEff)=N2XYZ(1,:); 
                           pHss{k1}.YData(IndEff)=N2XYZ(2,:); 
                           pHss{k1}.ZData(IndEff)=N2XYZ(3,:); 
                      end
                  end
             
            case 2
                [regParams1,~,Acc2]=absor(RightXYZ',LeftXYZ')     ;   
%                 regParams1.R = sqrtm(regParams1.R) ; %-------
                RightBundle= unique(TableList2(TableInd,2));
%                 if RightBundle==1
%                     continue;
%                 end

                
                
                LeftBundle=setdiff(randEdge,RightBundle);
                IndScafRight=scaffolBelongM==RightBundle;
                MovingSacf=CurrentXYZInScaf(IndScafRight,:) ;
                MovingSacf=regParams1.R*MovingSacf' +  regParams1.t ;
                CurrentXYZInScaf(IndScafRight,:)=MovingSacf';   %use absor fix to hinge axis, update        
                 xyz=LeftXYZ;
                 r0=mean(xyz);
                 xyz=bsxfun(@minus,xyz,r0);
                 [~,~,V]=svd(xyz,0);
                 
                 IndRight=union(find(TableList2(:,2)==RightBundle),find(TableList2(:,1)==RightBundle)) ;
                 RightAtRight=TableList2(IndRight,1)==RightBundle;
                 chainR0= [CurrentXYZInScaf( TableList2(IndRight(find(RightAtRight==1)) ,3)  ,:) ; CurrentXYZInScaf( TableList2(IndRight(find(RightAtRight==0)) ,4)  ,:)] ;
                 chainL0= [CurrentXYZInScaf( TableList2(IndRight(find(RightAtRight==1)) ,4)  ,:) ; CurrentXYZInScaf( TableList2(IndRight(find(RightAtRight==0)) ,3)  ,:)] ;
                 %----------          
                 theta = 0.1:1:360 ;
                 Eng=zeros(size(theta));
                 for k_theta=1:length(theta)
                      M_Trans_theta=AxelRot(theta(k_theta), V(:,1)',r0) ;
                      k_ChainL= M_Trans_theta(1:3,1:3)*chainR0' + M_Trans_theta(1:3,4) ;
                      chainL_k=k_ChainL'  ;
                      DXYZ_k=chainL_k-chainL0;
                       Eng(k_theta)= sumsqr(DXYZ_k);
                 end
                 BestAng= theta((Eng==min(Eng)));
                 %---------------
                 M_Trans_Best=AxelRot(BestAng(1), V(:,1)',r0) ;
                 
%                  M_Trans_Best=eye(4)
                 
                 QQ2= M_Trans_Best(1:3,1:3)*CurrentXYZInScaf(IndScafRight,:)' + M_Trans_Best(1:3,4) ;
                 CurrentXYZInScaf(IndScafRight,:)=QQ2';

                  pHss{1}.XData=CurrentXYZInScaf(:,1);
                  pHss{1}.YData=CurrentXYZInScaf(:,2);
                  pHss{1}.ZData=CurrentXYZInScaf(:,3);
                  for k1=2:length(pHss)
                      stapGlobalInd=pHss{k1}.UserData ;
                      IndEff=  find(BelongTransM(stapGlobalInd)==RightBundle );
                      if ~isempty(IndEff)
                          origXYZ=[ pHss{k1}.XData(IndEff) ;  pHss{k1}.YData(IndEff) ; pHss{k1}.ZData(IndEff) ]' ;
                          N1XYZ=  regParams1.R* origXYZ' +  regParams1.t ;
                          N2XYZ= M_Trans_Best(1:3,1:3)*N1XYZ + M_Trans_Best(1:3,4) ;
                           pHss{k1}.XData(IndEff)=N2XYZ(1,:); 
                           pHss{k1}.YData(IndEff)=N2XYZ(2,:); 
                           pHss{k1}.ZData(IndEff)=N2XYZ(3,:); 
                      end
                  end
                
                
                
                
        end
%           pHss{1}.XData(IndUpdate)= NewXYZAll(1,:)  +dX(jBun,1) 
        
        
        end
        
    end
%---------------------------------------------------------------------------------------------------------
    function FMAssembly(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand ,StramdIndTable  )
        
         ScaffoldChain=find(StramdIndTable(1,:)==max( StramdIndTable(1,:)));
%         sdfsf22=234
        ScafoldBundleCenter=zeros(max(BelongTransM),3) ;
         alpha=-0.05;Amp=0.1;  alpha2=0.05;Amp2=0.1;
        nBundle2=max(BelongTransM);
        LastScaf=sum(TStrand==1) ;
        
%         %find G-Center
%         for k1=1:max(BelongTransM)
%             IndBundleScaf=BelongTransM==k1 ;
%             IndBundleScaf( sum(TStrand==1)+1 :end)=0;
%             
%             ScafoldBundleCenter(k1,1)= mean(pHss{ScaffoldChain}.XData(IndBundleScaf)) ;
%             ScafoldBundleCenter(k1,2)= mean(pHss{ScaffoldChain}.YData(IndBundleScaf)) ;
%             ScafoldBundleCenter(k1,3)= mean(pHss{ScaffoldChain}.ZData(IndBundleScaf)) ;
% 
% %             sdfsf=234
%         end
%         sdfsf=234
        scaffolBelongM=BelongTransM(1:LastScaf);



        MaxnLoop=50; nloop2=1;  BackbondThres=1.2 ;
         while nloop2<MaxnLoop 
             
                     if nloop2<20
                     Amp2=Amp2*(nloop2+54)/(nloop2+55);
                     Amp=Amp*(nloop2+81)/(nloop2+82);
                     end
             
              dX=zeros(nBundle2,3); %3D translation of each iteration, like summation of all force one each bundle
              dM=zeros(nBundle2,3); %3D Rotation of each iteration
%               XYZ=[ pHss{ScaffoldChain}.XData ; pHss{ScaffoldChain}.YData ; pHss{ScaffoldChain}.ZData];
%               dxyz= XYZ(:,2:end)- XYZ(:,1:end-1);
%               backbondD=sqrt(dxyz(1,:).^2+dxyz(2,:).^2+dxyz(3,:).^2) ;
%               IndForceConnect=find(backbondD>BackbondThres) ;       
            IndForceConnect= find(BelongTransM(1:LastScaf-1)~= BelongTransM(2:LastScaf)) ;
            TableList=[BelongTransM(IndForceConnect) ,BelongTransM(IndForceConnect+1) ,IndForceConnect];
            %[ bundle i,  bundle j  , Dist-Index= Pi+1 -P1  ]
            for iT=1:size(TableList,1)
                %--------translation
                PA=[pHss{1}.XData(TableList(iT,3)) , pHss{1}.YData(TableList(iT,3)), pHss{1}.ZData(TableList(iT,3)) ];
                PB=[pHss{1}.XData(1+TableList(iT,3)) , pHss{1}.YData(1+TableList(iT,3)), pHss{1}.ZData(1+TableList(iT,3)) ];
                distance=norm(PA-PB)  ;
                Ydist=0.1*Amp*(1-exp(alpha*distance))  ;
                VV=(PA-PB)/(norm((PA-PB)) -10) ;
                 XX=Ydist*VV*distance/3;   %  internal force
                 dX(TableList(iT,1),:)= dX(TableList(iT,1),:) -XX;
                 dX(TableList(iT,2),:)=  dX(TableList(iT,2),:) +XX; 
%                 sdf=2;
                %----------rotation
%                 ddgdfg=2
                %----each bundle current center
                for k1=1:max(BelongTransM)
                    IndBundleScaf=scaffolBelongM==k1 ;
                    ScafoldBundleCenter(k1,1)= mean(pHss{ScaffoldChain}.XData(IndBundleScaf)) ;
                    ScafoldBundleCenter(k1,2)= mean(pHss{ScaffoldChain}.YData(IndBundleScaf)) ;
                    ScafoldBundleCenter(k1,3)= mean(pHss{ScaffoldChain}.ZData(IndBundleScaf)) ;
                end 
                %
                GCcenterLeft= ScafoldBundleCenter(TableList(iT,1) ,:);
                GCcenterRight= ScafoldBundleCenter(TableList(iT,2) ,:);
                  RVecLeft=PA-GCcenterLeft;
                  RVecRight=PB-GCcenterRight;
                OneMommentOnLeft=cross(RVecLeft,-XX);
                OneMommentOnRight=cross(RVecRight,XX);
                 %----transform Moment to small angular rotation
                 SmallRotationRight=zeros(1,3);
                 SmallRotationLeft=zeros(1,3);
                        for signi=1:3
                            if OneMommentOnRight(signi)>0
                             SmallRotationRight(signi)=Amp2*(1-exp(alpha2*OneMommentOnRight(signi) )); 
                            else 
                             SmallRotationRight(signi)=-Amp2*(1-exp(-alpha2*OneMommentOnRight(signi) ));        
                            end
                            if OneMommentOnLeft(signi)>0
                             SmallRotationLeft(signi)=Amp2*(1-exp(alpha2*OneMommentOnLeft(signi) )); 
                            else 
                             SmallRotationLeft(signi)=-Amp2*(1-exp(-alpha2*OneMommentOnLeft(signi) ));        
                            end 
                        end        
                        dM(TableList(iT,1),:)= dM(TableList(iT,1),:) +SmallRotationLeft;
                        dM(TableList(iT,2),:)= dM(TableList(iT,2),:) +SmallRotationRight;    
                        
%                  sdfsf=234    ;   
                        
                        
                
            end
            
            %---------------
            
%             sdfsdf=3;
            
            
           %%% --------------scaf update
            for jBun=1:nBundle2
               IndUpdate= scaffolBelongM==jBun ;
               AllXYZ=[pHss{1}.XData(IndUpdate); pHss{1}.YData(IndUpdate) ;pHss{1}.ZData(IndUpdate)] ;
                        Rx=[1,0,0;0 cosd(dM(jBun,1)), sind(dM(jBun,1)); 0 -sind(dM(jBun,1)), cosd(dM(jBun,1))];
                        Ry=[cosd(dM(jBun,2)),0,-sind(dM(jBun,2)) ;0 ,1 ,0 ; sind(dM(jBun,2)), 0 ,cosd(dM(jBun,2))];
                        Rz=[cosd(dM(jBun,3)),sind(dM(jBun,3)),0 ;  -sind(dM(jBun,3)) cosd(dM(jBun,3))  0 ; 0 ,0 ,1];
                NewXYZAll=Rx*Ry*Rz*AllXYZ ;
                        
                        
                pHss{1}.XData(IndUpdate)= NewXYZAll(1,:)  +dX(jBun,1) ;
                pHss{1}.YData(IndUpdate)=  NewXYZAll(2,:) +dX(jBun,2) ;
                pHss{1}.ZData(IndUpdate)=  NewXYZAll(3,:) +dX(jBun,3) ;
%                   pHss{1}.XData(IndUpdate)= pHss{1}.XData(IndUpdate) +dX(jBun,1) ;
%                 pHss{1}.YData(IndUpdate)= pHss{1}.YData(IndUpdate) +dX(jBun,2) ;
%                 pHss{1}.ZData(IndUpdate)= pHss{1}.ZData(IndUpdate) +dX(jBun,3) ;
                
            end
            %------stap update
            
            for iBundle=1:nBundle2
                for kpHs=2:length(pHss)
                  BundleInd=  BelongTransM(pHss{kpHs}.UserData)==iBundle ;
                  AllXYZ2=[pHss{kpHs}.XData(BundleInd); pHss{kpHs}.YData(BundleInd) ;pHss{kpHs}.ZData(BundleInd)] ; 
                        Rx=[1,0,0;0 cosd(dM(iBundle,1)), sind(dM(iBundle,1)); 0 -sind(dM(iBundle,1)), cosd(dM(iBundle,1))];
                        Ry=[cosd(dM(iBundle,2)),0,-sind(dM(iBundle,2)) ;0 ,1 ,0 ; sind(dM(iBundle,2)), 0 ,cosd(dM(iBundle,2))];
                        Rz=[cosd(dM(iBundle,3)),sind(dM(iBundle,3)),0 ;  -sind(dM(iBundle,3)) cosd(dM(iBundle,3))  0 ; 0 ,0 ,1];
                NewXYZAll2=Rx*Ry*Rz*AllXYZ2 ;

             pHss{kpHs}.XData(BundleInd)=NewXYZAll2(1,:) + dX(iBundle,1) ;
             pHss{kpHs}.YData(BundleInd)=NewXYZAll2(2,:) + dX(iBundle,2) ;
             pHss{kpHs}.ZData(BundleInd)=NewXYZAll2(3,:) + dX(iBundle,3) ;
%                   pHss{kpHs}.XData=   pHss{kpHs}.XData + dX(iBundle,1)' ;
%                   pHss{kpHs}.YData=   pHss{kpHs}.YData + dX(iBundle,2)' ;
%                   pHss{kpHs}.ZData=   pHss{kpHs}.ZData + dX(iBundle,3)' ;
%                   sdfsf=2  ;
                end
            end
            
            
            
             nloop2=nloop2+1 ;
             drawnow
             pause(0.1)
         end  % end of while loop

        nloop2
        
        
        
        

    end

%------------
    function     DoAssemble(src,evn,ConfBB, fH,pHss,BelongTransM,T,TStrand ,StramdIndTable  )
 
        sdgfsdg=435 ;
          MovementCase={'q', 'w', 'e' ,'a', 's', 'd'};
          ScaffoldChain=find(StramdIndTable(1,:)==max( StramdIndTable(1,:)));
          nbetter=0;
        for k33=1:100
          
          chain1= [pHss{ScaffoldChain}.XData; pHss{ScaffoldChain}.YData ; pHss{ScaffoldChain}.ZData];
%           dxyz=chain(:,2:end) - chain(:,1:end-1) ;
          E1 = norm(chain1(:,2:end) - chain1(:,1:end-1)) ;
          
%         keyMove(src,evn,fH,pHss,BelongTransM,TStrand,T ,popupH ,checkH ,editH) 
        randBundle=randi(max(BelongTransM)) ; popupH2.Value=randBundle;
        randTransOrRotate=  randi(2);   % do translate first
        ccH.Value=0 ;
        
        
        Delta.String='0.5' ; % linear-> 5 nm,  , rotation0> 5 deg
        Delta2.String='-0.5' ; % linear-> -5 nm,  , rotation0> -5 deg
        
        
        evn2.Character=MovementCase{randi(6)} ;
        switch randTransOrRotate
            case 1  % translate
                 ccH.Value=0 ;
                 DDxyz=0.05*rand(1,3) ;
           Sx.String=num2str(  DDxyz(1)) ;  Sy.String=num2str(  DDxyz(2)) ;    Sz.String=num2str(  DDxyz(3)) ;     
           EE1x.Character='q' ;     EE2x.Character='w' ;       EE3x.Character='e' ;       
           keyMove([],EE1x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Sx) ;
           keyMove([],EE2x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Sy) ;
           keyMove([],EE3x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Sz) ;
%          keyMove([],evn2,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Delta) ;
            case 2
                 ccH.Value=1 ;
            RRxyz=0.05*rand(1,3) ;
           Rx.String=num2str(  RRxyz(1)) ;  Ry.String=num2str(  RRxyz(2)) ;    Rz.String=num2str(  RRxyz(3)) ;     
           EE1x.Character='q' ;     EE2x.Character='w' ;       EE3x.Character='e' ;       
            keyMove([],EE1x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Rx) ;
            keyMove([],EE2x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Ry) ;
            keyMove([],EE3x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Rz) ; 
%           keyMove([],evn2,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Delta) ;
        end
        
         chain2= [pHss{ScaffoldChain}.XData; pHss{ScaffoldChain}.YData ; pHss{ScaffoldChain}.ZData];
%           dxyz=chain(:,2:end) - chain(:,1:end-1) ;
          E2 = norm(chain2(:,2:end) - chain2(:,1:end-1)) ;

          if E2>E1
            switch randTransOrRotate
                case 1  % translate
                     ccH.Value=0 ;
              NDDxyz=-DDxyz;
           Sx.String=num2str(  NDDxyz(1)) ;  Sy.String=num2str(  NDDxyz(2)) ;    Sz.String=num2str(  NDDxyz(3)) ;    
           EE1x.Character='q' ;     EE2x.Character='w' ;       EE3x.Character='e' ;       
           keyMove([],EE1x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Sx) ;
           keyMove([],EE2x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Sy) ;
           keyMove([],EE3x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Sz) ;
%              keyMove([],evn2,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Delta2) ;
                case 2
                     ccH.Value=1 ;
              NRRxyz=-RRxyz;       
           Rx.String=num2str(  NRRxyz(1)) ;  Ry.String=num2str(  NRRxyz(2)) ;    Rz.String=num2str(  NRRxyz(3)) ;  
           EE1x.Character='q' ;     EE2x.Character='w' ;       EE3x.Character='e' ;    
             keyMove([],EE3x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Rz) ; 
             keyMove([],EE2x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Ry) ;
             keyMove([],EE1x,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Rx) ;
%               keyMove([],evn2,[],pHss,BelongTransM,[],[] ,popupH2 ,ccH ,Delta2) ;
            end
          else
              nbetter=nbetter+1 ;
          end
%           ssdfdfsf=234;
        k33;
        end
        [E1  , nbetter]
        
        
    end
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
            [regParams,q,w]=absor(A3byN,B3byN);
             w.errmax
    PosV=    regParams.R*A3byN + regParams.t*ones(1,size(B3byN,2)) ;
    BVecNew=    regParams.R*OldBVec ;
    NVecNew=    regParams.R*OldNVec ;
    
   NewT( IndexAll,1:9)=   [PosV'-Coeff2*BVecNew',   BVecNew',   NVecNew'];
   
   if sum(IndexAll([65,67,68]))==3
       [N,edges] = histcounts(TStrand,0.5:max(TStrand)+0.5)
       sdfsf=3
   end
   
   
end

if ScafInd~=1  %scaf strand is not the first in topology file------found bug in V3
    MoveAheadind= sum(StramdIndTable(2,1:ScafInd-1));
    
    NewT=[ NewT(ScafL+1:ScafL+MoveAheadind ,: ) ; NewT(1:ScafL,:); NewT(ScafL+MoveAheadind+1:end,:)];
    BelongTransM=[ BelongTransM(ScafL+1:ScafL+MoveAheadind ,: ) ; BelongTransM(1:ScafL,:); BelongTransM(ScafL+MoveAheadind+1:end,:)];
end

    fileTransInd_name=strcat(PathName2, 'BM.mat') ;
     save(fileTransInd_name, 'BelongTransM') ;





 mmNewT=min(NewT(:,1:3)) ;
 MMNewT=max(NewT(:,1:3)) ;

 boxsize= MMNewT-mmNewT +[30,30,30];
 
 
 NewT(:,1:3)= NewT(:,1:3) - ones(size(NewT,1),1)*( mmNewT-[60,60,60]) ;
 
file2_name='prova22_tt.conf'
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


        
        hweiuhf=1 ;
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
ssdf=34
                 for kstp=2: length(pHss)  
                    if ~isempty(intersect( pHss{kstp}.UserData,GindexNeedTomove))
                        sdfsfhhs=234;
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
                        pHss{kstp}.MarkerSize=12;
                        pHss{kstp}.Color=[ 0,0.5,0.6] ;
                        
%                        pHss{kstp}.Color=[ 0.8,0,0] ;
                    else
                        pHss{kstp}.MarkerSize=12;
                        pHss{kstp}.Color=[ 1,0,0] ;
                        
%                         pHss{kstp}.Color=[ 0.8,0,0] ;
                    end
                 end         
    end
  

    
    
  