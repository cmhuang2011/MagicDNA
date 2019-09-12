function  JsonEasier_v2
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% addpath(genpath(pwd));
%--------Last editted , July 5th 2017
%--------Last editted , May 4th 2018


[JSON_filename,PathName3,FilterIndex]= uigetfile({'*.json','JSONfile' },'Select the json file');

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

 
 dxy=2.6        ;

   NumList=[-1,-1,-1,-1];    %[k, num, col , row, GcX, GcY]
  Eff=[];
   for k=1:length(dat.vstrands)
   NumList=[NumList;   [k, dat.vstrands{k}.num,  dat.vstrands{k}.col,  dat.vstrands{k}.row ]    ];   
  end
 NumList=setdiff(NumList,[-1,-1,-1,-1],'rows');
  
  Scaf_Start=[];  %Cyl-CadIndex
  ColRowStart=[];
  ScafR_Skip=ones( 10000,1);

  for k=1:length(dat.vstrands)
    ScafInCyl= dat.vstrands{k}.scaf ;
    [~,ind1] =ismember(   ScafInCyl(:,1:2), [-1,-1],'rows') ;
    [~,ind2] =ismember(   ScafInCyl(:,3:4), [-1,-1],'rows') ;
    ind2=~ind2 ;
    All= and (ind1, ind2) ;
    if sum(All)~=0
       ClyIndx= NumList(NumList(:,1)==k,2) ;
       Base=find(All)-1;
       Scaf_Start= [ClyIndx ,Base]; 
       ColRowStart= [dat.vstrands{k}.col,dat.vstrands{k}.row]      ;   
       ScafR_Skip(1)=dat.vstrands{k}.skip(Base+1);
       break;
    end
  end
  
  ScafR_AllBase=zeros( 10000,2); ks=2; ColRow=zeros( 10000,2);
  ScafR_AllBase(1,:)=Scaf_Start;   ColRow(1,:)=ColRowStart;
  
  Current=Scaf_Start ;
  while 1
  Cyli= NumList(NumList(:,2)==Current(1),1) ;
  ScafR_AllBase(ks,: ) = dat.vstrands{Cyli}.scaf( Current(2)+1,3:4) ;
  
    ks=ks+1 ;

  ScafR_Skip(ks)= dat.vstrands{Cyli}.skip(Current(2)+1 );
  
  ColRow(ks,: ) = [dat.vstrands{Cyli}.col ,dat.vstrands{Cyli}.row] ;

  Current=dat.vstrands{Cyli}.scaf( Current(2)+1,3:4);
      if sum(Current==[-1,-1])==2
      ScafR_AllBase(sum(ScafR_AllBase,2)==0 ,:)=[] ;     
      ColRow(sum(ColRow,2)==0 ,:)=[] ;
     ScafR_Skip(ScafR_Skip==1)=[];
     
      ScafR_AllBase=ScafR_AllBase(1:end-1,:);
       ColRow=ColRow(2:end,:);
     ScafR_Skip=ScafR_Skip(2:end,:);
     
%            [ScafR_AllBase(1:20,:), ScafR_Skip(1:20) ,ColRow(1:20,:)]  
%       [ScafR_AllBase(end-5:end,:), ScafR_Skip(end-5:end) ,ColRow(end-5:end,:)]  

     
      break ;
      end
  end
  ScafR_AllBase_Mindex=ScafR_AllBase;
  [~,b2]=ismember(ScafR_AllBase(:,1) , NumList(:,2)) ;
  ScafR_AllBase_Mindex(:,1) = b2;
  
%   figure;  plot3(ColRow(:,1),ColRow(:,2),ScafR_AllBase(:,2) ,'-r')
  
  fH=figure(123);clf;
  
  plot(ScafR_AllBase_Mindex(:,2) ,ScafR_AllBase_Mindex(:,1));
  ax1=gca ; hold on ;
  %----------     
  drawnow;
  warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  jFig = get(handle(fH), 'JavaFrame'); 
   jFig.setMaximized(true);
   drawnow;
  
  
  
%   col = (1:size(ScafR_AllBase_Mindex,1))*1000;  % This is the color
%     surface([ScafR_AllBase_Mindex(:,2),ScafR_AllBase_Mindex(:,2)]',[ScafR_AllBase_Mindex(:,1),ScafR_AllBase_Mindex(:,1)]' ,[col;col],...
%         'facecol','no', 'edgecol','interp', 'linew',2)
  
%-------------------------------
   NBase=size(ScafR_AllBase_Mindex,1) ;   %including skipped bases
   BelongTransM=zeros(NBase,1);
   cc=0;
    ValidAll=NumList ;
   staplePlot=cell(  size(ColorCode,1) ,1) ;
    
  for stapSSi= 1:    size(ColorCode,1)    % fill ds part
      CylInd=ColorCode(stapSSi,1) ;
      BaseStart=ColorCode(stapSSi,2) ;
      ColorIndex= find(b==ColorCode(stapSSi,3)) ;  % also mean which transformation matrix
      Next= dat.vstrands{CylInd}.stap(BaseStart+1,3:4)  ; 
     stapAndSacf= intersect(find(ScafR_AllBase_Mindex(:,1)==CylInd),find(ScafR_AllBase_Mindex(:,2)==BaseStart));
     BelongTransM(stapAndSacf)=ColorIndex;
     cc=cc+length(stapAndSacf);
     nloop=0;
     ThisStap=zeros(200,2) ; ks=2;
     ThisStap(1,:)=[CylInd,BaseStart] ;
     
     while sum(Next==[-1,-1])~=2
        CylIndNext= ValidAll(ValidAll(:,2)==Next(1),1) ;
        BaseNext= Next(2);
        stapAndSacf= intersect(find(ScafR_AllBase_Mindex(:,1)==CylIndNext),find(ScafR_AllBase_Mindex(:,2)==BaseNext)) ;
         BelongTransM(stapAndSacf)=ColorIndex;
         Next= dat.vstrands{CylIndNext}.stap(BaseNext+1,3:4)  ; 
         nloop=nloop+1;
         
         cc=cc+length(stapAndSacf);
          ThisStap(ks,:)=[CylIndNext,BaseNext] ;   ks=ks+1 ;
           if  length(stapAndSacf)~=2
          end
     end
     ThisStap(ks:end, :)=[];
     staplePlot{stapSSi} = plot( ThisStap(:,2) , ThisStap(:,1)+0.2  ,'--');
     
     staplePlot{stapSSi}.HitTest='off' ;
%      color =ColorCode(stapSSi,3)
%      ColorCode(stapSSi,:);
%     Dec_Tri= dec2hex( hex2dec('ffffff')-color)  
%         Dec_Tri= dec2hex( color,6)  ;
        
    rgbDec=ColorCode(stapSSi,3) ;
    RGBHex=dec2hex(rgbDec,6 ) ;
    Rc=hex2dec(RGBHex(1:2)); R=Rc/256 ;
    Gc=hex2dec(RGBHex(3:4)) ;  G=Gc/256 ;
    Bc=hex2dec(RGBHex(5:6))  ; B=Bc/256 ;
%                         [R,G,B]=[Rc,Gc ,  Bc]/256 ;  
        
        
        
%     if length(Dec_Tri)== 4 
%     R= 0 ;
%     G=(256-hex2dec(Dec_Tri(end-3:end-2))) /256 ;
%     B=(256-hex2dec(Dec_Tri(end-1:end))) /256 ;
%     elseif length(Dec_Tri)== 5 
%     R= (256-hex2dec(Dec_Tri(1:end-4))) /256   ;
%     G=(256-hex2dec(Dec_Tri(end-3:end-2))) /256 ;
%     B=(256-hex2dec(Dec_Tri(end-1:end))) /256  ;
%     elseif length(Dec_Tri)== 2
%     R=0;
%     G=0;
%      B=(256-hex2dec(Dec_Tri(end-1:end))) /256 ;
%     else
%     R= (256-hex2dec(Dec_Tri(1:2))) /256 ;
%     G=(256-hex2dec(Dec_Tri(3:4))) /256 ;
%     B=(256-hex2dec(Dec_Tri(5:6))) /256 ;
%     end
%     R= (256-hex2dec(Dec_Tri(1:2))) /256 
%     G=(256-hex2dec(Dec_Tri(3:4))) /256 
%     B=(256-hex2dec(Dec_Tri(5:6))) /256 

    
    staplePlot{stapSSi}.Color=[R,G,B];
     staplePlot{stapSSi}.UserData.BelongM= ColorIndex*ones(size(staplePlot{stapSSi}.XData)) ;
     
%      ColorCode(stapSSi,3)
     staplePlot{stapSSi}.UserData.IniPosition=ThisStap;
      staplePlot{stapSSi}.UserData.CLRCode=  ColorCode(stapSSi,3) ;
     nloop;
  end   
%--------------------------
  set(ax1,'Ydir','reverse'); ax1.Position= [ 0.08, 0.6 , 0.8 , 0.35] ;
  ax1.XLim=ax1.XLim + [ -20, 20] ;  ax1.YLim=ax1.YLim + [ -5, 5] ;




   ax2=axes(); 
   ax2.Parent=fH; ax2.Position=[ 0.08, 0.1 , 0.8 , 0.35 ]; 
   pH=plot(ScafR_AllBase_Mindex(:,2) ,ScafR_AllBase_Mindex(:,1));  ax2.XLim=ax2.XLim + [ -20, 20] ;  ax2.YLim=ax2.YLim + [ -5, 5] ;
   hold on;set(ax2,'Ydir','reverse');


   
% h = surf(peaks);
% colormap hsv
% fig = figure;
% ax = axes;
% new_handle = copyobj(h,ax);
% colormap(fig,hsv)
% view(ax,3)
% grid(ax,'on')
for k=1:length(staplePlot)
   staplePlotAx2{k} = copyobj(staplePlot{k},ax2) ;
end
%     h = findobj(gcf,'type','axes')

%-------------
   RemainNotAssign=find(BelongTransM==0);   LRNA2=length(RemainNotAssign) ;
   ScafNoAss=  RemainNotAssign;
  UnL=length(ScafNoAss) ;
   nwh=1;
  while ~isempty(ScafNoAss) 
       for kk=1:length(ScafNoAss)
             CylBase=   ScafR_AllBase_Mindex(ScafNoAss(kk),:) ;
             CADNANOINDEX=  ValidAll( ValidAll(:,1)==CylBase(1),2);
             TwoSideConnec = dat.vstrands{CylBase(1)}.scaf(CylBase(2)+1,:)  ;%read cadnano scaf 
             [~,indLeft]= ismember( TwoSideConnec(1:2),ScafR_AllBase,'rows') ;
             if TwoSideConnec(1)~=CADNANOINDEX   || abs(CylBase(2)-TwoSideConnec(2))~=1
                 indLeft=0;   % left side cross cylinder or jump
             end
             
             [~,indRight]= ismember( TwoSideConnec(3:4),ScafR_AllBase,'rows');
             if TwoSideConnec(3)~=CADNANOINDEX   || abs(CylBase(2)-TwoSideConnec(4))~=1
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

                    Gindex= intersect(find(ScafR_AllBase_Mindex(:,1)==CylFirstRep),find(ScafR_AllBase_Mindex(:,2)==Trust(2)));

                    if  BelongTransM(Gindex)~=0

                      BelongTransM( ScafNoAss(kk))= BelongTransM(Gindex)   ;
                    else
                        nothinghappen=1 ;
                    end
            end
       end
      [nwh length(ScafNoAss) ];
       RemainNotAssign=find(BelongTransM==0);
       ScafNoAss=  RemainNotAssign ;
       nwh=nwh+1 ;
       if nwh>80
          break 
       end
  end
   [nwh length(ScafNoAss) ]
%----------------------------
     fH.UserData.ColRow=ColRow;
     fH.UserData.ScafR_Skip=ScafR_Skip;
         Str={'--'};  for k=1:max(BelongTransM) ;   Str{k}=strcat('Bun  ',num2str(k));   end
        %---------
         checkH_View = uicontrol('Style', 'checkbox','String', 'Axis auto/equal','Unit','normalized','Position', [0.9 0.55 0.05 0.05]); 
         checkH_View.Callback=@(src,evn)checkFcn2(src,evn,ax2)  ; 

        popupH = uicontrol('Style', 'popup','String', Str,'Unit','normalized','Position', [0.9 0.87 0.06 0.08]); 
        editH = uicontrol('Style', 'edit','String', '1','Unit','normalized','Position', [0.9 0.85 0.05 0.05]); 
        
        pH.ButtonDownFcn=@(src,evn)lineselect(src,evn,popupH,pH , staplePlotAx2, fH ,BelongTransM,ax2);

        set(fH,'KeyPressFcn',@(src,evn)keyMove(src,evn,fH,pH , staplePlotAx2,BelongTransM ,popupH  ,editH,ax2  )  ) ;
        popupH.Callback=@(src,evn) popupFcn(src,evn,pH , staplePlotAx2, fH ,BelongTransM,ax2  );
               
        xlabel('X-axis') ;ylabel('Y-axis')
        grid on;set(gca,'XMinorTick','on');
        
        btn1 = uicontrol('Style', 'pushbutton', 'String', 'Align','Unit','normalized', 'Position', [0.9 0.3 0.05 0.05] ,...
            'Callback', {@(src,evn)AlignRouting(src,evn,fH,pH , staplePlotAx2,BelongTransM ,popupH  ,editH,ax2,ScafR_AllBase_Mindex ,NumList )});       
        
        btn2 = uicontrol('Style', 'pushbutton', 'String', 'ExportNewJson','Unit','normalized', 'Position', [0.9 0.2 0.05 0.05] ,...
            'Callback', {@(src,evn)ExportJSON(src,evn,fH,pH , staplePlotAx2,BelongTransM ,popupH  ,editH,ax2,ScafR_AllBase_Mindex ,NumList,JSON_filename,PathName3 )});       
           
       fprintf(' Start moving \n ')    

%%  
    function ExportJSON(src,evn,fH,pH , staplePlotAx2,BelongTransM ,popupH  ,editH,ax2,ScafR_AllBase_Mindex ,NumList,JSON_filename,PathName3 )
        ColRowLocal=fH.UserData.ColRow ;
        CylsOld=unique(ScafR_AllBase_Mindex(:,1));
        
        XYAll=[pH.XData ;pH.YData]' ;
        CylsNew=unique(XYAll(:,2));
        
        New_NumList=zeros(size(CylsNew,1),4) ;
        New_NumList(:,1)=CylsNew;
        
        EEven=0;  OOdd=1;
        for cylN=1:size(New_NumList,1)
            Bases=XYAll(:,2)==New_NumList(cylN,1) ;
            BasesY= XYAll( Bases ,1) ;
            Neibors= find(abs(BasesY(1:end-1)-BasesY(2:end))==1) ;
            QQ=BasesY(Neibors+1)- BasesY( Neibors) ;         
             if sum(QQ==-1)~= length(QQ)
                  New_NumList(cylN,2)= EEven ; EEven=EEven+2;
             elseif sum(QQ==1)~= length(QQ)
                  New_NumList(cylN,2)= OOdd ; OOdd=OOdd+2;    
             else
                 sdfsf=234;
             end
             
             if size(unique(fH.UserData.ColRow(Bases,:) ,'rows') ,1) ==1
                 sdfsf=3;
                 New_NumList(cylN,3:4)=unique(fH.UserData.ColRow(Bases,:),'rows' );
             else
                 sfdfsf=3
             end
             
        end
        %----------
        fH.UserData.ColRow;
         fH.UserData.CylinderOrder;
        
         ScafBase_k_numCyl_Z=zeros(size(fH.UserData.ColRow));  % [showingOrder, numCyl, Z];
         
         ScafBase_k_numCyl_Z(:,3) =pH.XData'  ;
         for ii=1: size(ScafBase_k_numCyl_Z,1) 
            ColRow_ii=fH.UserData.ColRow(ii,:)  ;
            [~,num_ii_Ind ]= ismember(   ColRow_ii , New_NumList(:,3:4) ,'rows') ;
            ColRow_ii;
            num_ii_Ind;
            num_ii=   New_NumList(num_ii_Ind,2) ;
            ScafBase_k_numCyl_Z(ii,2)=num_ii;       
            ScafBase_k_numCyl_Z(ii,1)=num_ii_Ind;       
         end
         
         %---------
         Stap_Base_k_numCyl_Z =cell( length(staplePlotAx2), 1) ;
         staplePlotAx2;
         for sstap= 1 :length(staplePlotAx2)
%              sstap
         [~, QQ] =ismember(staplePlotAx2{sstap}.UserData.IniPosition,  ScafR_AllBase_Mindex ,'rows') ;
%          if sstap==22
%              sdfsf=3
%          end
         QQ(QQ==0) =[];     % fix overhang case , Oct 2,2018
         Stap_Base_k_numCyl_Z{sstap} = ScafBase_k_numCyl_Z(QQ ,: ) ;
         end
         
         %-----------
           MM_zzzz= max(ScafBase_k_numCyl_Z(:,3));
           MM_zzzz= floor(MM_zzzz/32)*32 +64 ;

           
            datT=loadjson('CadNano.json');
            NNdat=datT;
%             jj=randi([10,99],1,1);
%             JJ=strcat('TBund',int2str(jj));
%             JJ=strcat(JJ,'.json');
            
%             prompt = {'Enter File name:'};
%             dlg_title = 'Input';
%             num_lines = 1;
%             defaultans = {'No'};
%             answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
%             if strcmp(answer,'No')
%             return;
%             end

%  New_NumList;
% ScafBase_k_numCyl_Z;   % [cylInde,  num , Z];
% Scaf_Start;
% dat.vstrands{20}.scaf(179,:);
% dat.vstrands{20}.scaf(177:181,:);

            file_name=strcat( 'Testttt2','.json');         
         JSON_filename
       
                     file_name=strcat(   JSON_filename(1:end-5),'_moved '   ,'.json');         

         
            NNdat.name=file_name;
            vHelix=NNdat.vstrands{1} ;
            vHelix.stap_colors=[-999,-888];
            vHelix.stap=-1 * ones(MM_zzzz,4);
            vHelix.skip=zeros(1,MM_zzzz);
            vHelix.scaf=-1 * ones(MM_zzzz,4);
            vHelix.loop=zeros(1,MM_zzzz);
            
            for cyli=1:size(New_NumList,1)
             NNdat.vstrands{cyli}=vHelix;
             NNdat.vstrands{cyli}.num = New_NumList(cyli,2) ;
             NNdat.vstrands{cyli}.col = New_NumList(cyli,3) ;
             NNdat.vstrands{cyli}.row = New_NumList(cyli,4) ;
            end
   
  %------------------------------------------------------------------------------------------------scaf           
   
            CurrentPP= ScafBase_k_numCyl_Z(1,:);
            NextPP=ScafBase_k_numCyl_Z(2,:);
            NNdat.vstrands{CurrentPP(1)}.scaf( CurrentPP(3)+1,:) =[-1,-1, NextPP(2), NextPP(3)] ;
            NNdat.vstrands{CurrentPP(1)}.skip(CurrentPP(3)+1)=fH.UserData.ScafR_Skip(1) ;
            
            for bbi=2: size(ScafBase_k_numCyl_Z,1)-1
             NNdat.vstrands{NextPP(1)}.skip( NextPP(3)+1) =fH.UserData.ScafR_Skip(bbi) ;
   
                
             FuturePP=    ScafBase_k_numCyl_Z(bbi+1,:) ;
             NNdat.vstrands{NextPP(1)}.scaf( NextPP(3)+1,:) =[CurrentPP(2), CurrentPP(3), FuturePP(2), FuturePP(3)] ;   
             CurrentPP=  NextPP;
             NextPP= FuturePP;
            end
             NNdat.vstrands{NextPP(1)}.skip( NextPP(3)+1) =fH.UserData.ScafR_Skip(bbi+1) ;
            NNdat.vstrands{NextPP(1)}.scaf( NextPP(3)+1,:)=[ CurrentPP(2), CurrentPP(3),-1,-1];
            
  %------------------------------------------------------------------------------------------------stap   
       for stpi=1:length(staplePlotAx2)
            CurrentPP= Stap_Base_k_numCyl_Z{stpi}(1,:);
            NextPP=Stap_Base_k_numCyl_Z{stpi}(2,:);
            NNdat.vstrands{CurrentPP(1)}.stap( CurrentPP(3)+1,:) =[-1,-1, NextPP(2), NextPP(3)] ;
            
            NNdat.vstrands{CurrentPP(1)}.stap_colors=[ NNdat.vstrands{CurrentPP(1)}.stap_colors ; CurrentPP(3) ,  staplePlotAx2{stpi}.UserData.CLRCode] ; 
            
            for bbi=2: size(Stap_Base_k_numCyl_Z{stpi},1)-1
             FuturePP=    Stap_Base_k_numCyl_Z{stpi}(bbi+1,:) ;
             NNdat.vstrands{NextPP(1)}.stap( NextPP(3)+1,:) =[CurrentPP(2), CurrentPP(3), FuturePP(2), FuturePP(3)] ;   
             CurrentPP=  NextPP;
             NextPP= FuturePP;
            end
            NNdat.vstrands{NextPP(1)}.stap( NextPP(3)+1,:)=[ CurrentPP(2), CurrentPP(3),-1,-1];
            
       end    
            
 %------------------------------------------------------------------------------------------------           
         for ColorCheck=1:      length(NNdat.vstrands)
        %      NNdat.vstrands{ColorCheck}
            if  size(NNdat.vstrands{ColorCheck}.stap_colors ,1 )>1
                 NNdat.vstrands{ColorCheck}.stap_colors=   circshift(sortrows(NNdat.vstrands{ColorCheck}.stap_colors),-1) ;
            end
         end
 
%          for k3=1:length( NNdat.vstrands)
%         NNdat.vstrands{k3}.skip= zeros( size( NNdat.vstrands{k3}.skip));
%          end        
         
         
             TTtext=savejson('Title',NNdat,'ArrayIndent',0,'Compact',1 );
%             TTtext=savejson('Title',NNdat,'ArrayIndent',0,'Compact',1 );
            TTtext(1:10)=[];
            TTtext(end-1:end)=[];
%             TTtext(TTtext==' ')=[];
%             TTtext(TTtext=='	')=[];
            IOfSC2=strfind(TTtext, ',[-999,-888]');  %for color json export
%             Cop=IOfSC2;
            for removedd=1:length(IOfSC2)
                 UUdataPosittion=strfind(TTtext, ',[-999,-888]');  %for color json export
                 UUdataPosittion;
               Exxtraindex= UUdataPosittion(1);
                 TTtext(Exxtraindex:Exxtraindex+11)=[];
%                  Cop=strfind(TTtext, ',[-999,-888]');
            end            
            
            SecTerm999888=strfind(TTtext, '-999,-888');  %for color json export
            for removedd=1:length(SecTerm999888)
                 UUdataPosittion=strfind(TTtext, '-999,-888');  %for color json export
                 if ~isempty(UUdataPosittion)
                 Exxtraindex= UUdataPosittion(1);
                 end
                  TTtext(Exxtraindex:Exxtraindex+8)=[];
            end
            
           JsonFolder=[pwd '\'];
%            PathName3
           
%             fileID = fopen(file_name,'w');
%              fileID = fopen([JsonFolder file_name],'w');
              fileID = fopen([PathName3 file_name],'w');
            fprintf(fileID,TTtext);
            fclose(fileID);
%             ReadTest=loadjson(file_name)      
%             ReadTest=loadjson([JsonFolder file_name])                 
            ReadTest=loadjson([PathName3 file_name])                 
            
            
            
         sdgfg=234;
         
%         G_CylTable=[-1,-1];
%         for kgi=1:length(BelongTransM)
%         G_CylTable=union(G_CylTable, [ScafR_AllBase_Mindex(kgi,1) ,XYAll(kgi,2)  ]   ,'rows') ;
%         end
%         G_CylTable=setdiff(G_CylTable, [-1,-1],'rows') ;
%         %------^^^^^^^^^
%        
%         InplaneM= max(NumList(:,3:4))+ [1,1] ;
%         BW=zeros(InplaneM) ;
%         for ka=1:size(NumList,1)
%         BW(NumList(ka,3)+1 ,NumList(ka,4)+1 )=1;
%         end
%         L = bwlabel(BW); 
%         
%         AddNumList=NumList;
%         for kb=1:size(NumList,1)
%          AddNumList(kb,5)= L(    AddNumList(kb,3)+1, AddNumList(kb,4)+1);
%         end
%         
%         RealBundles=unique(AddNumList(G_CylTable(:,1),5));   % has sacf
% %         Col_Row_RandMove= zeros( length(RealBundles),2) ;
%         
%         OneDBW=zeros(1, max(G_CylTable(:,2))) 
%         OneDBW(G_CylTable(:,2))=1 ;
%         L_one = bwlabel(OneDBW); 
%         
%         InplaneBuns=zeros(size(G_CylTable));  %[ Old-2Dbundle,  New-1DBundle]
%          for kv=1:size(InplaneBuns,1)
%           InplaneBuns(kv,1)= AddNumList(G_CylTable(kv,1) ,5);
%           InplaneBuns(kv,2)= L_one(G_CylTable(kv,2) );
%          end
%          InplaneBuns=unique(InplaneBuns,'rows');
%          leftOldRef=zeros( length(unique(InplaneBuns(:,1))), 3); %-----[OldBund, RefCol, RefRow];
%          for Ri=1:size(leftOldRef,1)
%              Old2DBund=InplaneBuns(Ri,1);
%              AllColRowInthis=AddNumList( AddNumList(:,5)==Old2DBund,3:4) ;
%              AllColRowInthis=sortrows(AllColRowInthis);
%              leftOldRef(Ri,:)=[InplaneBuns(Ri,1)  AllColRowInthis(1,:)];
%          end
%          
%          RightNewRef=zeros(length(unique(InplaneBuns(:,2))), 3); %-----[NewBund, RefCol, RefRow];
%          
%          for Rj=1:size(RightNewRef,1)
%             New1DbundInd= find(InplaneBuns(:,2)==Rj); New1DbundInd=New1DbundInd(1);
%             RightNewRef(Rj,:)=  [Rj, leftOldRef( New1DbundInd,2:3)] ;
%          end
%          RightNewRef=sortrows(RightNewRef);
%          
%          New_NumList;  G_CylTable;
%          L_one;AddNumList;  InplaneBuns;  RightNewRef;   leftOldRef;
%          
%          
%          New_NumList  ;
%          used=[];
%          for rri=1:size(InplaneBuns,1)
%              if ~ismember(InplaneBuns(rri,2),used)
%                  iiind=InplaneBuns(rri,1);
%                  sss=(leftOldRef(:,1)==iiind);
%                  [~,a]=ismember( leftOldRef(sss,2:3), AddNumList(:,3:4) ,'rows') ;
%                  Og= AddNumList(a,1) ;
%                  Ng= G_CylTable(G_CylTable(:,1)==Og ,2) 
%                  if isempty(Ng)
%                      sdf=3
%                  end
%                  
%                  IndInNew_NumList=(New_NumList(:,1)==Ng);
%                  New_NumList(IndInNew_NumList,3:4)= RightNewRef(InplaneBuns(rri ,2) ,2:3) ;
%                  used=union(used ,InplaneBuns(rri,2) );
%              end
%          end
%          sdff=3
%          for c=1: size(New_NumList,1)
%              Gindex_new= New_NumList(c,1) ;
%              OldGCylinders= G_CylTable((G_CylTable(:,2)==Gindex_new),1)  ;
%              
%              [~,Iid ] =ismember( OldGCylinders , G_CylTable(:,1)) ;
%              NewGCyl= G_CylTable(Iid,2) ;
%              %maybe multiple;
%              
%              RefBunOld=AddNumList(OldGCylinders,5);
%              
% %              InplaneBuns
% %              RefBunOld
%              [~,Indxx]=ismember( RefBunOld , leftOldRef(:,1)) ;
%              RefXY=leftOldRef(Indxx, 2:3) ;
%              
%              
%              OldColRow=AddNumList(OldGCylinders,3:4) ;
%              
%              [RefXYExist,  IndIn_NN] = ismember( RefXY , New_NumList(:,3:4) ,'rows') ;
%              
% %              Newd_ColRo=
%              dColRow= OldColRow(RefXYExist ,:) -RefXY(RefXYExist ,:) ;
%              if isempty(  IndIn_NN(RefXYExist))
% %                  continue
%              RefBunNew= InplaneBuns(InplaneBuns(:,1)==RefBunOld,2);
%               RefBunNew_XY= RightNewRef(RightNewRef(:,1)==RefBunNew,2:3) ;
%                 [RefXY_Inplane,  Ind22] = ismember( RefXY , leftOldRef(:,2:3) ,'rows') ;  
%                  Number=  leftOldRef(Ind22,1) ;
%                  
% % %                  
% %                  InplaneBuns
% %                  
%                   esf=3;
%                  New_NumList(c,3:4)=  OldColRow-RefXY +RefBunNew_XY;
%                  continue;
%              end
%              IndIn_NN(RefXYExist)
%              New_NumList(c,3:4)= dColRow+ New_NumList( IndIn_NN(RefXYExist),3:4) ;
%          end
% 
%           New_NumList;  G_CylTable;
%          L_one;AddNumList;  InplaneBuns;  RightNewRef;   leftOldRef;
% 
% 
%          
%         sdfsf=234
        
        
    end

%     function [c,ceq] = ContactConstraint(x,Centers,Radius,RequiredEdge,IntialPCell,EdgeNotContact)
%         ceq=[];
%         c=zeros(size(RequiredEdge,1)+size(EdgeNotContact,1)   ,1) ;
%         a=5;
%         for k_r=1:  size(RequiredEdge,1)
%             Bi=RequiredEdge(k_r,1) ;         Bj=RequiredEdge(k_r,2) ;
%         c(k_r)  = norm(Centers(Bi,:)+x(Bi,:) - Centers(Bj,:)-x(Bj,:) )- Radius(Bi)- Radius(Bj) ;
%         end
%         
%         for k_r2=1:size(EdgeNotContact,1)
%             Bi=EdgeNotContact(k_r2,1) ;         Bj=EdgeNotContact(k_r2,2) ;
%         c(k_r+k_r2)  = -(  norm(Centers(Bi,:)+x(Bi,:) - Centers(Bj,:)-x(Bj,:) ) -a* Radius(Bi)-a* Radius(Bj) ) ;
%         end        
%              
%     end

    function [c,ceq] = ContactConstraint(x,Centers,Radius,RequiredEdge,IntialPCell,EdgeNotContact)
        ceq=[];
        c=zeros(size(EdgeNotContact,1)   ,1) ;
        a=5;
%         for k_r=1:  size(RequiredEdge,1)
%             Bi=RequiredEdge(k_r,1) ;         Bj=RequiredEdge(k_r,2) ;
%         c(k_r)  = norm(Centers(Bi,:)+x(Bi,:) - Centers(Bj,:)-x(Bj,:) )- Radius(Bi)- Radius(Bj) ;
%         end
        
        for k_r2=1:size(EdgeNotContact,1)
            Bi=EdgeNotContact(k_r2,1) ;         Bj=EdgeNotContact(k_r2,2) ;
        NewBi=IntialPCell{Bi} +  ones(size(IntialPCell{Bi},1),1)*x(Bi,:) ;  
              NewBi(:,1)= NewBi(:,1) -min(NewBi(:,1))+2 ;  NewBi(:,2)= NewBi(:,2) -min(NewBi(:,2))+2 ;   %add 5/21
        NewBj=IntialPCell{Bj} +  ones(size(IntialPCell{Bj},1),1)*x(Bj,:) ;
         NewBi=round(NewBi) ;
         NewBj=round(NewBj) ;
        BW=zeros(max(   [NewBi;  NewBj])) ;
         BW(  NewBi(:,1), NewBi(:,2))=1   ;  
         BW(  NewBj(:,1), NewBj(:,2))=1    ;
          [~,num] = bwlabel(BW);             
        c(k_r2)  = 1.5 - num ;
        
        end        
             
    end




    function f=InPlaneMoveFunc(x ,Edges, IntialPCell,RequiredEdge)
        f=0;
        x;
        for Ei=1:size(RequiredEdge,1)
        Bi=RequiredEdge(Ei,1) ;
        Bj=RequiredEdge(Ei,2) ;
       
        
        NewBi=IntialPCell{Bi} +  ones(size(IntialPCell{Bi},1),1)*x(Bi,:) ;
        NewBj=IntialPCell{Bj} +  ones(size(IntialPCell{Bj},1),1)*x(Bj,:) ;
         NewBi=round(NewBi) ;
         NewBj=round(NewBj) ;
        Val= size( intersect(NewBi,NewBj,'rows') ,1)  ;
          
        %----
%           q=0;
%            for cell_i=1: size(IntialPCell{Bi},1)
%                e=Inf;
%             for cell_j=1: size(IntialPCell{Bj},1)
% %                 Den=size(IntialPCell{Bi},1)*size(IntialPCell{Bj},1);
%                col_E=   abs( (IntialPCell{Bi}(cell_i,1)+ x(Bi,1)) - ( IntialPCell{Bj}(cell_j,1)+ x(Bj,1) ) ) ;
%                row_E=   abs( (IntialPCell{Bi}(cell_i,2)+ x(Bi,2)) - ( IntialPCell{Bj}(cell_j,2)+ x(Bj,2) ) ) ;
% %                 f= f+ abs(col_E) + abs( row_E) ;
%                 e=min( e, abs(col_E-q)  +  abs(row_E-q)  );
%             end
%             
%               f=f - e;
%            end
           %----
           f= f- Val;
           
        end   
%         f= f  ;
            
%         sdfsf=3;
    end

% function f=InPlaneMoveFunc(x ,Edges, IntialPCell)
% 
%         a1=1;
%         a2=1;
%         f=0;
%         for Ei=1:size(Edges,1)
%         Bi=Edges(Ei,1) ;         Bj=Edges(Ei,2) ;
%         e=Inf;  q=2;
%            for cell_i=1: size(IntialPCell{Bi},1)
%             for cell_j=1: size(IntialPCell{Bj},1)
% %                 Den=size(IntialPCell{Bi},1)*size(IntialPCell{Bj},1);
%                col_E=   abs( (IntialPCell{Bi}(cell_i,1)+ x(Bi,1)) - ( IntialPCell{Bj}(cell_j,1)+ x(Bj,1) ) ) ;
%                row_E=   abs( (IntialPCell{Bi}(cell_i,2)+ x(Bi,2)) - ( IntialPCell{Bj}(cell_j,2)+ x(Bj,2) ) ) ;
% %                 f= f+ abs(col_E) + abs( row_E) ;
%                 e=min( e, abs(col_E-q)  +  abs(row_E-q)  );
%             end
%            end
%            f=f+e;
%         end   
%             
% %         sdfsf=3;
%     end

    function  Out=CM_round(inputs)
       Out=zeros(size(inputs));
       for jj=1:   size(inputs,1)
          P= inputs(jj,:) ;
          x =round(P(1))-2:1:round(P(1))+2;
          y = round(P(2))-2:1:round(P(2))+2;
         [X,Y] = meshgrid(x,y) ;
         xx=reshape(X, [25,1]) ;         yy=reshape(Y, [25,1]);
         Candi=[xx,yy];
         
         Candi(mod(sum(Candi,2),2)==1,:)=[];
         
         Dxy=Candi- ones(size(Candi,1),1)*P ;
         d= Dxy(:,1).^2 + Dxy(:,2).^2;
%          d
%          Candi
         Ind=find(d==min(d)); 
         Ind=Ind(1);
        First= Candi(Ind,:);
        DD=First-P ;
        
        Out(jj,:)=First ;
       end       
%        Out= inputs + ones(size(inputs,1),1) *DD;
  
        
    end

    function AlignRouting(src,evn,fH,pH , staplePlotAx2,BelongTransM ,popupH  ,editH,ax2 ,ScafR_AllBase_Mindex ,NumList )  
        ColRowLocal=fH.UserData.ColRow ;
%         CurrentAllXY=[pH.XData',  pH.YData'];
%         BundlesMove=zeros(max(BelongTransM) ,2) ;
%         for kA=1:max(BelongTransM)
%             Ind=BelongTransM==kA;
%             CXY=    CurrentAllXY(Ind,:) ;
%             O_XY= flip(ScafR_AllBase_Mindex(Ind,:) ,2);
%             Movement= unique(CXY-O_XY,'rows') ;
%             BundlesMove(kA,:)=Movement;
%         end
        %-------------
        nBundle=max(BelongTransM) ;
%         unique(ColRow(BelongTransM==2,:),'rows')
        BundlesPlanarAdjM=zeros(nBundle,nBundle);
       for i=1:nBundle
           for j=i+1:nBundle
              Yplot_i =unique(pH.YData(BelongTransM==i));
              Yplot_j =unique(pH.YData(BelongTransM==j));
             if ~isempty( intersect(Yplot_i,Yplot_j ) )
              BundlesPlanarAdjM(j,i)= 1; BundlesPlanarAdjM(i,j)= 1;
             end
           end
       end
        BundlesPlanarAdjM
        PlanarConnGraph=graph(BundlesPlanarAdjM);
        Group= conncomp(PlanarConnGraph) ;   N_Group=max(Group);
        
        BundleInit_ColRow=cell(nBundle,1); centers=zeros(nBundle,2); Radius=zeros(nBundle,1);
        for kq=1:nBundle   
            BundleInit_ColRow{kq}= unique(ColRowLocal(BelongTransM==kq,:),'rows') ;
           centers(kq,:)= mean(  BundleInit_ColRow{kq}) ;
           Radius(kq)= min(max(BundleInit_ColRow{kq})-min(BundleInit_ColRow{kq}));
        end
        
        x_opt=zeros(nBundle,2);
        lb=zeros(nBundle,2); ub= 30*ones(nBundle,2);
        for kw=1:nBundle
        lb(kw,:)= -min(BundleInit_ColRow{kw}) + [2 , 2];
        end
        
        ub=lb+30 ;
        
%         x_opt=lb;
%         
%         x_opt=20*rand(size(x_opt))
%          x = fmincon(fun,x0,A,b,[],[],lb,ub)
%         x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
        [u,v]=find(BundlesPlanarAdjM ==0) ;  Ind= u<v ;
        u=u(Ind) ; v=v(Ind); EdgeNotContact=[u,v] ;
        [u,v]=find(BundlesPlanarAdjM ==1) ;  Ind= u<v ;
        u=u(Ind) ; v=v(Ind); RequiredEdge=[u,v] ;
       
        
%         fun=@(x)InPlaneMoveFunc(x ,EdgeNotContact, BundleInit_ColRow ,RequiredEdge  ) ;
% %         fun(x_opt ,EdgeNotContact, BundleInit_ColRow)
%        options = optimoptions('fmincon','Algorithm','interior-point','MaxFunctionEvaluations',3000  ,'Display','iter' ,'DiffMinChange',1 );
% %        options=[];
% %         [x,fval]= fmincon(fun,x0,Aine,bine,Aeq,beq,lb,ub,nonlcon,options) ;
%        nonlcon=@(x)ContactConstraint(x,centers,Radius,RequiredEdge,BundleInit_ColRow,EdgeNotContact) ;
%          [x,fval]= fmincon(fun,x_opt,[],[],[],[],lb,ub,nonlcon, options);
         
         nloopW=1; nMax=200;
         BestV=Inf; BestNx=Inf;
         Bestx=[]; 
         
%          x=lb;
         TargetBundleCenter= [6:6:14+N_Group*6  ;  6:6:14+N_Group*6  ]' ;
        

         x=zeros(size(lb));
         for kzz=1:size(x,1)
           x(kzz,:)=    TargetBundleCenter(Group(kzz),:) - centers(kzz,:) ;
         end
%           for kcc= 1: size(lb,1)
%            x(kcc,1)= randi([lb(kcc,1), ub(kcc,1) ]) ;
%            x(kcc,2)= randi([lb(kcc,2), ub(kcc,2) ]) ;
%          end
           dx=zeros(size(lb));
           dx2=zeros(size(lb));
           x
           
           List=[-2,0; 2,0 ; 0,0; 0, 2 ;0 , -2 ; 1,-1 ; 1,1;-1,-1;-1,1] ;
         while nloopW<nMax

          if     nloopW<10
              dx=zeros(size(x));
          else
              dx=  List(randi(size(List,1)  ,size(lb,1),1),:)  ;
          end
    %              dx2=zeros(size(lb));
    %              for dF=1:size(            RequiredEdge,1)
    %                  OriCenter =centers(RequiredEdge(dF,:) ,:) ;
    %                  Mov=  x(RequiredEdge(dF,:) ,:) + dx(RequiredEdge(dF,:) ,:)   ;
    %                  CurrP=OriCenter+Mov;
    %                  Vec= CurrP(2,:) - CurrP(1,:); 
    %                  dx2(RequiredEdge(dF,1),:)= dx2(RequiredEdge(dF,1),:) + round(0.5*Vec)  ;
    %                  dx2(RequiredEdge(dF,2),:)= dx2(RequiredEdge(dF,2),:) - round(0.5*Vec)  ;
    %              end
    %             dx=dx

              val=   InPlaneMoveFunc(CM_round(x+dx) ,EdgeNotContact, BundleInit_ColRow ,RequiredEdge  ) ;
              nonC=ContactConstraint(x+dx,centers,Radius,RequiredEdge,BundleInit_ColRow,EdgeNotContact)   ;
              if sum(nonC>0)==0
                if val<BestV
    %                 if  norm(x)<= BestNx
                  Bestx=x+dx ;
                  BestV=  val;
                  x=CM_round(x+dx );
                  [nloopW ,BestV]
    %               x=Bestx- dx ;
    %               BestNx=norm(x);
    %                 end
                end 
              end

             nloopW=nloopW+1    ;
         end   % end of while 
         
         val=   InPlaneMoveFunc(Bestx ,EdgeNotContact, BundleInit_ColRow ,RequiredEdge  )
         ContactConstraint(Bestx,centers,Radius,RequiredEdge,BundleInit_ColRow,EdgeNotContact) 
         x=Bestx ;
         BestV
         
         ColRowMove=CM_round(x)
%        ColRowMove=  [round(x(:,1)/2)*2,  round(sum(x,2)/2)*2- round(x(:,1)/2)*2]
         
%          ColRowMove=x;

         figure(108);clf;subplot(1,2,1)
         hold on;set(gca,'Ydir','reverse');
         for ke=1:length(BundleInit_ColRow)
         scatter3(BundleInit_ColRow{ke}(:,1) ,BundleInit_ColRow{ke}(:,2) ,ke*ones(size(BundleInit_ColRow{ke},1),1)  ) ;
         XY=mean(BundleInit_ColRow{ke} );
         text(XY(1),XY(2) ,ke, num2str(ke) ) 
         end
         subplot(1,2,2); hold on; set(gca,'Ydir','reverse');
         for kr=1:length(BundleInit_ColRow)
         scatter3(BundleInit_ColRow{kr}(:,1)+ColRowMove(kr,1) ,BundleInit_ColRow{kr}(:,2)+ColRowMove(kr,2),kr*ones(size(BundleInit_ColRow{kr},1),1)  )
          XY=mean(BundleInit_ColRow{kr});
         text(XY(1)+ColRowMove(kr,1),XY(2)+ColRowMove(kr,2) ,kr, num2str(kr) ) 
         end
        
%          x2=zeros(size(x)); x2(5,:)=[2,2];
         
         New_BundleInit_ColRow=BundleInit_ColRow;
         
         New_ColRaw=[-1,-1,-1,-1,-1];
         for kt=1:length(BundleInit_ColRow)
         New_BundleInit_ColRow{kt}= New_BundleInit_ColRow{kt} +ColRowMove(kt,:) ;
         
         New_ColRaw= union(New_ColRaw, [  kt*ones(size(BundleInit_ColRow{kt}(:,1))) ,BundleInit_ColRow{kt}, New_BundleInit_ColRow{kt} ],'rows' ) ;
         end
         New_ColRaw=setdiff(New_ColRaw,[-1,-1,-1,-1,-1],'rows');
         
%          
         ColRowLocal;
         Old_Y= pH.YData' ; NewY=zeros(size(Old_Y));
         NewColRawByBase=   zeros(size(ColRowLocal)) ;
         for basei=1:size(ColRowLocal,1) 
             OldColRowi= ColRowLocal(basei,:); 
%              [~ , inddd] = ismember(OldColRowi, New_ColRaw(:,1:2), 'rows') ;
             inddd= intersect(find(OldColRowi(1)==New_ColRaw(:,2)) ,find(OldColRowi(2)==New_ColRaw(:,3)));
             if length(inddd)>1
                 BelongBundles= New_ColRaw(inddd,1) ;
%                  if intersect(BelongBundles, BelongTransM(basei) )  > length(inddd)
%                      dfsfg=21
%                  end
%                  inddd
                 
                 inddd= inddd(  find(BelongBundles==BelongTransM(basei) )   );
             end
             NewColRawByBase(basei,:)= New_ColRaw(inddd, 4:5);
%              [~,IndY]=ismember(NewColRawByBase(basei,:),New_ColRaw(:, 4:5), 'rows') ;
%              NewY(basei)=IndY ;
             New_ColRaw(inddd, 6) = Old_Y(basei);
         end
         
         GroupHeight=zeros(max(Group),1);
         for k12=1:length(GroupHeight)
            Budnles=  find(Group==k12) ;
             ys=0; ns=0;
             for jj=1:length(Budnles)
             ys= ys + sum(Old_Y(BelongTransM==Budnles(jj))) ;
                 ns=ns + sum(BelongTransM==Budnles(jj)) ;
             end
             GroupHeight(k12) = ys/ns ;
         end
         [~,bGH]=sort(GroupHeight) ;
         
         
         NewYOrderList=New_ColRaw(:, [1 4 5]); 
         for kk4=1:size(NewYOrderList,1)
         NewYOrderList(kk4,4)=  find(Group(NewYOrderList(kk4,1)) ==bGH)  ;
         end
         
         NewYOrderList=unique(NewYOrderList(:,2:4),'rows') ;
         NewYOrderList=sortrows(NewYOrderList,[3 2 1])  ;
         
%          New_ColRaw= sortrows(New_ColRaw,[1 4 5])  ;   % respect to user assign Y
         for  basei=1:size(ColRowLocal,1) 
          [~,IndY]=ismember(  NewColRawByBase(basei,:),NewYOrderList(:,1:2) , 'rows') ;
          NewY(basei)=IndY ;
         end
         
         BundlesPlanarAdjM
         Group;

         
         if isfield(popupH.UserData,'HighLight' )
         delete(popupH.UserData.HighLight);
         popupH.UserData = rmfield(popupH.UserData,'HighLight');
         end
          
         Old_Y;
         pH.YData=NewY' ;
         fH.UserData.ColRow=NewColRawByBase;
         fH.UserData.CylinderOrder=NewYOrderList;
         
         staplePlotAx2;
         for stpi=1: length(staplePlotAx2)
            ScafXY= [ pH.XData' , Old_Y];
            ThisStapBeforeAlign= [ staplePlotAx2{stpi}.XData ; staplePlotAx2{stpi}.YData-0.2 ]' ;
             [~, IndInScaf] = ismember(round(ThisStapBeforeAlign), ScafXY,'rows') ;
%              IndInScaf=setdiff(IndInScaf,0)
%             if  sum(IndInScaf==0)~=0
%                 sdfsfgdgd=4
%             end
            IndInScaf_Tep=IndInScaf  ;  
%             IndInScaf_Tep(IndInScaf_Tep==0)=1 ;
% stpi
% if stpi==22;
%     sdg=3
% end
            IndInScaf_Tep(IndInScaf_Tep==0 )=[];    %fixed overhang cases, Oct 2 2018
            Movement=  NewY(IndInScaf_Tep)' - Old_Y(IndInScaf_Tep)' ;
            Movement(IndInScaf'==0)=0;
             
             staplePlotAx2{stpi}.YData= staplePlotAx2{stpi}.YData +  Movement ;
         end
         
         

         mmX=floor(min(pH.XData/32))*32-32 ;
         pH.XData=pH.XData-mmX ;
         
          for stpi=1: length(staplePlotAx2)
         staplePlotAx2{stpi}.XData=  staplePlotAx2{stpi}.XData-mmX ;
         end
%          XYAll=[pH.XData ;pH.YData]' ;
%          Ymovement=zeros(size(BundlesMove,1),1) ; sumY=Ymovement;
%          XYAll_temp=XYAll ;
%          nloop22=0;
%          while 1
%              %--------check diretionality is same
%              Effect=1;
%              Cyls=unique(XYAll_temp(:,2));
%              for cyli=1:length(Cyls)
%                  Bases= XYAll_temp( XYAll_temp(:,2)==Cyls(cyli) ,1) ;
%                  Neibors= find(abs(Bases(1:end-1)-Bases(2:end))==1) ;
%                  QQ=Bases(Neibors+1)- Bases( Neibors) ;
%                  if sum(QQ==-1)~= length(QQ) && sum(QQ==1)~= length(QQ)
%                  Effect=0;ss
%                   break;
%                  end
%              end
%              %---------
%              if Effect==1;    break;   end
%              R_Move=randi([0,2],size(Ymovement,1),1) ;
%              R_Move(R_Move==2)=-1;
%              
%              Ymovement= Ymovement +R_Move ;
%              Ymovement= mod(Ymovement,2) ;
%              for kq=1:length(Ymovement)
%                XYAll_temp(BelongTransM==kq,2)=  XYAll_temp(BelongTransM==kq,2) + Ymovement(kq);
%              end
%              sumY=Ymovement+sumY;
%              nloop22=nloop22+1;
%          end
%          nloop22;
%          mmX=min(XYAll_temp(:,1)) ; mmY=min(XYAll_temp(:,2))-2 ;
%          Xmove=  floor(mmX/32)*32 -32 ;
%          
%          
         
%          
%          pH.XData=XYAll_temp(:,1)'-Xmove ; 
%          pH.YData=XYAll_temp(:,2)' -mmY ;        
%          popupH.UserData.HighLight.YData=  pH.YData(popupH.UserData.HighLight.UserData)';
%          popupH.UserData.HighLight.XData=  pH.XData(popupH.UserData.HighLight.UserData)';
%         %-----------------------
%          CylsNew=unique(XYAll_temp(:,2));
%         
%         New_NumList=zeros(size(CylsNew,1),4) ;
%         New_NumList(:,1)=CylsNew;

         ssfsdf=124;
         
    end


    function lineselect(src,evn,popupH,pH , staplePlotAx2, fH ,BelongTransM,ax2)
  
        xy=evn.IntersectionPoint(1:2);
        
        XYAll=[pH.XData ;pH.YData]' ;
        
        dXY= XYAll-  ones(size(XYAll,1),1)*xy ;
        d=dXY(:,1).^2 +dXY(:,2).^2 ;
        Ind= find(d==min(d));
        BelongTransM(Ind)
        popupH.Value= BelongTransM(Ind) ;
        popupFcn(popupH,[],pH, staplePlotAx2, fH ,BelongTransM,ax2  );
    end
  
      function keyMove(src,evn,fH,pH , staplePlotAx2,BelongTransM,popupH  ,editH, ax2) 
            GindexNeedTomove= find(BelongTransM ==   popupH.Value);
                    TransVec=[];
                     LinearIncrement= round(str2double(editH.String));    
                     if isnan(LinearIncrement)
                         return
                     end
                     
                        switch evn.Character
                          case 'a'
                          TransVec= -[32*LinearIncrement,0,] ;
                          case 'd'
                          TransVec=[32*LinearIncrement,0,]   ;    
                          case 'w'
                          TransVec=[0,-LinearIncrement]      ;         
                          case 's'
                           TransVec=[0,LinearIncrement]    ;      
                          otherwise
                                return
                        end
                   if ~isempty(TransVec)
                %----scaf
%                     Inscaf=intersect(pH.UserData,GindexNeedTomove) ;
                     pH.XData(GindexNeedTomove)= pH.XData(GindexNeedTomove)+ TransVec(1);
                     pH.YData(GindexNeedTomove)= pH.YData(GindexNeedTomove)+ TransVec(2);
                     if isfield(popupH.UserData,'HighLight' )
                     popupH.UserData.HighLight.XData=popupH.UserData.HighLight.XData+ TransVec(1);
                     popupH.UserData.HighLight.YData=popupH.UserData.HighLight.YData+ TransVec(2);
                     end
                   end

                %-----
                   for stpi=1:length(staplePlotAx2)
                       Indstp= popupH.Value ==   staplePlotAx2{stpi}.UserData.BelongM ;
                       
                     staplePlotAx2{stpi}.XData(Indstp)= staplePlotAx2{stpi}.XData(Indstp)+ TransVec(1);
                     staplePlotAx2{stpi}.YData(Indstp)= staplePlotAx2{stpi}.YData(Indstp)+ TransVec(2);                       
                   end
                
                
                   
%          axis equal; 
         xlim auto; ylim auto;
              switch evn.Character
                  case 'r'
                  ax2.XLimMode='manual';
                   XXX= ax2.XLim +2;    xlim(XXX);
                  case 'f'
                   ax2.XLimMode='manual';
                  XXX= ax2.XLim -2;     xlim(XXX)
                  case 't'
                  ax2.YLimMode='manual';
                   YYY= ax2.YLim +2;    ylim(YYY);
                  case 'g'
                  ax2.YLimMode='manual';
                  YYY= ax2.YLim -2;    ylim(YYY);

              end     
      end
  
    function popupFcn(src,evn,pH , staplePlotAx2, fH ,BelongTransM ,ax2 )
         GindexNeedTomove= find(BelongTransM ==   src.Value);
%         axes=ax2
         if ~isfield(src.UserData,'HighLight' )
         XX=pH.XData(GindexNeedTomove);
         YY=pH.YData(GindexNeedTomove);
             src.UserData.HighLight=scatter(ax2,XX,YY,'.r');
             src.UserData.HighLight.UserData=GindexNeedTomove;
         else
         XX=pH.XData(GindexNeedTomove);
         YY=pH.YData(GindexNeedTomove);
         src.UserData.HighLight.XData=XX;
         src.UserData.HighLight.YData=YY;
         
         src.UserData.HighLight.UserData=GindexNeedTomove;
         
         end
         
%                     if ~isempty(intersect( pHss{kstp}.UserData,GindexNeedTomove))
%                         pHss{1}.MarkerSize=6;
%                         pHss{1}.Color=[ 0,0.5,0.6] ;
%                     else
%                         pHss{1}.MarkerSize=10;
%                         pHss{1}.Color=[ 1,0,0] ;
%                     end
%                    
    end
  
          function checkFcn2(src,~,aH)     
         % used to change axis equal/normal with chechbox
            if src.Value==1
                axis(aH,'normal')
            else
                 axis(aH,'equal')
            end
        end
  
end

