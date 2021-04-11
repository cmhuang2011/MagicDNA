function  oxDNAPatInitial(src,evn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rotate3d off ;
ss_Assembly= findobj(0,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;

ax= findobj(0,'Tag','oxDNA') ; 
axPat= findobj(0,'Tag','oxDNA2Pat') ; 
popOxdnaPattern= findobj(0,'Tag','popOxdnaPattern') ; 
btn_oxDNAPat2 =  findobj(0,'Tag','btn_oxDNAPat2') ; 
btn_oxDNAPat3 =  findobj(0,'Tag','btn_oxDNAPat3') ; 

t_json=findobj(0,'Tag','t_json') ;
tData = t_json.Data ;



axes(axPat);cltab; hold on ;axis equal;

popOxdnaPattern.Value=1;

% Pattern_choice = choosePattern ;
% Pattern_choice.mode='Chain';
% Circular
Pattern_choice.mode='Circular';

Pattern_choice.n=4 ;

str2 = cellfun(@num2str,num2cell(1:Pattern_choice.n),'uniformoutput',0);
popOxdnaPattern.String=str2;

popOxdnaPattern.Max=Pattern_choice.n ;


PrecScaf = ax.UserData.scaf ;
PrecStap = ax.UserData.stap ;
PrecCS = ax.UserData.CS ;

% new_handle = copyobj(h,p)

NewScaf =cell(Pattern_choice.n ,3) ;  % [backboneHelx , BasesCenter , NVec ] ; 
NewStap =cell(Pattern_choice.n ,3) ;
NewCS= cell(Pattern_choice.n ,3) ;
%--------
mergeStapleH = mergeGO( PrecStap.pStap , axPat ) ; %mergeStapleH.Color=[1,0,0];
mergeStapleHCenter = mergeGOscatter(  PrecStap.pStap_center , axPat ) ; %mergeStapleH.Color=[1,0,0];

mergeScafH = mergeGO( PrecScaf.pScaf2 , axPat ) ; %mergeStapleH.Color=[1,0,0];
mergeScafHCenter = mergeGOscatter(  PrecScaf.pScaf_center , axPat ) ; %mergeStapleH.Color=[1,0,0];



mergeCSH =  mergeGO( PrecCS.pCS , axPat ) ; %mergeStapleH.Color=[1,0,0];

for k=1:      length(PrecCS.pCS_center)         
PrecCS.pCS_center{k}.UserData.Seq=  tData{length(GetHyperB.StapList3)+1+k ,3} ;
end


mergeCSHCenter = mergeGOscatter(   PrecCS.pCS_center , axPat ) ; %mergeStapleH.Color=[1,0,0];

NewScaf{1,1}= mergeScafH; 
NewScaf{1,2}= mergeScafHCenter; 


NewStap{1,1}= mergeStapleH; 
NewStap{1,2}= mergeStapleHCenter; 

NewCS{1,1} = mergeCSH;
NewCS{1,2} = mergeCSHCenter;

Ccolor = get(0, 'DefaultAxesColorOrder') ;
 
for i_go=1: size(NewScaf,1)
%     ScafFields= fieldnames(PrecScaf)    ; % retrieve field 1 ,3 ,5  
%     NewScaf{i_go,1} = copyobj( PrecScaf.(ScafFields{1}){1},axPat) ;
%     NewScaf{i_go,2} = copyobj( PrecScaf.(ScafFields{3}){1},axPat) ;
    
    
    if i_go~=1     
        NewScaf{i_go,1} = copyobj( mergeScafH,axPat) ;    
        NewScaf{i_go,2} = copyobj( mergeScafHCenter,axPat) ;    
        
        NewStap{i_go,1} = copyobj( mergeStapleH,axPat) ;    
        NewStap{i_go,2} = copyobj( mergeStapleHCenter,axPat) ;    

        
        NewCS{i_go,1} = copyobj( mergeCSH,axPat) ;    
        NewCS{i_go,2} = copyobj( mergeCSHCenter,axPat) ;    

    end
        NewScaf{i_go,2}.Visible = 'off'; 
        NewStap{i_go,2}.Visible = 'off'; 
        
      NewScaf{i_go,1}.Color =  Ccolor(mod(i_go, size(Ccolor,1))+1 ,:) ;
      NewStap{i_go,1}.Color =  Ccolor(mod(i_go, size(Ccolor,1))+1 ,:) ;
end

for k =1: size(NewScaf,1)
     NewScaf{k,1}.ButtonDownFcn =@(src,evn )SelectScaf(src,evn ,NewScaf ) ;
     NewStap{k,1}.ButtonDownFcn =@(src,evn )SelectScaf(src,evn ,NewStap ) ;
     
end


mergeStpNVec = zeros(sum(cellfun('length',PrecStap.NVecstap))+length(PrecStap.NVecstap)-1 ,3) ;
nStpBase = 1 ;
for k=1:length(PrecStap.NVecstap)
    Inds=nStpBase:nStpBase+length(PrecStap.NVecstap{k})-1 ;
    mergeStpNVec(Inds,: ) =    PrecStap.NVecstap{k} ;
    nStpBase=nStpBase+length(PrecStap.NVecstap{k})+1 ;
end

mergeCSNVec = zeros(sum(cellfun('length',PrecCS.NVecCS))+length(PrecCS.NVecCS)-1 ,3) ;
nStpBase = 1 ;
for k=1:length(PrecCS.NVecCS)
    Inds=nStpBase:nStpBase+length(PrecCS.NVecCS{k})-1 ;
    mergeCSNVec(Inds,: ) =    PrecCS.NVecCS{k} ;
    nStpBase=nStpBase+length(PrecCS.NVecCS{k})+1 ;
end

Offset = [0,0 ,0] ;

for movei=1:size(NewScaf,1)
    NewScaf{movei,1}.XData =  NewScaf{movei,1}.XData + movei*Offset(1) ;
    NewScaf{movei,2}.XData =  NewScaf{movei,2}.XData + movei*Offset(1) ;
    NewScaf{movei,2}.UserData.NVec =PrecScaf.NVecscaf{1} ;
%     NewScaf{movei,3} =PrecScaf.NVecscaf ;

    
    NewStap{movei,1}.XData =  NewStap{movei,1}.XData + movei*Offset(1) ;
    NewStap{movei,2}.XData =  NewStap{movei,2}.XData + movei*Offset(1) ;
    NewStap{movei,2}.UserData.NVec =mergeStpNVec ;
%     NewStap{movei,3} =mergeStpNVec ;


    NewCS{movei,1}.XData =  NewCS{movei,1}.XData + movei*Offset(1) ;
    NewCS{movei,2}.XData =  NewCS{movei,2}.XData + movei*Offset(1) ;
%     NewCS{movei,3}=PrecCS.NVecCS;
    
    
end


mergeCSHAll =  mergeGO( NewCS(:,1), axPat ) ; mergeCSHAll.Color=[0,1,0];
mergeCSHAll_Center =  mergeGOscatter( NewCS(:,2), axPat ) ; mergeCSHAll_Center.CData=[0,1,0];

mergeCSNVecAll=repmat([mergeCSNVec; [0,0,0]] ,size(NewScaf,1),1) ;
mergeCSNVecAll(end-1,:)=[];
for k=1:size(NewCS,1)
delete(NewCS{k,1}) ; delete(NewCS{k,2}) ;
end
% All closing strands across instances are merged into one GO(graphic object). 

sdfsf=3;
   
   
% s=1 ;PrecCS.CSLengths;
% q=1;
% nTotalLin = size(PrecCS.CSLengths,1) * Pattern_choice.n  ;
Interval = sum(sum(PrecCS.CSLengths))+1 +size(PrecCS.CSLengths,1)-1 ;

CSBaseBM = -1*ones(size(mergeCSHAll.XData)) ;
for k=1:Pattern_choice.n   
    ins=(1:Interval)+Interval*(k-1);
    CSBaseBM(ins)=k ;
end
w=1;Savearr= cell(size(PrecCS.CSLengths,1)*size(PrecCS.CSLengths,2),1 ) ;m=1;
for jRow=1: size(PrecCS.CSLengths,1) 
     for jCol=1: size(PrecCS.CSLengths,2) 
          a= PrecCS.CSLengths(jRow,jCol) ;
          arrs = w:a+w-1 ;
          Savearr{m}=arrs ;m=m+1;
          w=w+a;
     end
     w=w+1;
end

% Savearr([3,4])=Savearr([4,3]) ; % hard

CSBaseAll = -1*ones(1,length(mergeCSHAll.XData)) ;
kline=1 ;
for i=1: Pattern_choice.n-1
    for j=1:size(PrecCS.CSLengths,1)
        L_inds=Savearr{2*j-1} +(i-1)*Interval ;
        Cellind=mod(2*j,2*size(PrecCS.CSLengths,1));
        if Cellind==0;Cellind=2*size(PrecCS.CSLengths,1);end
        R_inds=Savearr{Cellind  }; R_inds=R_inds+1*Interval +(i-1)*Interval ;
        
%         if i ==1
        CSBaseAll(L_inds)=kline;
        CSBaseAll(R_inds)=kline;
%         else
            
%         end
        
        kline=kline+1 ;
    end
end
% IsRing=1 ;

if  strcmp( Pattern_choice.mode,'Circular')==1
   for jring=1:size(PrecCS.CSLengths,1)
      L_inds=Savearr{2*jring}  ;
      Cellind=mod(2*jring-1,2*size(PrecCS.CSLengths,1));
       if Cellind==0;Cellind=2*size(PrecCS.CSLengths,1);end 
        R_inds=Savearr{Cellind  }; R_inds=R_inds+(i)*Interval  ;
        CSBaseAll(L_inds)=kline;
        CSBaseAll(R_inds)=kline;
        kline=kline+1 ;
   end    
end

% uv = unique(CSBaseAll);
% nCount  = histc(CSBaseAll,uv) ;
 
CollectXYZ = [mergeCSHAll.XData' , mergeCSHAll.YData' ,mergeCSHAll.ZData' ];
CollectXYZ_Center = [mergeCSHAll_Center.XData' , mergeCSHAll_Center.YData' ,mergeCSHAll_Center.ZData' ];
CollectSeq =mergeCSHAll_Center.UserData.Seq ;
% nP1 = find(
% mergeCSHAll.XData
% NewCollectXYZ=zeros(size(CollectXYZ)) ;
% sn=1;
% for k=1:max(CSBaseBM)
%    Ind =  CSBaseBM==k ;
%     
%    NewCollectXYZ(sn:sn+sum(Ind)-1 , :) = CollectXYZ(Ind,:) ;
%    sn=sn+sum(Ind)+1 ;
% end
% 
% NewCollectXYZ(NewCollectXYZ==0) =NaN ;
% 
% mergeCSHAll.XData=NewCollectXYZ(:,1);
% mergeCSHAll.YData=NewCollectXYZ(:,2) ;
% mergeCSHAll.ZData=NewCollectXYZ(:,3) ;


NewCollectXYZ=zeros(size(CollectXYZ)) ; 
NewCollectXYZ_Center=zeros(size(CollectXYZ_Center)) ;
NewCSSeq = char(size(CollectSeq)) ;

ShuffleBM=  -1*ones(size(CSBaseBM));
% CSBaseBM(CSBaseBM==Pattern_choice.n+1) =1;
sn=1;
for k=1:max(CSBaseAll)
   Ind =  CSBaseAll==k ;    
   NewCollectXYZ(sn:sn+sum(Ind)-1 , :) = CollectXYZ(Ind,:) ;
   NewCollectXYZ_Center(sn:sn+sum(Ind)-1 , :) = CollectXYZ_Center(Ind,:) ;
   mergeCSNVecAll(sn:sn+sum(Ind)-1 , :) = mergeCSNVecAll(Ind,:) ;
   NewCSSeq(sn:sn+sum(Ind)-1 ) = CollectSeq(Ind) ;
   
   
   
   ShuffleBM(sn:sn+sum(Ind)-1) = CSBaseBM(Ind) ;
   sn=sn+sum(Ind)+1 ;
   
end

NewCollectXYZ(NewCollectXYZ==0) =NaN ;
NewCollectXYZ_Center(NewCollectXYZ_Center==0) =NaN ;

% 
% uv = unique(ShuffleBM)
% nCount  = histc(ShuffleBM,uv) 



mergeCSHAll.XData=NewCollectXYZ(:,1);
mergeCSHAll.YData=NewCollectXYZ(:,2) ;
mergeCSHAll.ZData=NewCollectXYZ(:,3) ;
mergeCSHAll.UserData.BM =ShuffleBM ;

mergeCSHAll_Center.XData=NewCollectXYZ_Center(:,1);
mergeCSHAll_Center.YData=NewCollectXYZ_Center(:,2) ;
mergeCSHAll_Center.ZData=NewCollectXYZ_Center(:,3) ;
mergeCSHAll_Center.UserData.BM =ShuffleBM ;
mergeCSHAll_Center.UserData.NVec =mergeCSNVecAll ;

mergeCSHAll_Center.UserData.Seq =NewCSSeq;


sdfsf=3 ;
% uv = unique(ShuffleBM);
% nCount  = histc(ShuffleBM,uv) ;

% figure;plot(ShuffleBM,mergeCSHAll.XData);  %hold on;yyaxis right; plot(mergeCSHAll.XData)


fH=gcf ;

fH.KeyPressFcn=@(src,evn) keypPressForPattern(src,evn,NewScaf,NewStap,mergeCSHAll,mergeCSHAll_Center)   ; % temporary change to other keypress fcn
fH.UserData.saveKeyMove2 = @(src,evn) keypPressForPattern(src,evn,NewScaf,NewStap,mergeCSHAll,mergeCSHAll_Center ) ;


% figure(66);clf;yyaxis left; plot(CSBaseBM);  hold on;yyaxis right; plot(mergeCSHAll.XData)

sss_pattern= findobj(0,'Tag','sss_pattern') ;
sss_pattern.ButtonDownFcn=@(src,evn)RecoverKeyMove(src,evn,2) ; 
xlabel('X'); ylabel('Y') ;zlabel('Z') ;

btn_oxDNAPat2.Callback= @(src,evn)ExportPatOxDNA(src,evn,NewScaf ,NewStap ,mergeCSHAll,mergeCSHAll_Center ,Pattern_choice) ;


end


function ExportPatOxDNA(src,evn,NewScaf ,NewStap ,mergeCSHAll,mergeCSHAll_Center ,Pattern_choice) 
fprintf('start printing pattern oxDNA \n ')
tic
ss_Assembly= findobj(0,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;

t_json=findobj(0,'Tag','t_json') ;
tData = t_json.Data ;
n_pattern =  size(NewScaf,1 ) ;
n_closingstrand =  size(tData,1)-1 -length(GetHyperB.StapList3) ;

scafL = repmat(length( NewScaf{1,1}.XData ) ,1,n_pattern) ;
StpL = repmat(cellfun('length',t_json.Data(1:length(GetHyperB.StapList3),3))' ,1,n_pattern) ;

CSall = bwlabel(~isnan( mergeCSHAll.XData) );

scafSeq =GetHyperB.pSeq{1}(1:scafL(1)) ;
    saveallseq = char(ones(1,sum(scafL)+sum(StpL)+sum(CSall~=0) )) ;

%     fprintf('%u %u\n',sum(scafL)+sum(StpL)+sum(CSall~=0),length(scafL)+length(StpL)+ max(CSall) );  % header
% 
%     fprintf('%u %u %u %u %u %u\n',sum(scafL),sum(StpL),sum(CSall),length(scafL),length(StpL), max(CSall) );  % header

%---------------------------- topology file
    file_name='provaPat.top' ;
    fileID2 = fopen(file_name,'w');
    fprintf(fileID2,'%u %u\n',sum(scafL)+sum(StpL),length(scafL)+length(StpL));  % header

%     fprintf(fileID2,'%u %u\n',sum(scafL)+sum(StpL)+sum(CSall~=0),length(scafL)+length(StpL)+ max(CSall) );  % header
    %----scaf---------
    gb =1 ; gs =1 ;
    for rp=1: n_pattern
        fprintf(fileID2,'%u %c %d %u\n', gs, scafSeq(1), -1   , gb)  ;  gb=gb+1 ; 
        for scafi=2:length(scafSeq)-1   
            fprintf(fileID2,'%u %c %u %u\n', gs, scafSeq(scafi),gb-2,gb )  ; gb=gb+1 ;
        end
        fprintf(fileID2,'%u %c %u %d\n', gs, scafSeq(scafi+1),gb-2, -1) ;  gb=gb+1 ; gs=gs+1;
        saveallseq(gb-length(scafSeq):gb-1) =scafSeq ;
% sdfsf=3
%     end    
    %------staple-----------------
%     for rp=1: n_pattern
        for stpi2=1:length(GetHyperB.StapList3)
            StapiSeq=tData{stpi2,3};

                fprintf(fileID2,'%u %c %d %u\n', gs, StapiSeq(1),-1, gb)  ; gb=gb+1 ;
                for stapiline=2:length(StapiSeq)-1   
                fprintf(fileID2,'%u %c %u %u\n', gs, StapiSeq(stapiline),gb-2 ,gb)  ; gb=gb+1 ;
                end
                fprintf(fileID2,'%u %c %u %d\n', gs, StapiSeq(end),gb-2, -1) ;  gb=gb+1 ;  gs=gs+1;
                
            saveallseq(gb-length(StapiSeq):gb-1) =StapiSeq ;
            if ~isempty(findstr(StapiSeq,'?'))
            showerror=1 
            end   
        end   
    end
     %------closing strand-----------
% 
%      for csi = 1 :max(CSall)
%      Inds =CSall ==csi ;
%      SeqsCS=mergeCSHAll_Center.UserData.Seq(Inds) ;
% %      Csmod = mod(csi,n_closingstrand) ; Csmod(Csmod==0)=n_closingstrand ;
% %      SeqsCS = tData{length(GetHyperB.StapList3)+1+Csmod , 3} ;
%      
%             fprintf(fileID2,'%u %c %d %u\n', gs, SeqsCS(1),-1, gb)  ; gb=gb+1 ;
%             for stapiline=2:length(SeqsCS)-1   
%             fprintf(fileID2,'%u %c %u %u\n', gs, SeqsCS(stapiline),gb-2 ,gb)  ; gb=gb+1 ;
%             end
%             fprintf(fileID2,'%u %c %u %d\n', gs, SeqsCS(end),gb-2, -1) ;  gb=gb+1 ;  gs=gs+1;     
%      saveallseq(gb-length(SeqsCS):gb-1) =SeqsCS ;           
%      end
fclose(fileID2);

%------------------------------- configuration 
% NewScaf ,NewStap ,mergeCSHAll,mergeCSHAll_Center 
sdfsf=3  ;

    CentersVec = zeros(sum(scafL)+sum(StpL) , 3); 
    BVec = zeros(sum(scafL)+sum(StpL) , 3); 
    NVec = zeros(sum(scafL)+sum(StpL), 3);     


%     CentersVec = zeros(sum(scafL)+sum(StpL)+sum(CSall~=0) , 3); 
%     BVec = zeros(sum(scafL)+sum(StpL)+sum(CSall~=0) , 3); 
%     NVec = zeros(sum(scafL)+sum(StpL)+sum(CSall~=0), 3);     
    %-------scaffold 
     gb2 =1 ;
   
        for rp=1: n_pattern
            CenterXYZ =[NewScaf{rp,2}.XData' , NewScaf{rp,2}.YData' , NewScaf{rp,2}.ZData'] ;
            CurrentBaseInd = gb2:gb2-1+size(CenterXYZ,1) ;
            CentersVec( CurrentBaseInd,: ) =CenterXYZ ;
            BackBone = [NewScaf{rp,1}.XData' , NewScaf{rp,1}.YData' , NewScaf{rp,1}.ZData'] ;

            QQ= CenterXYZ- BackBone ;
            d=  sqrt(QQ(:,1).^2 +  QQ(:,2).^2 +  QQ(:,3).^2 ) ;
            BVec(CurrentBaseInd ,:) = QQ/mean(d)   ; 

            NVec(CurrentBaseInd ,:) =NewScaf{rp,2}.UserData.NVec ;
            gb2=gb2+size(CenterXYZ,1) ;
%         end
    %------------- staple
%          for rp=1: n_pattern
            CenterXYZ =[NewStap{rp,2}.XData' , NewStap{rp,2}.YData' , NewStap{rp,2}.ZData'] ; 
            CenterXYZ(isnan(CenterXYZ(:,1)),: )=[]   ;% remove NaN ;
            CurrentBaseInd = gb2:gb2-1+size(CenterXYZ,1) ;
            CentersVec( CurrentBaseInd,: ) =CenterXYZ ;
            
            BackBone = [NewStap{rp,1}.XData' , NewStap{rp,1}.YData' , NewStap{rp,1}.ZData'] ;
            IndKeep =~isnan(BackBone(:,1)) ;
             BackBone(isnan(BackBone(:,1)),: )=[]   ;% remove NaN ;
            QQ= CenterXYZ- BackBone ;
            d=  sqrt(QQ(:,1).^2 +  QQ(:,2).^2 +  QQ(:,3).^2 ) ;
            BVec(CurrentBaseInd ,:) = QQ/mean(d)   ; 
            
            Ncc=NewStap{rp,2}.UserData.NVec ;
            NVec(CurrentBaseInd ,:) = Ncc(IndKeep ,:) ;   %remove [0,0,0];
       
            gb2=gb2+size(CenterXYZ,1) ;
         end
     %----------------closing strand   
%      mergeCSHAll_Center; mergeCSHAll;
% count =0;
%         for csi = 1 :max(CSall)
%             Ind = CSall==csi ;
%             CenterXYZ =[mergeCSHAll_Center.XData(Ind)' , mergeCSHAll_Center.YData(Ind)' , mergeCSHAll_Center.ZData(Ind)'] ; 
%             CurrentBaseInd = gb2:gb2-1+size(CenterXYZ,1) ;
%             CentersVec( CurrentBaseInd,: ) =CenterXYZ ;
%             
%             BackBone = [mergeCSHAll.XData(Ind)' , mergeCSHAll.YData(Ind)' , mergeCSHAll.ZData(Ind)'] ; 
%             QQ= CenterXYZ- BackBone ;
%             d=  sqrt(QQ(:,1).^2 +  QQ(:,2).^2 +  QQ(:,3).^2 ) ;
%             BVec(CurrentBaseInd ,:) = QQ/mean(d)   ;
%             
%             Ncc=mergeCSHAll_Center.UserData.NVec( Ind,:) ;
%             if sum(sum(Ncc,2)==0)>1
%                 sfsff=33 ;
%             end
%             
%             NVec(CurrentBaseInd ,:) = Ncc ;
%             
%             gb2=gb2+size(CenterXYZ,1) ;
%             count=count+size(CenterXYZ,1) ;
%         end
        
    bounds = [ min(CentersVec) ;  max(CentersVec)]  ;
    boxCenter =0.5*(bounds(1,:) +bounds(2,:) )  ;
    boxCenter = mean(CentersVec ) ;
    
    Dbox = bounds(2,:) -  bounds(1,:) ;
    KK=1;
    boxsize =max( KK*50*ceil(Dbox/50) ) ; 
    cc= 0.5*[boxsize,boxsize,boxsize] ;
%     CentersVec = CentersVec  -boxCenter  ;
    CentersVec = CentersVec  -boxCenter  ;
%     CentersVec = CentersVec  -boxCenter +cc ;

    file2_name='provaPat33.conf';
    fileID = fopen(file2_name,'w');
    fprintf(fileID,'t = 0\n'); E0=0;
    fprintf(fileID,'b = %9.6f %9.6f %9.6f\n',1.8*max(boxsize),1.8*max(boxsize),1.8*max(boxsize));
    fprintf(fileID,'E = %8.6f %8.6f %8.6f\n',E0,E0,E0);

    for k=1:  size(CentersVec,1)
    fprintf(fileID,'%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f\n',CentersVec(k,1:3)/0.85,BVec(k,1:3),NVec(k,1:3),zeros(1,6) );
    end
    fclose(fileID);
    fprintf(' finished printing conf and topology in pattern \n '  )       

    fprintf(' Make sure monomer had been exported for dsRemain update in class ...\n')
    GetHyperB.ExportExtradsRemain(n_pattern ,Pattern_choice.mode) ;
    
    %------------------dsRemain
    %   Coeff=0.8 ;delta=0.01;
%   trials = 0.5:0.001:0.7;
%   results=zeros(length(trials),2) ;
%  for k=1:length(trials)
%     Coeff= trials(k);
%     CentersVec2=CentersVec+Coeff*BVec ;
%     C = uniquetol(CentersVec2,0.00005,'ByRows',true); 
% %  
% %    CentersVecRR=CentersVec+(Coeff+delta)*BVec ;
% %    CRR = uniquetol(CentersVec2,0.0001,'ByRows',true); size(CRR)
% %    
% %     CentersVecLL=CentersVec+(Coeff-delta)*BVec ;
% %     CLL = uniquetol(CentersVec2,0.0001,'ByRows',true); size(CLL)
%    results(k,:) = [trials(k), size(C,1)] ;
% 
%  end
%  figure(77882);clf;plot( results(:,1) ,results(:,2));
%  sdfsf=3
%  Result show :[0.589,tol=0.00005]
    
    
    
% sdfs=3 


%     %----------
%     Coeff=0.589 ;
%     CentersVec2=CentersVec+Coeff*BVec ;
%     [C,IA,IC] = uniquetol(CentersVec2,0.00005,'ByRows',true,'OutputAllIndices',true ); 
%     LL_IA= cellfun('length',IA) ;
%     ScafStapCorrelate2 =zeros(sum(LL_IA ==2) ,2 ); cj=1;
%     Iamwroong=0;
%      for j=1:length(LL_IA)
%          if LL_IA(j)==2
%              ScafStapCorrelate2(cj,:) =(IA{j})' ; cj=cj+1;
%     %          checkSeq =(IA{j})'
%     %          [
%              SeqPair=saveallseq((IA{j})');
%              if strcmp(SeqPair,'AT')||strcmp(SeqPair,'TA')||strcmp(SeqPair,'CG')||strcmp(SeqPair,'GC')
%                  sdf=3;
%              else
%                 fprintf('%i %i %s %i \n',IA{j} ,SeqPair,Iamwroong)
%                  Iamwroong=Iamwroong+1;
%              end
% 
%          end    
%      end
%      ScafStapCorrelate2=ScafStapCorrelate2-1 ; % index difference between Matlab and python
%      fprintf('printing deRemain \n')
%      file3_name='dSRemainPat.conf'   ;
%      fileID2 = fopen(file3_name,'w');
%      
%      for iF=1:size(ScafStapCorrelate2,1)
%             fprintf(fileID2,'{\n' );
%             fprintf(fileID2,'type = mutual_trap\n' );
%             fprintf(fileID2,'particle = %u\n' ,ScafStapCorrelate2(iF,1));
%             fprintf(fileID2,'ref_particle  = %u\n' ,ScafStapCorrelate2(iF,2));
%             fprintf(fileID2,'stiff = %u \n' ,100 );
%             fprintf(fileID2,'r0 = 1.2 \n'  );
%             fprintf(fileID2,'}\n' );  
%             fprintf(fileID2,'{\n' );
%             fprintf(fileID2,'type = mutual_trap\n' );
%             fprintf(fileID2,'particle = %u\n' ,ScafStapCorrelate2(iF,2));
%             fprintf(fileID2,'ref_particle  = %u\n' ,ScafStapCorrelate2(iF,1));
%             fprintf(fileID2,'stiff = %u \n' ,100 );
%             fprintf(fileID2,'r0 = 1.2 \n' );    
%             fprintf(fileID2,'}\n' );  
%      end
%      fclose(fileID2);

toc
    fprintf(' finished printing pattern oxdna formats \n '  )       

end


function SelectScaf(src,evn ,NewScaforStap )

for k =1: size(NewScaforStap,1)
    NewScaforStap{k,1}.LineWidth = 0.5 ; 
    if isequal(  NewScaforStap{k,1} ,src)
        SelectOne = k ;
    end
end
%  src.LineWidth = 2 ; 

%  sdsdf=3
popOxdnaPattern= findobj(0,'Tag','popOxdnaPattern') ; 

if evn.Button==1
    popOxdnaPattern.Value = SelectOne ;
elseif evn.Button==3
    popOxdnaPattern.Value = union( popOxdnaPattern.Value,  SelectOne) ;
end

for k =1: size(NewScaforStap,1)
   if ismember(k,popOxdnaPattern.Value )
        NewScaforStap{k,1}.LineWidth = 2; 
   end    
end

end



% 
function  keypPressForPattern(src,evn,NewScaf,NewStap,mergeCSHAll, mergeCSHAll_Center )
popOxdnaPattern= findobj(0,'Tag','popOxdnaPattern') ; 
checkH_oxDNAPat= findobj(0,'Tag','checkH_oxDNAPat') ; 
editH_oxDNAPat= findobj(0,'Tag','editH_oxDNAPat') ; 
stepsize=str2num(editH_oxDNAPat.String) ;
SelectGO =popOxdnaPattern.Value ;

if checkH_oxDNAPat.Value==1   % translation

    switch evn.Key
        case 'q'  % +X
            field='XData' ; InCrement= stepsize ;
        case 'w' % +Y
            field='YData' ; InCrement=stepsize ;
        case 'e' % +Z
            field='ZData' ; InCrement=stepsize ;
        case 'a' % -X
            field='XData' ; InCrement=-stepsize ;
        case 's' % -Y
            field='YData' ; InCrement=-stepsize ;
        case 'd'  % -Z      
            field='ZData' ; InCrement=-stepsize ;
        otherwise
            return
    end
            
            Inds= ismember(mergeCSHAll.UserData.BM , SelectGO ) ;% Inds(end-1)=[];
%             find(Inds);
            mergeCSHAll.(field)(Inds)=  mergeCSHAll.(field)(Inds)+InCrement ;
            mergeCSHAll_Center.(field)(Inds)=  mergeCSHAll_Center.(field)(Inds)+InCrement ;
            
            for k = 1: length(SelectGO)
                SGO = SelectGO(k) ;
                NewScaf{SGO,1}.(field) =   NewScaf{SGO,1}.(field)+InCrement ;
                NewScaf{SGO,2}.(field) =   NewScaf{SGO,2}.(field)+InCrement ;
                
                NewStap{SGO,1}.(field) =   NewStap{SGO,1}.(field)+InCrement ;
                NewStap{SGO,2}.(field) =   NewStap{SGO,2}.(field)+InCrement ;
            end
else
    theta= stepsize*pi/180 ;
     switch evn.Key
        case 'q'  % +X
                       RMat=[1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
        case 'w' % +Y
                        RMat=[cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
        case 'e' % +Z
                       RMat=[cos(theta) sin(theta) 0; -sin(theta)  cos(theta)  0; 0 0 1];
        case 'a' % -X
                       RMat=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
        case 's' % -Y
                       RMat=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
        case 'd'  % -Z      
                       RMat=[cos(theta) -sin(theta) 0; sin(theta)  cos(theta)  0; 0 0 1];
        otherwise
            return               
     end
    
    
    
     GCindv = zeros(length(SelectGO) ,3 ) ;
     for jj = 1: length(SelectGO)
     GCindv(jj,:) =  mean([NewScaf{SelectGO(jj),1}.XData' , NewScaf{SelectGO(jj),1}.YData', NewScaf{SelectGO(jj),1}.ZData']);
     end
     RotationalCenter = mean(GCindv ,1) ;
    
%     RotationalCenter = mean([NewScaf{SelectGO,1}.XData' , NewScaf{SelectGO,1}.YData', NewScaf{SelectGO,1}.ZData']);
    for jj2 = 1: length(SelectGO)
        SGO = SelectGO(jj2) ;
        for k=1:2
            XYZ1 =[NewScaf{SGO,k}.XData' , NewScaf{SGO,k}.YData', NewScaf{SGO,k}.ZData'] ;
            XYZ1 = XYZ1 - ones(size(XYZ1,1),1)*RotationalCenter;  XYZ1=XYZ1* RMat ;
            XYZ1=XYZ1 + ones(size(XYZ1,1),1)*RotationalCenter ;
            NewScaf{SGO,k}.XData= XYZ1(:,1)' ;     NewScaf{SGO,k}.YData= XYZ1(:,2)'  ;   NewScaf{SGO,k}.ZData= XYZ1(:,3)';
        end
        NewScaf{SGO,2}.UserData.NVec=NewScaf{SGO,2}.UserData.NVec*RMat ;
    end
    
    
    
    for jj2 = 1: length(SelectGO)
        SGO = SelectGO(jj2) ;
        for k=1:2
            XYZ1 =[NewStap{SGO,k}.XData' , NewStap{SGO,k}.YData', NewStap{SGO,k}.ZData'] ;
            XYZ1 = XYZ1 - ones(size(XYZ1,1),1)*RotationalCenter;  XYZ1=XYZ1* RMat ;
            XYZ1=XYZ1 + ones(size(XYZ1,1),1)*RotationalCenter ;
            NewStap{SGO,k}.XData= XYZ1(:,1)' ;     NewStap{SGO,k}.YData= XYZ1(:,2)'  ;   NewStap{SGO,k}.ZData= XYZ1(:,3)';
        end
        NewStap{SGO,2}.UserData.NVec=NewStap{SGO,2}.UserData.NVec*RMat ;
        
    end
%   

    
    Inds= ismember(mergeCSHAll.UserData.BM , SelectGO ) ;% Inds(end-1)=[];
    
%     Inds= mergeCSHAll.UserData.BM==SelectGO ;% Inds(end-1)=[];
    XYZCS= [ mergeCSHAll.XData(Inds)' ,mergeCSHAll.YData(Inds)' ,mergeCSHAll.ZData(Inds)' ,];
    XYZCS = XYZCS - ones(size(XYZCS,1),1)*RotationalCenter;  XYZCS=XYZCS* RMat ;  
    XYZCS=XYZCS + ones(size(XYZCS,1),1)*RotationalCenter ;
    mergeCSHAll.XData(Inds)= XYZCS(:,1)' ;       mergeCSHAll.YData(Inds)= XYZCS(:,2)' ;     mergeCSHAll.ZData(Inds)= XYZCS(:,3)' ;   
   
%     Inds= mergeCSHAll_Center.UserData.BM==SelectGO ;% Inds(end-1)=[];
    XYZCS= [ mergeCSHAll_Center.XData(Inds)' ,mergeCSHAll_Center.YData(Inds)' ,mergeCSHAll_Center.ZData(Inds)' ,];
    XYZCS = XYZCS - ones(size(XYZCS,1),1)*RotationalCenter;  XYZCS=XYZCS* RMat ;  
    XYZCS=XYZCS + ones(size(XYZCS,1),1)*RotationalCenter ;
    mergeCSHAll_Center.XData(Inds)= XYZCS(:,1)' ;       mergeCSHAll_Center.YData(Inds)= XYZCS(:,2)' ;     mergeCSHAll_Center.ZData(Inds)= XYZCS(:,3)' ;   
    mergeCSHAll_Center.UserData.NVec(Inds,:)=mergeCSHAll_Center.UserData.NVec(Inds,:)*RMat ;
    
end
%           NewScaf{SelectGO,2}.UserData.NVec
%          NewScaf{SelectGO,3}{1}
        


end




