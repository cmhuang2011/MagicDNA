
f41H=figure(41);clf ;
% subplot(1,2,1);
hold on; axis  equal;

AllGXYZ=cell(1,length(obj.containBundle));
for k=1:length(obj.containBundle)
AllGXYZ{k}=obj.containBundle{k}.CylinderXYZGlobal ;
C1Center=mean([AllGXYZ{k}(1,1:3);AllGXYZ{k}(1,1:3);AllGXYZ{k}(1,4:6)]);
f41H.UserData.textH{k}=text( C1Center(1),C1Center(2),C1Center(3)+5, num2str(k),'FontSize',36);
end 

 SaveGHelix=cell(1,length(obj.containBundle));
for Bundlei=1:length(obj.containBundle)
 QQWithPM10Bases =obj.containBundle{Bundlei}.HelixXYZG;     
SaveGHelix{Bundlei}= QQWithPM10Bases;
end

 SaveGHelixBVec=cell(1,length(obj.containBundle)); %in old, before-simul coordinate
for Bundlei=1:length(obj.containBundle)
 QQWithPM10Bases =obj.containBundle{Bundlei}.HelixXYZGBVec;     
SaveGHelixBVec{Bundlei}= QQWithPM10Bases;
 end
 

SacfR=obj.ScafRouting ;
SacfR;
plotXYZ=zeros(size(SacfR));
for k=1:size(SacfR,1)
   bundle=SacfR(k,1);  Cyl=SacfR(k,2);
   alpha=SacfR(k,3)- obj.containBundle{bundle}.Zbase1(Cyl);
   beta=obj.containBundle{bundle}.Zbase2(Cyl)-SacfR(k,3);
   P= AllGXYZ{bundle}(Cyl,1:3);
   Q=AllGXYZ{bundle}(Cyl,4:6);
   XYZ=(beta*P + alpha*Q )/(alpha+beta);
   plotXYZ(k,:)=XYZ;   
end

SSplotXYZ=size(plotXYZ) ;
% plot3(plotXYZ(:,1), plotXYZ(:,2), plotXYZ(:,3) )
x = plotXYZ(:,1)';
y = plotXYZ(:,2)';
z = plotXYZ(:,3)';
col = (1:size(plotXYZ,1))*1000;  % This is the color
surface([x;x],[y;y],[z;z],[col;col], 'facecol','no', 'edgecol','interp', 'linew',2);
%     DSh=plot3(XYZdata(:,1),-XYZdata(:,2),XYZdata(:,3),'Linewidth',2);   %draw route
scatter3(plotXYZ(1,1),plotXYZ(1,2),plotXYZ(1,3),'s','filled','SizeData',100);   %mark head
scatter3(plotXYZ(end,1),plotXYZ(end,2),plotXYZ(end,3),'d','filled','SizeData',100);  %mark tail
%  xlim auto ; ylim auto; zlim auto;
axis equal;grid on;  xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;


  %------------------------------------------------------------------  
 f3H=figure(43);clf ;
 hold on; axis  equal;
 SacfR;
MaxBase=max(SacfR(:,3));
skipPattern1=9:60:MaxBase;   % skip mod2 =0   %test 4/20
skipPattern2=39:60:MaxBase;
%             skipPattern1=[];   % skip mod2 =0
%             skipPattern2=[];
 %-----------Jan 9,   plotting simulation helix result
 AllHelixRouting=zeros(10000,3); kc=1;
 AllHelixRouting2=zeros(10000,3);
 AppedixAllHelR=zeros(10000,3);  % [Bundle , Cyl, Base];
 AppedixBVec=zeros(10000,3);
 AppedixNVec=zeros(10000,3);
 skipinScaf=0;
for edgeinSCR=1:2:size(SacfR,1)-1
%     edgeinSCR
     bundle=SacfR(edgeinSCR,1);  Cyl=SacfR(edgeinSCR+1,2); Bundle=obj.containBundle{bundle};
     BaseStart=SacfR(edgeinSCR,3);   BaseEnd=SacfR(edgeinSCR+1,3);
     RelativeBS=BaseStart-obj.containBundle{bundle}.Zbase1(Cyl)+11; 
     RelativeBE=BaseEnd-obj.containBundle{bundle}.Zbase1(Cyl)+11; 
     QQ=linspace(RelativeBS,RelativeBE,abs(RelativeBE-RelativeBS)+1);
     
     if strcmp(Bundle.Lattice, 'Square') 
         if BaseStart>BaseEnd  %go to left------------------
            skipinScaf=skipinScaf+length(intersect(QQ, skipPattern2-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable'));  
            QQ=setdiff(QQ, skipPattern2-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable');
            skipP=skipPattern2;      
         else
           skipinScaf=skipinScaf+length(intersect(QQ, skipPattern1-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable'));   
           QQ=setdiff(QQ, skipPattern1-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable');  
           skipP=skipPattern1;
         end
     else
         skipP=[];        %only square lattice needs skip
     end     
%      
     PartHelix=SaveGHelix{bundle}{Cyl}(QQ,:);
     PartBvec= SaveGHelixBVec{bundle}{Cyl}(QQ,:);
     if ~isempty(obj.containBundle{bundle}.SimulateTransMFromTM2)  %if had done simulation -----------------------
         TAll=obj.containBundle{bundle}.SimulateTransMFromTM2;
%          QQQ=TAll(1:3,1:3)'*PartHelix'-TAll(1:3,4)*ones(1, size(PartHelix,1)  );  %-old
         QQQ=TAll(1:3,1:3)*PartHelix'+TAll(1:3,4)*ones(1, size(PartHelix,1)  );
         AllHelixRouting2(kc:kc+size(PartHelix,1)-1,: )=QQQ' ;
%          AppedixAllHelR(kc:kc+size(PartHelix,1)-1,: )= [ones(size(PartHelix,1),1)*bundle, ones(size(PartHelix,1),1)*Cyl , linspace(BaseStart,BaseEnd,1+abs(BaseStart- BaseEnd))']   ;
         AppedixAllHelR(kc:kc+size(PartHelix,1)-1,: )= [ones(size(PartHelix,1),1)*bundle, ones(size(PartHelix,1),1)*Cyl ,setdiff( linspace(BaseStart,BaseEnd,1+abs(BaseStart- BaseEnd))', skipP ,'stable' ) ]   ;
         AppedixBVec(kc:kc+size(PartBvec,1)-1,: )=(TAll(1:3,1:3)'*PartBvec')';   %rotate only, relative Vec
         if BaseStart>BaseEnd
         NVecinleft=AllGXYZ{bundle}(Cyl,4:6)-AllGXYZ{bundle}(Cyl,1:3);NVecinleft=NVecinleft/norm(NVecinleft);
         else
         NVecinleft=AllGXYZ{bundle}(Cyl,1:3)-AllGXYZ{bundle}(Cyl,4:6) ; NVecinleft=NVecinleft/norm(NVecinleft);
         end
         NVecinleftArray=ones(size(PartBvec,1),1)*NVecinleft;
         AppedixNVec(kc:kc+size(PartBvec,1)-1,: )=(TAll(1:3,1:3)'*NVecinleftArray')';   %rotate only, relative Vec
         
     end
     PartHelixV2=PartHelix;
     PartHelixV2(:,1:2)=1.3*PartHelix(:,1:2);%---------debug
     AllHelixRouting(kc:kc+size(PartHelix,1)-1,: )=PartHelixV2 ;
     kc=kc+size(PartHelix,1);
%      SaveGHelix{bundle}        
     PreVbumdle=bundle ;
end
        AllHelixRouting( sum(AllHelixRouting,2)==0,:)=[];
        RemovedInd=sum(AllHelixRouting2,2)==0;
        AllHelixRouting2( RemovedInd,:)=[];
        AppedixAllHelR(RemovedInd,:)=[];
        AppedixBVec(RemovedInd,:)=[];
        AppedixNVec(RemovedInd,:)=[];
        
        if size(AppedixAllHelR,1)<8064
        p8064=p8064Seq   ;
        ScafSeqATCG=p8064(1:size(AppedixAllHelR,1));
        else
        ScafSeqATCG=randseq(size(AppedixAllHelR,1));
        end
subplot(1,2,1) ;hold on; axis  equal;

plot3(  AllHelixRouting(:,1),AllHelixRouting(:,2),AllHelixRouting(:,3) ,'b' );
    x = AllHelixRouting(:,1)'; y = AllHelixRouting(:,2)'; z = AllHelixRouting(:,3)';
%     col = (1:size(AllHelixRouting,1))*1000;  % This is the color
%     surface([x;x],[y;y],[z;z],[col;col],'facecol','no', 'edgecol','interp', 'linew',2);      
     subplot(1,2,2);hold on; axis  equal;
%       scatter3(  AllHelixRouting2(:,1),AllHelixRouting2(:,2),AllHelixRouting2(:,3),'.');
     plot3(  AllHelixRouting2(:,1),AllHelixRouting2(:,2),AllHelixRouting2(:,3) ,'.-b' );  
    x = AllHelixRouting2(:,1)'; y = AllHelixRouting2(:,2)'; z = AllHelixRouting2(:,3)';
    col = (1:size(AllHelixRouting2,1))*1000;  % This is the color
%     surface([x;x],[y;y],[z;z],[col;col],'facecol','no', 'edgecol','interp','linew',2);
       
    %-----------------------------staple part
    
    
 SaveGHelixStap=cell(1,length(obj.containBundle));
for Bundlei=1:length(obj.containBundle)
    QQWithPM10Bases =obj.containBundle{Bundlei}.HelixXYZGStap;     
    SaveGHelixStap{Bundlei}= QQWithPM10Bases;
end
  
       
SaveGHelix;
  
StapHelixCell=cell(1,length(obj.StapList3) );
StapHelix2Cell=cell(1,length(obj.StapList3) );
AppedixStapHelRCell=cell(1,length(obj.StapList3) );
StapSeqATCGCell=cell(1,length(obj.StapList3) );
AppedixAllHelR;
AppedixBVecCell=cell(1,length(obj.StapList3) );
AppedixNVecCell=cell(1,length(obj.StapList3) );
 skipinStap=0;

  for stai=1:length(obj.StapList3)
 StapAll=obj.StapList3{stai};
 StapHelix=zeros(10000,3); kc=1;   % pre-allocate
 StapHelix2=zeros(10000,3);
 AppedixStapHelR=zeros(10000,3);  % [Bundle , Cyl, Base];
    for edgeinSCR=1:2:size(StapAll,1)
         C5Cyl=StapAll(edgeinSCR,1);
         bundle=obj.RelateTable(obj.RelateTable(:,5)==C5Cyl,1);   %culti-section needs be cautious
         Bundle=obj.containBundle{bundle};
         Cyl=obj.RelateTable(obj.RelateTable(:,5)==C5Cyl,2);

    %      bundle=StapAll(edgeinSCR,1);  Cyl=StapAll(edgeinSCR+1,2);
         BaseStart=StapAll(edgeinSCR,2);   BaseEnd=StapAll(edgeinSCR+1,2);
         if BaseStart==BaseEnd
             sdfsf=23434;
         end
            if   length(Cyl)~=1
                bundle=unique(bundle);
                for whicCyl=1:length(Cyl)
                lB1=obj.containBundle{bundle}.Zbase1(Cyl(whicCyl));
                lB2=obj.containBundle{bundle}.Zbase2(Cyl(whicCyl)) ;   
                [lB1,lB2];
                if BaseStart>=lB1 && BaseStart<=lB2
                    Cyl=Cyl(whicCyl);
                    break;
                end
                end
            end
         RelativeBS=BaseStart-obj.containBundle{bundle}.Zbase1(Cyl)+11 ;
         RelativeBE=BaseEnd-obj.containBundle{bundle}.Zbase1(Cyl)+11;
         QQ=linspace(RelativeBS,RelativeBE,abs(RelativeBE-RelativeBS)+1);
         
         if  strcmp(Bundle.Lattice, 'Square') 
             if BaseStart<BaseEnd  %go to left------------------stap
                 skipinStap=skipinStap+length(intersect(QQ, skipPattern2-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable'));
                QQ=setdiff(QQ, skipPattern2-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable');
                skipP=skipPattern2;
             else
                  skipinStap=skipinStap+length(intersect(QQ, skipPattern1-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable')); 
                 QQ=setdiff(QQ, skipPattern1-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable');  
                   skipP=skipPattern1;
             end
         else
             skipP=[];        %only square lattice needs skip
         end     
        
         
         PartHelix2=SaveGHelixStap{bundle}{Cyl}(QQ,:);
%         skipP=[];
         if ~isempty(obj.containBundle{bundle}.SimulateTransMFromTM2)  %if had done simulation ----------------
             TAll=obj.containBundle{bundle}.SimulateTransMFromTM2;
    %          QQQ=TAll(1:3,1:3)'*PartHelix2'-TAll(1:3,4)*ones(1, size(PartHelix2,1)  );%---old
             QQQ=TAll(1:3,1:3)*PartHelix2'+TAll(1:3,4)*ones(1, size(PartHelix2,1)  );
             StapHelix2(kc:kc+size(PartHelix2,1)-1,: )=QQQ' ;
             AppedixStapHelR(kc:kc+size(PartHelix2,1)-1,: )= [ones(size(PartHelix2,1),1)*bundle, ones(size(PartHelix2,1),1)*Cyl , setdiff(linspace(BaseStart,BaseEnd,1+abs(BaseStart- BaseEnd))', skipP ,'stable')  ]   ;
         end
         
%          setdiff( linspace(BaseStart,BaseEnd,1+abs(BaseStart- BaseEnd))', skipP ,'stable' )
         
         PartHelix2V2=PartHelix2;
         PartHelix2V2(:,1:2)=1.3*PartHelix2V2(:,1:2);       %---------debug
         StapHelix(kc:kc+size(PartHelix2,1)-1,: )=PartHelix2V2 ;     
         kc=kc+size(PartHelix2,1);
    %      SaveGHelix{bundle}        
    end
        StapHelix( sum(StapHelix,2)==0,:)=[];  %delete extra
        RemovedInd=sum(StapHelix2,2)==0;
        StapHelix2( RemovedInd,:)=[];
        AppedixStapHelR(RemovedInd,:)=[];  
        
        StapHelixCell{stai}=StapHelix;
        StapHelix2Cell{stai}=StapHelix2;
        AppedixStapHelRCell{stai}=AppedixStapHelR;
        
        [Incase,SeqInd]=ismember(AppedixStapHelR,AppedixAllHelR,'rows');
        if nnz(Incase)~=length(Incase)   %has single-stranded staple
            sdfsf=234;
        end
        
%         QWE=ScafSeqATCG(SeqInd)
        StapSeqATCGCell{stai}=  seqcomplement(ScafSeqATCG(SeqInd));   % complementary, A<->T, C<->G.------------
        AppedixBVecCell{stai}=-AppedixBVec(SeqInd,:);   %find scaffold BVec, add negative
        AppedixNVecCell{stai}=-AppedixNVec(SeqInd,:);   %find scaffold NVec, add negative
       [stai, length(StapSeqATCGCell{stai})];
         subplot(1,2,1);
  plot3(  StapHelix(:,1),StapHelix(:,2),StapHelix(:,3));
     subplot(1,2,2);
    plot3(  StapHelix2(:,1),StapHelix2(:,2),StapHelix2(:,3));
     
  end
%   toc
AllHelixRouting;
AllHelixRouting2;
AppedixAllHelR;  % [Bundle , Cyl, Base];
AppedixBVec;
AppedixNVec;
 
%  AllHelixRouting2=AllHelixRouting
 
 StapHelix2Cell;AppedixBVecCell;AppedixNVecCell;
 
 ScafSeqATCG;
 

skipinStap ;
skipinScaf ;
 
%------------- Add extra bp on  bridging points -Mar 8
%  ForCompensate1=AllHelixRouting2(2:end,:)-AllHelixRouting2(1:end-1,:)  ;
%  ddFComp= sqrt(ForCompensate1(:,1).^2+ForCompensate1(:,2).^2+ForCompensate1(:,3).^2) ;
%  thresholdToAddExtraBase=0.7  ;
%  IncremToAdd=1;
%  IndTooFar=ddFComp>thresholdToAddExtraBase ;
%  TooFar=find(IndTooFar);
%  qq=0 ;
%  UnitedGiven2bps=2;
%  FilterTooFar=zeros(size(TooFar));
%  for cc=1:sum(IndTooFar)
%      if AppedixAllHelR(TooFar(cc),1)~= AppedixAllHelR(1+TooFar(cc),1)  % cross bundle
%     WW= AllHelixRouting2(TooFar(cc),:)- AllHelixRouting2(1+TooFar(cc),:) ;
%     norm(WW);
%     qq=qq+1 ;
%        FilterTooFar(cc)=  TooFar(cc);
%      end
%  end
%  
%  
%  FilterTooFar(FilterTooFar==0)=[];   %Cross-Bundle ForcedConnection
%  for k=1:length(FilterTooFar)   %either points is end bridge point
%      UP3= AppedixAllHelR( FilterTooFar(k),:);
%     ALLup3Base=  AppedixAllHelR(and(AppedixAllHelR(:,1)==UP3(1),AppedixAllHelR(:,2)==UP3(2)),3);
%     if  UP3(3)==min(ALLup3Base) || UP3(3)==max(ALLup3Base)
%         FilterTooFar(k)=0;
%     end
%      Next3= AppedixAllHelR( FilterTooFar(k)+1,:);
%     ALLNext3Base=  AppedixAllHelR(and(AppedixAllHelR(:,1)==Next3(1),AppedixAllHelR(:,2)==Next3(2)),3);
%      if  Next3(3)==min(ALLNext3Base) || Next3(3)==max(ALLNext3Base)
%         FilterTooFar(k)=0;   
%     end
%  end
%  FilterTooFar(FilterTooFar==0)=[];
%  
%  OriInterval=union(FilterTooFar, FilterTooFar+ones(size(FilterTooFar))) ;
%  EE=[1;size(AllHelixRouting2,1)];
%   OriInterval=union(OriInterval,EE) ;  %used to reassemble data,  since need to insert multiple sections
% 
%   NewAppedixPVec=zeros( size(AllHelixRouting2,1) +  UnitedGiven2bps*length(FilterTooFar) ,3) ;
%   NewAppedixBVec=NewAppedixPVec ;
%   NewAppedixNVec=NewAppedixPVec ;
% 
%   newSeq(1, length(ScafSeqATCG)+ UnitedGiven2bps*length(FilterTooFar)) = char(0);
%   for plugin=1:2: length(OriInterval)
%     NewAppedixPVec( plugin-1+OriInterval(plugin):plugin-1+OriInterval(plugin+1),:)= AllHelixRouting2(  OriInterval(plugin):OriInterval(plugin+1),:) ;
%     NewAppedixBVec(  plugin-1+OriInterval(plugin):plugin-1+OriInterval(plugin+1),:)= AppedixBVec(  OriInterval(plugin):OriInterval(plugin+1),:) ;
%     NewAppedixNVec(  plugin-1+OriInterval(plugin):plugin-1+OriInterval(plugin+1),:)= AppedixNVec(  OriInterval(plugin):OriInterval(plugin+1),:) ;
%     newSeq( plugin-1+OriInterval(plugin):plugin-1+OriInterval(plugin+1))=ScafSeqATCG( OriInterval(plugin):OriInterval(plugin+1));
%   end
% ININ=find(NewAppedixPVec(:,1)==0);   %  IndNeedInsertNew
% for k=1:2:length(ININ)
%     
%    OVecAll=[NewAppedixPVec(ININ(k)-1,:) , NewAppedixBVec(ININ(k)-1,:) , NewAppedixNVec(ININ(k)-1,:) ;...
%             NewAppedixPVec(ININ(k)+2,:) , NewAppedixBVec(ININ(k)+2,:) , NewAppedixNVec(ININ(k)+2,:) ] ;
%    Insert1=2/3*OVecAll(1,:) + 1/3*OVecAll(2,:)  ;
%    Insert2=1/3*OVecAll(1,:) + 2/3*OVecAll(2,:) ;
%    Insert=[Insert1;Insert2] ;
%    Insert(1,4:6)=  Insert(1,4:6)/norm(Insert(1,4:6)) ;  % unit vec
%    Insert(1,7:9)=  Insert(1,7:9)/norm(Insert(1,7:9)) ;
%    Insert(2,4:6)=  Insert(2,4:6)/norm(Insert(2,4:6)) ;
%    Insert(2,7:9)=  Insert(1,7:9)/norm(Insert(1,7:9)) ;
%     
%    NewAppedixPVec(ININ(k):ININ(k)+1,1:3)= Insert(:,1:3);
%    NewAppedixBVec(ININ(k):ININ(k)+1,1:3)= Insert(:,4:6);
%    NewAppedixNVec(ININ(k):ININ(k)+1,1:3)= Insert(:,7:9);
%    newSeq(ININ(k))='T' ;
%    newSeq(ININ(k)+1)='T' ;
% end
  %++++++++++end of  extra bp on  bridging points -Mar 8
%   AllHelixRouting2  ;   AppedixBVec  ; AppedixNVec; ScafSeqATCG;
%   AllHelixRouting2=NewAppedixPVec ;
%   AppedixBVec= NewAppedixBVec;
%   AppedixNVec = NewAppedixNVec ;
%   ScafSeqATCG = newSeq ;
 %-------------------end of March 8 adding
 
 
 return
 
 mmAll=min(AllHelixRouting)-[10 ,10,10] ;
 MMAll=max(AllHelixRouting);
 
 %----------------
 UnitCoeff=0.85; CoB=0.33;
 file_name=strcat('T',  num2str(randi([10,99])));
%  file2_name=strcat(file_name,'.conf');
 file2_name='prova.conf';
fileID = fopen(file2_name,'w');
fprintf(fileID,'t = 0\n');
boxsize=2*abs(mmAll-MMAll);
E0=0;
fprintf(fileID,'b = %9.6f %9.6f %9.6f\n',max(boxsize),max(boxsize),max(boxsize));
fprintf(fileID,'E = %8.6f %8.6f %8.6f\n',E0,E0,E0);
%------------% scaffold

% AllHelixRouting2=AllHelixRouting;%----------export left plot
% StapHelix2Cell=StapHelixCell;

QQ=ones(size(AllHelixRouting2,1),1); count=size(AllHelixRouting2,1);
fprintf(fileID,'%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f\n',[ ((AllHelixRouting2-QQ*mmAll )'-CoB*AppedixBVec')/UnitCoeff ;-AppedixBVec';-AppedixNVec' ;zeros(size(AllHelixRouting2))';zeros(size(AllHelixRouting2))'  ] );
 %---------% staple
for stpi=1: length(StapHelix2Cell)
      BVecM= AppedixBVecCell{stpi};

  PosiM= ( StapHelix2Cell{stpi}-ones(size(StapHelix2Cell{stpi},1),1)*mmAll -CoB*BVecM )/UnitCoeff   ;
  NVecM= AppedixNVecCell{stpi}; 
fprintf(fileID,'%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f\n',[PosiM' ;-BVecM';-NVecM' ;zeros(size(PosiM))';zeros(size(PosiM))'  ] );
   count=count+size(PosiM,1);
end
%--------
fclose(fileID);

%------------------topology file
%  stpi=0  ;  %------
%  ScafSeqATCG
 
count;
StapSeqATCGCell;

%  file3_name=strcat('Topo',file_name,'.dat');
 file3_name='prova.top' ;
fileID2 = fopen(file3_name,'w');
fprintf(fileID2,'%u %u\n',count,stpi+1 );
%-----initial
%----scaf
fprintf(fileID2,'%u %c %d %u\n', 1, ScafSeqATCG(1),length(ScafSeqATCG)-1   , 1)  ;
for scafi=2:length(ScafSeqATCG)-1   
fprintf(fileID2,'%u %c %u %u\n', 1, ScafSeqATCG(scafi),scafi-2,scafi )  ;
end
fprintf(fileID2,'%u %c %u %d\n', 1, ScafSeqATCG(1),scafi-1, 0) ; 
%------staple

ss=scafi-1;
for stpi2=1: length(StapHelix2Cell)
StapiSeq=StapSeqATCGCell{stpi2};
    
    fprintf(fileID2,'%u %c %d %u\n', stpi2+1, StapiSeq(1),-1, ss+3)  ;
    for stapiline=2:length(StapiSeq)-1   
    fprintf(fileID2,'%u %c %u %u\n', stpi2+1, StapiSeq(stapiline),ss+stapiline ,ss+2+stapiline )  ;
    end
    fprintf(fileID2,'%u %c %u %d\n', stpi2+1, StapiSeq(end),ss+stapiline+1, -1) ; 
ss=ss+stapiline+1;
end
fclose(fileID2);

  dfgdg=23425;
        
        