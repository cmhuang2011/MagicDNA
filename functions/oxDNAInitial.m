function oxDNAInitial(src,evn)
% src-> btn of initial oxDNA, main-> Mech-> oxDNA -> Conf and force
%
%   Detailed explanation goes here

tic
fH=gcf;
ss_Assembly= findobj(fH,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;
% skipBase= GetHyperB.skipBase ;
ax= findobj(fH,'Tag','oxDNA') ;
btn2= findobj(fH,'Tag','btn_oxDNA2') ;   %export configuration and topology
btn3= findobj(fH,'Tag','btn_oxDNA3') ;  %export build
oxDNASlider =  findobj(fH,'Tag','slider_oxdna') ;
oxDNACheck =  findobj(fH,'Tag','check_oxdna') ;

SimulateAssembly(GetHyperB,[],[])   ;
axes(ax);cltab; hold on ;axis equal;


%%

% pStapleH=plotStaple_Helix( GetHyperB  ) ;   % plot staple strand
% pScaf=plotScaf2_Helix( GetHyperB,{GetHyperB.scafC5}  ) ;     % plot scaf strand

% for k=1:length(pStapleH);  pStapleH{k}.Color=[1,0,0]; end
% pScaf{1}.Color =[0,0,1] ;



% pScaf2=plotScaf2_Helix_V2( GetHyperB,{GetHyperB.scafC5}  ) ;     % plot scaf strand
Isstap = 0; TM=2 ;
[pScaf2,ScafHelix,pScaf_center,ScafBaseCenterHelix ,NVecscaf,scafBundleRout ]=plotScaf2_Helix_V2( GetHyperB,GetHyperB.scafC5,Isstap ,[0,0,1] ,TM ) ;     % plot scaf strands
for k=1:length(pScaf2) % MagicDNA2 , insert adjust
 AdjustInsert( pScaf2{k} , pScaf_center{k} ) ;
end


% return
[pStap,StapHelix,pStap_center,StapBaseCenterHelix, NVecstap,stapBundleRout ]=plotScaf2_Helix_V2( GetHyperB,GetHyperB.StapList3,1 ,[1,0,0],TM ) ;     % plot staple strand
for k=1:length(pStap)  % MagicDNA2 , insert adjust
 AdjustInsert( pStap{k} , pStap_center{k} ) ;
end


% sdfs=3

scaf.pScaf2=pScaf2 ; scaf.ScafHelix=ScafHelix; scaf.pScaf_center=pScaf_center ; scaf.ScafBaseCenterHelix=ScafBaseCenterHelix; scaf.NVecscaf=NVecscaf ;
stap.pStap=pStap ; stap.StapHelix=StapHelix; stap.pStap_center=pStap_center ; stap.StapBaseCenterHelix=StapBaseCenterHelix; stap.NVecstap=NVecstap ;


% [pScaf2,ScafHelix,pScaf_center,ScafBaseCenterHelix]=plotScaf2_Helix_V2( GetHyperB,GetHyperB.ClosingCornerNotation,Isstap ,[0.2,0.2,0.2]  ) ;     % plot scaf strands

[pCS,CShelix,pCS_center,CSBaseCenterHelix ,NVecCS,CSBundleRout ]=plotScaf2_Helix_V2( GetHyperB,GetHyperB.ClosingCornerNotation,Isstap ,[0.2,0.2,0.2],TM  ) ;     % plot scaf strands

CS.pCS=pCS ; CS.CShelix=CShelix ; CS.pCS_center=pCS_center ; CS.CSBaseCenterHelix=CSBaseCenterHelix ; CS.NVecCS=NVecCS ;

% GetHyperB.ClosingCornerNotation
CSLengths= zeros(size(GetHyperB.ClosingCornerNotation,1) ,2)  ;
for k=1:size(CSLengths,1)
    if isempty( GetHyperB.ClosingCornerNotation{k})
        continue
    end
    
    CSLengths(k,1) =  size(interpolateBase(  GetHyperB.ClosingCornerNotation{k}(1:2,:) )  ,1 ) ;
    CSLengths(k,2) =  size(interpolateBase(  GetHyperB.ClosingCornerNotation{k}(3:4,:) )  ,1 );
end
CS.CSLengths=CSLengths ;



% sdf=3

oxDNASlider.Callback=@(src,evn)silderexplode(src,evn, GetHyperB ,scaf, stap ,GetHyperB.scafC5,GetHyperB.StapList3,CS,GetHyperB.ClosingCornerNotation ) ;

btn2.Callback= @(src,evn)ExportOxDNA(src,evn,scaf ,stap ,GetHyperB,CS, oxDNACheck,scafBundleRout,stapBundleRout,CSBundleRout) ;
btn3.Callback= @(src,evn)ExportBUILD(src,evn,scaf ,stap ,GetHyperB,CS, oxDNACheck) ;

oxDNACheck.Callback =@(src,evn)IncludeCS(src,evn,CS) ;


ax.UserData.scaf=scaf;
ax.UserData.stap=stap;
ax.UserData.CS=CS;

toc

% scafL = cellfun('length',scafBundleRout)  ;
% scafSeq =cellfun('length',stapBundleRout)  ;
% ClosingStrandL = cellfun('length',CSBundleRout)  ;
% n_BM= sum(scafL)+sum(scafSeq)+sum(ClosingStrandL) ;
% BM =zeros( n_BM,1) ;
% BM(1:sum(scafL)) =scafBundleRout{1} ;
% BM=[ cell2mat(scafBundleRout); cell2mat(stapBundleRout); cell2mat(CSBundleRout)] ;

% sdaf=3


end

function  ExportBUILD(src,evn,scaf ,stap ,GetHyperB,CS, oxDNACheck)

prompt = {'Enter BILD File name:'};    dlg_title = 'Input';  num_lines = 1;
defaultans = {'oxDNA_unrelaxed'};    answer = inputdlg(prompt,dlg_title,num_lines,defaultans) ;

opts.Interpreter = 'tex'; opts.Default = 'coarsed-grained';

if oxDNACheck.Value ==1
    note='Note: with closing strands ';
else
    note='Note: without closing strands ';
end
answer_CG = questdlg({'\fontsize{15} Export CG or ribbong model ?' ,note}  , ...
    'chimera options', 'coarsed-grained','ribbon',opts);

switch answer_CG
    case 'coarsed-grained'
        file2_name=strcat(answer{1},'_CG_', '.bild')     ;
    case 'ribbon'
        file2_name=strcat(answer{1},'_ribbon_' , '.bild')      ;
end


%      c = uisetcolor
tic
fileID = fopen(file2_name,'w');

%      sdfsdf=3
ColorScaf =[0,0,1];
fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorScaf  )    ;
for k=1:length(scaf.pScaf2{1}.XData)-1
    BackBoneA=[scaf.pScaf2{1}.XData(k),scaf.pScaf2{1}.YData(k),scaf.pScaf2{1}.ZData(k)] ;
    BackBoneB=[scaf.pScaf2{1}.XData(k+1),scaf.pScaf2{1}.YData(k+1),scaf.pScaf2{1}.ZData(k+1)] ;
    fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,0.1 )    ;
end

ColorStap =[1,0,0];
fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorStap  )    ;
for k=1: length(stap.pStap)
    for j=1:length(stap.pStap{k}.XData)-1
        BackBoneA=[ stap.pStap{k}.XData(j), stap.pStap{k}.YData(j), stap.pStap{k}.ZData(j)] ;
        BackBoneB=[ stap.pStap{k}.XData(j+1), stap.pStap{k}.YData(j+1), stap.pStap{k}.ZData(j+1)] ;
        fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,0.1 )    ;
    end
end

if oxDNACheck.Value==1
    ColorCS = [0.5,0.5,0.5] ;
    fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorCS  )    ;
    for k=1: length(CS.pCS)
        for j=1:length(CS.pCS{k}.XData)-1
            BackBoneA=[ CS.pCS{k}.XData(j), CS.pCS{k}.YData(j), CS.pCS{k}.ZData(j)] ;
            BackBoneB=[ CS.pCS{k}.XData(j+1), CS.pCS{k}.YData(j+1),CS.pCS{k}.ZData(j+1)] ;
            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,BackBoneB,0.1 )    ;
        end
    end
end

if strcmp(answer_CG, 'coarsed-grained')
    ColorBackBoneToBase =[0.1,0.9,0.9];
    fprintf(fileID , '.color %8.6f %8.6f %8.6f \n',ColorBackBoneToBase  )    ;
    for k=1:length(scaf.pScaf_center{1}.XData)
        Center=[scaf.pScaf_center{1}.XData(k),scaf.pScaf_center{1}.YData(k),scaf.pScaf_center{1}.ZData(k)] ;
        Backbone=[scaf.pScaf2{1}.XData(k),scaf.pScaf2{1}.YData(k),scaf.pScaf2{1}.ZData(k)] ;
        fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',Center,Backbone,0.1 )    ;
    end
    
    for k=1: length(stap.pStap)
        for j=1:length(stap.pStap{k}.XData)
            BackBoneA=[ stap.pStap{k}.XData(j), stap.pStap{k}.YData(j), stap.pStap{k}.ZData(j)] ;
            Center=[  stap.pStap_center{k}.XData(j),  stap.pStap_center{k}.YData(j),  stap.pStap_center{k}.ZData(j)] ;
            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,Center,0.1 )    ;
        end
    end
    
    
    if oxDNACheck.Value ==1
        for k=1: length(CS.pCS)
            for j=1:length(CS.pCS{k}.XData)-1
                BackBoneA=[ CS.pCS{k}.XData(j), CS.pCS{k}.YData(j), CS.pCS{k}.ZData(j)] ;
                Center=[ CS.pCS_center{k}.XData(j), CS.pCS_center{k}.YData(j),CS.pCS_center{k}.ZData(j)] ;
                fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',BackBoneA,Center,0.1 )    ;
            end
        end
        
    end
    
end



fclose(fileID);
fprintf('Exporting BILD takes %4.2f sec \n' ,toc)

end

function IncludeCS(src,evn,CS)
if  src.Value==1
    for k= 1:length(CS.pCS)
        CS.pCS{k}.LineStyle ='-' ;
    end
else
    for k= 1:length(CS.pCS)
        CS.pCS{k}.LineStyle =':' ;
    end
end

end


function ExportOxDNA(src,evn,scaf ,stap ,GetHyperB,CS, oxDNACheck,scafBundleRout,stapBundleRout,CSBundleRout )
%% In fact, this function also export mrDNA file
t_json=findobj(gcf,'Tag','t_json') ;
tData = t_json.Data ;
StpL = cellfun('length',t_json.Data(1:length(GetHyperB.StapList3),3))  ;

% cellfun('length',stap.StapHelix)
% scafL = size(scaf.ScafHelix{1} ,1)  ;
% scafSeq =GetHyperB.pSeq(1:size(scaf.ScafHelix{1} ,1)) ;
% ClosingStrandL = cellfun('length',t_json.Data(length(GetHyperB.StapList3)+2:end ,3))  ;

%     if ~isempty(findstr(scafSeq,'?'))
%     showerroe=1
%     end
%-------------------
if oxDNACheck.Value==1
    BM2=[ cell2mat(scafBundleRout); cell2mat(stapBundleRout); cell2mat(CSBundleRout )] ;
    %      BM2=[ cell2mat(scafBundleRout); cell2mat(stapBundleRout);cell2mat(CSBundleRout(1:12) )] ;
else
    BM2=[ cell2mat(scafBundleRout); cell2mat(stapBundleRout)] ;
end


% %--------------mrDNA_stack
%     mrDNA_stack=[ cell2mat( GetHyperB.ScafAllBase); cell2mat( GetHyperB.StapAllBase)] ;
%     scafsL =size(cell2mat( GetHyperB.ScafAllBase),1 ) ;
%     mrDNA_stack(:,end+1)=-1 ;
%
%     for k=1:size(mrDNA_stack,1)-1
%         if mrDNA_stack(k,1)==mrDNA_stack(k+1,1) && abs(diff(mrDNA_stack(k:k+1,2) ))<=2 % due to skip within one base
%          mrDNA_stack(k,3) =k;
%         end
%
%     end
%     Inds=find(mrDNA_stack(:,3) ==-1) ; % check crossover
%     for k=1:length(Inds)
%     C5_Base_Stack = mrDNA_stack(Inds(k),:)      ;
%         if mod(GetHyperB.RelateVec(C5_Base_Stack(1)),2) ==0  && Inds(k)<=scafsL    % check helical direction: cylinder direction ,  isScaf
%              TargetStack = C5_Base_Stack(1:2) + [0 ,1] ;
%         elseif mod(GetHyperB.RelateVec(C5_Base_Stack(1)),2) ==1  && Inds(k)<=scafsL
%               TargetStack = C5_Base_Stack(1:2) + [0 ,-1] ;
%         elseif mod(GetHyperB.RelateVec(C5_Base_Stack(1)),2) ==0  && Inds(k)>scafsL
%               TargetStack = C5_Base_Stack(1:2) + [0 ,-1] ;
%         elseif mod(GetHyperB.RelateVec(C5_Base_Stack(1)),2) ==1  && Inds(k)>scafsL
%               TargetStack = C5_Base_Stack(1:2) + [0 ,1] ;
%         end
%        [Lia,Locb] = ismember(TargetStack, mrDNA_stack(:,1:2),'rows' ) ;
%        if  Lia==1
%           QQ=find(and(mrDNA_stack(:,1)==TargetStack(1),mrDNA_stack(:,2)==TargetStack(2))) ; % should be two: one is scaffold base, the other is staple
%           for cc=1:length(QQ)
%            if ~xor(Inds(k)>scafsL ,QQ(cc)-1>scafsL) % find same domain one
%                [Inds(k) ,  QQ(cc)];
%                mrDNA_stack(Inds(k) ,3 ) = QQ(cc)-1 ;
%            end
%           end
% %           [Inds(k),  Locb-1]
%        end
%     end
% % %--------------mrDNA_stack



BM=BM2(:,1) ;
fileTransInd_name=strcat( 'BM.mat') ;
save(fileTransInd_name, 'BM') ;

%-------------------for bending structures,  hard code
BM3=BM ;
prompt = {'Enter Curved Bundles:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {''};
answer = inputdlg(prompt,dlgtitle,dims,definput) ;
numArr2=str2num(answer{1}) ; numArr2(numArr2<0)=[]; numArr2(numArr2>max(BM))=[];
% SpiltBunde=[];
SpiltBunde = numArr2;
fprintf('export BM3 which break transMs of bundle [ %s ] into 3 piece \n' ,num2str(SpiltBunde) ) ;
cc=1;
for k=1:length(SpiltBunde)
    AA=BM2(BM2(:,1)==SpiltBunde(k ),2) ;
    %     meanVal = mean(AA) ;
    %     k
    minV=min(AA) ;maxV=max(AA) ;
    ThreePoints= linspace(minV , maxV,4) ;
    
    IndMove1 =  and(and(BM2(:,1)==SpiltBunde(k ) , BM2(:,2)>ThreePoints(2)) ,   BM2(:,2)<=ThreePoints(3) )  ;
    BM3(IndMove1) =max(BM) +cc ; cc=cc+1;
    
    IndMove2 =  and(and(BM2(:,1)==SpiltBunde(k ) , BM2(:,2)>ThreePoints(3)) ,   BM2(:,2)<=ThreePoints(4) )  ;
    BM3(IndMove2) =max(BM) +cc ; cc=cc+1;
    
end


fileTransInd_name=strcat( 'BM3.mat') ;
save(fileTransInd_name, 'BM3') ;
%------------------------


%-------------for grouping bundles
% Grouping={1:3,4:6,7:9,10:12} ;
% HaveAss=[];
% BM3=zeros(size(BM)) ;
% for k=1:length(Grouping)
%  OriInds = Grouping{k} ;
%   IndsInBM= ismember(BM ,OriInds) ;
%    BM3( IndsInBM) =k ;
%    HaveAss=union(OriInds,HaveAss) ;
% end
%
% RestOfBundles = setdiff(1:max(BM), HaveAss) ;
% for cc=1:length(RestOfBundles)
% Inds =ismember(BM ,RestOfBundles(cc)) ;
% BM3( Inds) =k+cc ;
% end
%
% fprintf('grouping bundles  into single pieces \n'  ) ;
% % % -----------------
% SpiltBunde=[];
% BMnew=BM3;
% SpiltBunde = [13 14 15];
% fprintf('export BM3 which break transMs of bundle [ %s ] into 3 piece \n' ,num2str(SpiltBunde) ) ;
% cc=1;
% for k=1:length(SpiltBunde)
%     AA=BM2(BM2(:,1)==SpiltBunde(k ),2) ;
%     meanVal = mean(AA) ;
% %     k
%     minV=min(AA) ;maxV=max(AA) ;
%     ThreePoints= linspace(minV , maxV,4) ;
%
%     IndMove1 =  and(and(BM2(:,1)==SpiltBunde(k ) , BM2(:,2)>ThreePoints(2)) ,   BM2(:,2)<=ThreePoints(3) )  ;
%     BM3(IndMove1) =max(BMnew) +cc ;
% %     max(BM3) +cc
%     cc=cc+1;
%     IndMove2 =  and(and(BM2(:,1)==SpiltBunde(k ) , BM2(:,2)>ThreePoints(3)) ,   BM2(:,2)<=ThreePoints(4) )  ;
%     BM3(IndMove2) = max(BMnew) +cc ; cc=cc+1;
% %     max(BM3) +cc
%
% end
% fileTransInd_name=strcat( 'BM3.mat') ;
% save(fileTransInd_name, 'BM3') ;
%-------------------

% return

%% topology
[scafL,~] = cellfun(@size, scaf.ScafHelix) ;
scafSeq=GetHyperB.pSeq ;
% scafL = size(scaf.ScafHelix{1} ,1)  ;
% scafSeq =GetHyperB.pSeq(1:size(scaf.ScafHelix{1} ,1)) ;
ClosingStrandL = cellfun('length',t_json.Data(length(GetHyperB.StapList3)+2:end ,3))  ;
%      stepAssign=5 ;Ex_DecoratePlatesWithOverhangs ; % selective add closing strands
file3_name='prova.top' ;
fileID2 = fopen(file3_name,'w');
%-----initial---------
if oxDNACheck.Value==0
    fprintf(fileID2,'%u %u\n',sum(scafL)+sum(StpL),length(StpL)+length(scafL) );  % header
    saveallseq = char(ones(1,sum(scafL)+sum(StpL))) ;
    
else
    fprintf(fileID2,'%u %u\n',sum(scafL)+sum(StpL)+sum(ClosingStrandL),length(ClosingStrandL)+length(StpL)+length(scafL) );  % header
    saveallseq = char(ones(1,sum(scafL)+sum(StpL)+sum(ClosingStrandL)  )) ;
end
%         fprintf(fileID2,'%u %u\n',scafL,1 );  % header

%----scaf---------
q=1 ;
for scaf_j = 1 : length(scafSeq)
    saveallseq(q:scafL(scaf_j)-1+q) =scafSeq{scaf_j}  ;
    fprintf(fileID2,'%u %c %d %i\n', scaf_j, scafSeq{scaf_j}(1), q   , -1)  ;  q=q+1;
    %         fprintf(fileID2,'%u %c %d %u\n', 1, scafSeq(1), -1   , 1)  ;
    for scafi=2:length(scafSeq{scaf_j})-1
        fprintf(fileID2,'%u %c %u %u\n', scaf_j, scafSeq{scaf_j}(scafi),q,q-2 )  ; q=q+1;
        %                 fprintf(fileID2,'%u %c %u %u\n', 1, scafSeq(scafi),scafi-2,scafi )  ;
    end
    fprintf(fileID2,'%u %c %i %d\n', scaf_j, scafSeq{scaf_j}(scafi+1),-1, q-2) ;  q=q+1;
    %        fprintf(fileID2,'%u %c %u %d\n', 1, scafSeq(scafi+1),scafi-1, -1) ;    
end
%------staple-----------------
ss=q-1-2;
for stpi2=1:length(GetHyperB.StapList3)
    StapiSeq=tData{stpi2,3};
    %        if strcmp(GetHyperB.ScafOption.prevStack ,'polyT')
    StapiSeq(findstr(StapiSeq,'?') ) ='T' ;   %   poly T design
    %        end    
    fprintf(fileID2,'%u %c %d %i\n', stpi2+length(scafL), StapiSeq(1), ss+3,-1)  ;         %     fprintf(fileID2,'%u %c %d %u\n', stpi2+1, StapiSeq(1),-1, ss+3)  ;
    
    for stapiline=2:length(StapiSeq)-1
        fprintf(fileID2,'%u %c %u %u\n', stpi2+length(scafL), StapiSeq(stapiline),ss+stapiline+2 ,ss+stapiline )  ; %            fprintf(fileID2,'%u %c %u %u\n', stpi2+1, StapiSeq(stapiline),ss+stapiline ,ss+2+stapiline )  ;
    end
    fprintf(fileID2,'%u %c %i %d\n', stpi2+length(scafL), StapiSeq(end),-1,ss+stapiline+1) ;    %    fprintf(fileID2,'%u %c %u %d\n', stpi2+1, StapiSeq(end),ss+stapiline+1, -1) ;
    
    saveallseq(ss+3:ss+3+stapiline) =StapiSeq;
    ss=ss+stapiline+1;
    if ~isempty(findstr(StapiSeq,'?'))
%         showerroe=1
        fprintf('Staple %i has unkwnown sequence which has been substituded as ''T''', stpi2)
    end
end
%------closing strand-----------
if oxDNACheck.Value==1
    for csi2=1:length(ClosingStrandL)
        csiSeq=tData{length(GetHyperB.StapList3)+1+csi2,3};
        
        fprintf(fileID2,'%u %c %d %i\n', csi2+stpi2+scaf_j, csiSeq(1), ss+3,-1)  ;  % fprintf(fileID2,'%u %c %d %u\n', csi2+stpi2+1, csiSeq(1),-1, ss+3)  ;
        for csiline=2:length(csiSeq)-1
            fprintf(fileID2,'%u %c %u %u\n', csi2+stpi2+scaf_j, csiSeq(csiline),ss+csiline+2 ,ss+csiline )  ;     %      fprintf(fileID2,'%u %c %u %u\n', csi2+stpi2+1, csiSeq(csiline),ss+csiline ,ss+2+csiline )  ;
        end
        fprintf(fileID2,'%u %c %i %d\n', csi2+stpi2+scaf_j, csiSeq(end),-1, ss+csiline+1  ) ;          %   fprintf(fileID2,'%u %c %u %d\n', csi2+stpi2+1, csiSeq(end),ss+csiline+1, -1) ;
        
        saveallseq(ss+3:ss+3+csiline) =csiSeq;
        ss=ss+csiline+1;
        
        if ~isempty(findstr(csiSeq,'?'))
%             showerroe=1
            fprintf('Closing strand %i has unkwnown sequence which has been substituded as ''T''', csi2)
        end
    end
end

fclose(fileID2);
%% configuration file
%----------------
% findstr(scafSeq,'?')
if oxDNACheck.Value==0
    CentersVec = zeros(sum(scafL)+sum(StpL) , 3);
    BVec = zeros(sum(scafL)+sum(StpL) , 3);
    NVec = zeros(sum(scafL)+sum(StpL) , 3);
else
    CentersVec = zeros(sum(scafL)+sum(StpL)+sum(ClosingStrandL) , 3);
    BVec = zeros(sum(scafL)+sum(StpL)+sum(ClosingStrandL) , 3);
    NVec = zeros(sum(scafL)+sum(StpL)+sum(ClosingStrandL) , 3);
end

cc=1;
for scaf_j = 1 :length(scafL)
    CentersVec(cc:cc+scafL(scaf_j)-1  ,:) = [scaf.pScaf_center{scaf_j}.XData ;  scaf.pScaf_center{scaf_j}.YData; scaf.pScaf_center{scaf_j}.ZData]'  ;  % slider can change the exported configuration
    BackBone =[scaf.pScaf2{scaf_j}.XData ; scaf.pScaf2{scaf_j}.YData ;scaf.pScaf2{scaf_j}.ZData  ]' ;
    QQ= CentersVec(cc:cc+scafL(scaf_j)-1   ,:) - BackBone ;
    d=  sqrt(QQ(:,1).^2 +  QQ(:,2).^2 +  QQ(:,3).^2 ) ;
    BVec(cc:cc+scafL(scaf_j)-1  ,:) = QQ/mean(d)   ;
    NVec(cc:cc+scafL(scaf_j)-1   ,:) =  scaf.NVecscaf{scaf_j}   ;
    
    cc=cc+scafL(scaf_j) ;
end
%----------
c=0;
for k = 1:  length(stap.pStap)
    nB_stap =StpL(k) ;
    %     if k==18
    %         sdfsdf=3
    %     end
    %     k
    CentersVec(sum(scafL)+c+1:sum(scafL)+c+nB_stap  ,:) =[stap.pStap_center{k}.XData ; stap.pStap_center{k}.YData ; stap.pStap_center{k}.ZData ]' ;
    
    Centeri_stap =    [stap.pStap_center{k}.XData ; stap.pStap_center{k}.YData ; stap.pStap_center{k}.ZData ]' ;
    BackBonei_stap =    [stap.pStap{k}.XData ; stap.pStap{k}.YData ; stap.pStap{k}.ZData ]' ;
    QQ2= Centeri_stap - BackBonei_stap ;
    
    BVec(sum(scafL)+c+1:sum(scafL)+c+nB_stap  ,:) = QQ2/mean(d)    ;
    NVec(sum(scafL)+c+1:sum(scafL)+c+nB_stap  ,:) = stap.NVecstap{k}  ;
    c=c+nB_stap ;
end

if oxDNACheck.Value==1
    
    
    for k = 1:  length(ClosingStrandL)
        nB_cs =ClosingStrandL(k) ;
        CentersVec(sum(scafL)+c+1:sum(scafL)+c+nB_cs  ,:) =[CS.pCS_center{k}.XData ; CS.pCS_center{k}.YData ; CS.pCS_center{k}.ZData ]' ;
        
        Centeri_cs =  [CS.pCS_center{k}.XData ; CS.pCS_center{k}.YData ; CS.pCS_center{k}.ZData ]' ;
        BackBonei_stap =   [CS.pCS{k}.XData ; CS.pCS{k}.YData ; CS.pCS{k}.ZData ]'  ;
        QQ3= Centeri_cs - BackBonei_stap ;
        
        BVec(sum(scafL)+c+1:sum(scafL)+c+nB_cs   ,:) = QQ3/mean(d)    ;
        NVec(sum(scafL)+c+1:sum(scafL)+c+nB_cs   ,:) = CS.NVecCS{k}  ;
        c=c+nB_cs ;
        
    end
end

% CentersVec=CentersVec;
bounds = [ min(CentersVec) ;  max(CentersVec)]  ;
boxCenter =0.5*(bounds(1,:) +bounds(2,:) )  ;
Dbox = bounds(2,:) -  bounds(1,:) ;
KK=1.5;
boxsize =max( KK*50*ceil(Dbox/50) ) ;  cc= 0.5*[boxsize,boxsize,boxsize] ;
CentersVec = CentersVec  -boxCenter +cc ;
% % % %------initial twist configuration 
% ExcInd = [1:981, 7462:7560,  size(CentersVec,1)-1079:size(CentersVec,1)];
% IndAccount=ones(size(CentersVec,1),1)==1 ;
% IndAccount(ExcInd) =false ;
% % 
% CentersVec_ori =CentersVec ;
% ConfToTwist = CentersVec(IndAccount ,:) ;
% BVecToTwist = BVec(IndAccount ,:) ;
% NVecToTwist = NVec(IndAccount ,:) ;
% % 
% N=100 ;
% TotalTwistAngle = -1080  % degree, right-handed as positive
% 
% TotalTwistAngle=-TotalTwistAngle ;
% XSlice=  linspace( min(ConfToTwist(:,1)),max(ConfToTwist(:,1))+0.1, N) ;
% AngArr = linspace(-TotalTwistAngle/2,TotalTwistAngle/2,N-1) ;
% % 
% % [XY , XYcyl]= mappingBend(max(XSlice)-min(XSlice) , 2*(max(ConfToTwist(:,2))-min(ConfToTwist(:,2))) , 2*N-2 , 9 ,120) ;
% % % AngArr=AngArr+min(AngArr);
% % % figure; hold on;
% % for k = 1 : length(XSlice)-1
% % Inds  =   and(ConfToTwist(:,1)>=XSlice(k) ,ConfToTwist(:,1)<=XSlice(k+1)  ) ;
% % RotationAngle = AngArr(k) ;
% % % RotationAngle = 30;
% % % RMat = RotationAxisToRotaionMatrix( [1 0 0],RotationAngle )     ;
% % % RMat = RotationAxisToRotaionMatrix( [0 1 0],RotationAngle )   ;  
% % % RMat=inv(RMat);
% % 
% % RefXY = XY(18*k-17:18*k ,:) ; GcRefXY = mean(RefXY) ;  RefXY = RefXY -GcRefXY  ;   RefXY(:,3)=0 ;
% % RefXYcyl = XYcyl(18*k-17:18*k , :) ;  RefXYcyl = RefXYcyl -GcRefXY  ;  RefXYcyl(:,3)=0 ;
% % Gcenter = mean(ConfToTwist(Inds,:)) ;
% %  ConfToTwist(Inds,:) = ConfToTwist(Inds,:) -Gcenter ;
% %  
% % %    A3byNaLL=  transpose(TrajOri(:,1:3,frI )) ;  % which will be transformed.
% %    [regParams,~,~]=absor(RefXY',RefXYcyl');
% % %   A_prime = transpose(regParams.R*A3byNaLL + regParams.t  ); % Base centers
% % 
% % ConfToTwist(Inds,:) = transpose(regParams.R*ConfToTwist(Inds,:)' + regParams.t  ); % Base centers
% % 
% % 
% % %     ConfToTwist(Inds,:) =   (ConfToTwist(Inds,:))*RMat  ;
% % %         ConfToTwist(Inds,:) =   (ConfToTwist(Inds,:)- Gcenter)*RMat + Gcenter ;
% % % scatter3( ConfToTwist(Inds,1),ConfToTwist(Inds,2),ConfToTwist(Inds,3));
% % 
% % %   ConfToTwist(Inds,2)=  ConfToTwist(Inds,2)+100 ;
% %     BVecToTwist(Inds,:) =  transpose(regParams.R*BVecToTwist(Inds,:)'  ); 
% %     NVecToTwist(Inds,:) =  transpose(regParams.R*NVecToTwist(Inds,:)'   ); 
% %     
% % end
% for k = 1 : length(XSlice)-1
% Inds  =   and(ConfToTwist(:,1)>=XSlice(k) ,ConfToTwist(:,1)<XSlice(k+1)  ) ;
% RotationAngle = AngArr(k) ;
% RMat = RotationAxisToRotaionMatrix( [1 0 0],RotationAngle )     ;
% Gcenter = mean(ConfToTwist(Inds,:)) ;
%     ConfToTwist(Inds,:) =   (ConfToTwist(Inds,:)- Gcenter)*RMat + Gcenter ;
%     BVecToTwist(Inds,:) =   BVecToTwist(Inds,:)*RMat  ;
%     NVecToTwist(Inds,:) =   NVecToTwist(Inds,:)*RMat  ;
% end
% CentersVec(IndAccount ,:)  = ConfToTwist ;
% BVec(IndAccount ,:)=BVecToTwist ;
% NVec(IndAccount ,:)=NVecToTwist ;
% % % %-----------hard code

% CentersVec(:,1)=CentersVec(:,1)/7*5 ;
%------
file2_name='prova3322.conf';
fileID = fopen(file2_name,'w');
fprintf(fileID,'t = 0\n'); E0=0;
fprintf(fileID,'b = %9.6f %9.6f %9.6f\n',max(boxsize),max(boxsize),max(boxsize));
fprintf(fileID,'E = %8.6f %8.6f %8.6f\n',E0,E0,E0);

for k=1:  size(CentersVec,1)
    fprintf(fileID,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %2.1f %2.1f %2.1f %2.1f %2.1f %2.1f\n',CentersVec(k,1:3)/0.85,BVec(k,1:3),-NVec(k,1:3),zeros(1,6) );  % changed NVec for primes of oxDNA
end
fclose(fileID);

fprintf('finished printing oxdna formats\n'  )


% CentersVec2=CentersVec+0.8*BVec ;
%  figure(7788);cla; plot3(CentersVec2(:,1),CentersVec2(:,2),CentersVec2(:,3),'.-'  ) ;hold on;

% %   Coeff=0.8 ;delta=0.01;
% %   trials = 0.5:0.001:0.7;
% %   results=zeros(length(trials),2) ;
% %  for k=1:length(trials)
% %     Coeff= trials(k);
% %     CentersVec2=CentersVec+Coeff*BVec ;
% %     C = uniquetol(CentersVec2,0.001,'ByRows',true);
% % %
% % %    CentersVecRR=CentersVec+(Coeff+delta)*BVec ;
% % %    CRR = uniquetol(CentersVec2,0.0001,'ByRows',true); size(CRR)
% % %
% % %     CentersVecLL=CentersVec+(Coeff-delta)*BVec ;
% % %     CLL = uniquetol(CentersVec2,0.0001,'ByRows',true); size(CLL)
% %    results(k,:) = [trials(k), size(C,1)] ;
% %
% %  end
% %  figure(77882);clf;plot( results(:,1) ,results(:,2))
% %  Result show :[0.6,tol=0.001]
% %
% Coeff=0.57 ;  %0.6
% CentersVec2=CentersVec+Coeff*BVec ;   % used for mrDNA base position
%    figure;scatter3(CentersVec2(:,1) , CentersVec2(:,2),CentersVec2(:,3),'.')

DontUseConf = [ cell2mat( GetHyperB.ScafAllBase); cell2mat( GetHyperB.StapAllBase)] ;
% [~,IA,~] = uniquetol(DontUseConf,0.002,'ByRows',true,'OutputAllIndices',true );

[~,IA22,IC22] = unique(DontUseConf,'rows' );
IA = cell(size(IA22)) ;
for k =1: length(IA)
  IA{k} =    find(IA22(k) == IC22) ;
end
    
%  [C,IA,IC] = uniquetol(CentersVec2,0.002,'ByRows',true,'OutputAllIndices',true );
LL_IA= cellfun('length',IA) ;
ScafStapCorrelate2 =zeros(sum(LL_IA ==1) ,2 ); cj=1;
%   ScafStapCorrelate2 =zeros(sum(LL_IA ==2) ,2 ); cj=1;

Iamwroong=0;
for j=1:length(LL_IA)
    if LL_IA(j)==2
        ScafStapCorrelate2(cj,:) =(IA{j})' ; cj=cj+1;
        %          checkSeq =(IA{j})'
        %          [
        SeqPair=saveallseq((IA{j})');
        if strcmp(SeqPair,'AT')||strcmp(SeqPair,'TA')||strcmp(SeqPair,'CG')||strcmp(SeqPair,'GC')
%             sdf=3;
        else
            Iamwroong=Iamwroong+1 ;
            ScafStapCorrelate2(cj-1,:) =[0,0 ] ;
        end
        
    end
end
ScafStapCorrelate2( sum(ScafStapCorrelate2,2)==0 ,:)=[];

ScafStapCorrelate2=ScafStapCorrelate2-1 ; % index difference between Matlab and python

fprintf('printing deRemain\n')
file3_name='dSRemain.conf'   ;
fileID5 = fopen(file3_name,'w');

for iF=1:size(ScafStapCorrelate2,1)
    fprintf(fileID5,'{\n' );
    fprintf(fileID5,'type = mutual_trap\n' );
    fprintf(fileID5,'particle = %u\n' ,ScafStapCorrelate2(iF,1));
    fprintf(fileID5,'ref_particle  = %u\n' ,ScafStapCorrelate2(iF,2));
    fprintf(fileID5,'stiff = %u \n' ,0.2 );
    fprintf(fileID5,'r0 = 1.2 \n'  );
    fprintf(fileID5,'PBC = 1 \n'  ); % peirodic boundary condition, suggested from the forum
    fprintf(fileID5,'}\n' );
    fprintf(fileID5,'{\n' );
    fprintf(fileID5,'type = mutual_trap\n' );
    fprintf(fileID5,'particle = %u\n' ,ScafStapCorrelate2(iF,2));
    fprintf(fileID5,'ref_particle  = %u\n' ,ScafStapCorrelate2(iF,1));
    fprintf(fileID5,'stiff = %u \n' ,0.2 );
    fprintf(fileID5,'r0 = 1.2 \n' );
    fprintf(fileID5,'PBC = 1 \n'  );
    fprintf(fileID5,'}\n' );
end

fclose(fileID5);

GetHyperB.MonomerDsRemain  = ScafStapCorrelate2 ; % already substract 1 for oxDNA indexing

% % 
% % % % %----------MOE hard code
% Inds = CentersVec(:,2)<290.5  ;
% Indsub1 = and(CentersVec(:,1)<250 , Inds) ;
% Indsub2 = and(CentersVec(:,1)>380 , Inds) ;
% Indsub3 = and( and(CentersVec(:,1)<=380,CentersVec(:,1)>=250) , Inds) ;
% BM3(Indsub1) =2 ; BM3(Indsub2) =3 ;BM3(Indsub3) =4 ;
% fileTransInd_name=strcat( 'BM3.mat') ;
% save(fileTransInd_name, 'BM3') ;
% % %------------------



FineTuneOXDNA_Conf_wBM     ;

return
%---------------mrDNA Chrstopher Muffeo, UIUC
% header: x,y,z,bpindex ,stackindex, ntOn3'
%------------scaf
%--------------mrDNA_stack
mrDNA_stack=[ cell2mat( GetHyperB.ScafAllBase); cell2mat( GetHyperB.StapAllBase)] ; % C5 all-base routing, come from cadnanoInitial
scafsL =size(cell2mat( GetHyperB.ScafAllBase),1 ) ;
mrDNA_stack(:,end+1)=-1 ;

for k=1:size(mrDNA_stack,1)-1
    if mrDNA_stack(k,1)==mrDNA_stack(k+1,1) && abs(diff(mrDNA_stack(k:k+1,2) ))<=2 % due to skip within one base
        mrDNA_stack(k,3) =k;
    end
    
end
Inds=find(mrDNA_stack(:,3) ==-1) ; % check crossover
for k=1:length(Inds)
    C5_Base_Stack = mrDNA_stack(Inds(k),:)      ;
    if mod(GetHyperB.RelateVec(C5_Base_Stack(1)),2) ==0  && Inds(k)<=scafsL    % check helical direction: cylinder direction ,  isScaf
        TargetStack = C5_Base_Stack(1:2) + [0 ,1] ;
    elseif mod(GetHyperB.RelateVec(C5_Base_Stack(1)),2) ==1  && Inds(k)<=scafsL
        TargetStack = C5_Base_Stack(1:2) + [0 ,-1] ;
    elseif mod(GetHyperB.RelateVec(C5_Base_Stack(1)),2) ==0  && Inds(k)>scafsL
        TargetStack = C5_Base_Stack(1:2) + [0 ,-1] ;
    elseif mod(GetHyperB.RelateVec(C5_Base_Stack(1)),2) ==1  && Inds(k)>scafsL
        TargetStack = C5_Base_Stack(1:2) + [0 ,1] ;
    end
    [Lia,Locb] = ismember(TargetStack, mrDNA_stack(:,1:2),'rows' ) ;
    if  Lia==1
        QQ=find(and(mrDNA_stack(:,1)==TargetStack(1),mrDNA_stack(:,2)==TargetStack(2))) ; % should be two: one is scaffold base, the other is staple
        for cc=1:length(QQ)
            if ~xor(Inds(k)>scafsL ,QQ(cc)-1>scafsL) % find same domain one
                [Inds(k) ,  QQ(cc)];
                mrDNA_stack(Inds(k) ,3 ) = QQ(cc)-1 ;
            end
        end
        %           [Inds(k),  Locb-1]
    end
end
% %--------------mrDNA_stack



ScafStapCorrelate2_save =ScafStapCorrelate2 ;  % come from dsRemain for bpindex column
IndSiwtch =ScafStapCorrelate2(:,1)>  ScafStapCorrelate2(:,2) ;
ScafStapCorrelate2(IndSiwtch,:)=[ScafStapCorrelate2(IndSiwtch,2),ScafStapCorrelate2(IndSiwtch,1)];
bptable =sortrows(ScafStapCorrelate2) ;


fprintf('printing mrDNA \n')
file3_name='mrDNA.dat'   ;
fileID4 = fopen(file3_name,'w');
fprintf(fileID4,'x y z bpi sti nt3pi \n');   %header
s=0;
for mul_scafi= 1: length(scafL)
    
    [u,v]=find(bptable==s);   % find complementary base of one scaffold;
    if isempty(u)
        bpi= -1 ;
    elseif v==1
        bpi= bptable(u,2)  ;
    else
        bpi= bptable(u,1)  ;
    end
    %     [s , bpi]
    fprintf(fileID4,'%4.2f %4.2f %4.2f %i %i %i \n',  CentersVec2(1,1:3),bpi, mrDNA_stack(s+1,3) , s+1)  ; s=s+1 ;
    for scafi=2:scafL(mul_scafi)
        [u,v]=find(bptable==s);   % find complementary base ;
        if isempty(u)
            bpi= -1 ;
        elseif v==1
            bpi= bptable(u,2)  ;
        else
            bpi= bptable(u,1)  ;
        end
        %         [s,bpi]
        if scafi~=scafL(mul_scafi)
            fprintf(fileID4,'%4.2f %4.2f %4.2f %i %i %i \n',  CentersVec2(s+1,1:3),bpi, mrDNA_stack(s+1,3), s+1)  ; s=s+1 ;
        else
            fprintf(fileID4,'%4.2f %4.2f %4.2f %i %i %i \n',  CentersVec2(s+1,1:3),bpi, mrDNA_stack(s+1,3) , -1)  ; s=s+1 ;  % end of the scaffold
        end
        
    end
end
%------------staple
for stpi2=1:length(GetHyperB.StapList3)
    %     StapiSeq=tData{stpi2,3};    StpL
    [u,v]=find(bptable==s);
    if isempty(u)
        bpi= -1 ;
    elseif v==1
        bpi= bptable(u,2)  ;
    else
        bpi= bptable(u,1)  ;
    end
    fprintf(fileID4,'%4.2f %4.2f %4.2f %i %i %i \n',  CentersVec2(s+1,1:3),bpi, mrDNA_stack(s+1,3) , s+1)  ; s=s+1 ;
    
    for stapiline=2:StpL(stpi2)
        [u,v]=find(bptable==s);
        if isempty(u)
            bpi= -1 ;
        elseif v==1
            bpi= bptable(u,2)  ;
        else
            bpi= bptable(u,1)  ;
        end
        
        if stapiline~=StpL(stpi2)
            fprintf(fileID4,'%4.2f %4.2f %4.2f %i %i %i \n',  CentersVec2(s+1,1:3),bpi, mrDNA_stack(s+1,3) , s+1)  ; s=s+1 ;
        else
            fprintf(fileID4,'%4.2f %4.2f %4.2f %i %i %i \n',  CentersVec2(s+1,1:3),bpi, mrDNA_stack(s+1,3) , -1)  ; s=s+1 ;
        end
        
    end
end

fclose(fileID4); % printing  mrDNA file


%-

end




function silderexplode(src,evn, GetHyperB  ,scaf, stap ,CornerNotation , StapList3,CS ,ClosingCornerNotation)

pScaf2=scaf.pScaf2 ; ScafHelix=scaf.ScafHelix ; pScaf_center= scaf.pScaf_center ; ScafBaseCenterHelix= scaf.ScafBaseCenterHelix;

pStap=stap.pStap ; StapHelix=stap.StapHelix ; pStap_center= stap.pStap_center ; StapBaseCenterHelix= stap.StapBaseCenterHelix;


SacfR=GetHyperB.ScafRouting ;
% MaxBase=max(SacfR(:,3));
% skipPattern1=9:60:MaxBase;   % skip mod2 =0   %test 4/20
% skipPattern2=39:60:MaxBase;
skipBase= GetHyperB.skipBase ;

Val =src.Value ;

for stai=1:length(CornerNotation)   %actually mean scaffold, in case multi scaffolds in futures
    StapAll=CornerNotation{stai}    ;
    lineXYZ=ScafHelix{stai} ;
    scatterXYZ =  ScafBaseCenterHelix{stai} ;
    
    %      StapAll=GetHyperB.StapList3{stai}    ;
    nn= 1 ;
    for edgeinSCR=1:2:size(StapAll,1)
        C5Cyl=StapAll(edgeinSCR,1);
        bundle=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,1)   ;%multi-section needs be cautious
        bundle=bundle(1) ;
        Bundle=GetHyperB.containBundle{bundle};
        BaseStart=StapAll(edgeinSCR,2);   BaseEnd=StapAll(edgeinSCR+1,2);
        if  strcmp(Bundle.Lattice, 'Square')
            if BaseStart<BaseEnd  %go to left------------------stap
                %                 skipP=skipPattern1;
                skipP= skipBase(skipBase(:,1)==C5Cyl,2) ;
            else
                %                 skipP=skipPattern2;
                skipP= skipBase(skipBase(:,1)==C5Cyl,2) ;
            end
        else
            skipP=[];        %only square lattice needs skip
        end
        
        BasesArr =  setdiff(linspace(BaseStart,BaseEnd , abs(BaseStart-BaseEnd) +1  ) , skipP,'stable') ;
        nK =  length(BasesArr) ;
        
        %          if nn+nK-1>size(lineXYZ,1)
        %         sdsdf=3
        %          end
        PartHelix2 = lineXYZ(nn:nn+nK-1 ,:  ) ;
        PartScatter2 = scatterXYZ(nn:nn+nK-1 ,:  ) ;
        
        T_simulation=GetHyperB.containBundle{bundle}.SimulateTransMFromTM2;  % simulation matrix
        Vec = T_simulation(1:3,4)*Val ;
        
        QQQ= transpose(  eye(3)*( PartHelix2' -Vec*ones(1, size(PartHelix2,1)) )   );   % parallel position
        WWW= transpose(  eye(3)*( PartScatter2' -Vec*ones(1, size(PartScatter2,1)) )   );   % parallel position
        
        lineXYZ(nn:nn+nK-1 ,:  ) = QQQ;
        scatterXYZ(nn:nn+nK-1 ,:  ) = WWW;
        
        nn=nn+nK ;
    end
    
    pScaf2{stai}.XData =  lineXYZ(:,1)' ;pScaf2{stai}.YData =  lineXYZ(:,2)' ;pScaf2{stai}.ZData =  lineXYZ(:,3)' ;
    pScaf_center{stai}.XData =  scatterXYZ(:,1)' ;pScaf_center{stai}.YData =  scatterXYZ(:,2)' ;pScaf_center{stai}.ZData =  scatterXYZ(:,3)' ;
    
end
%% stap
nnK=0;
for stai=1:  length(StapList3)   % staples
    StapAll=StapList3{stai}    ;
    lineXYZ=StapHelix{stai} ;
    scatterXYZ =  StapBaseCenterHelix{stai} ;
    
    %      StapAll=GetHyperB.StapList3{stai}    ;
    nn= 1 ;
    for edgeinSCR=1:2:size(StapAll,1)
        C5Cyl=StapAll(edgeinSCR,1);
        bundle=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,1)   ;%multi-section needs be cautious
        bundle=bundle(1) ;
        Bundle=GetHyperB.containBundle{bundle};
        Cyl=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,2);   %C2 express
        Cyl=Cyl(1) ;
        BaseStart=StapAll(edgeinSCR,2);   BaseEnd=StapAll(edgeinSCR+1,2);
        if  strcmp(Bundle.Lattice, 'Square')  && Cyl~=-1   % not on overhang
            if BaseStart<BaseEnd  %go to left------------------stap
                %                 skipP=skipPattern2;
                skipP= skipBase(skipBase(:,1)==C5Cyl,2) ;
            else
                %                 skipP=skipPattern1;
                skipP= skipBase(skipBase(:,1)==C5Cyl,2) ;
            end
        else
            skipP=[];        %only square lattice needs skip
        end
        
        BasesArr =  setdiff(linspace(BaseStart,BaseEnd , abs(BaseStart-BaseEnd) +1  ) , skipP,'stable')    ;
        nK =  length(BasesArr) ;
        
        PartHelix2 = lineXYZ(nn:nn+nK-1 ,:  ) ;
        PartScatter2 = scatterXYZ(nn:nn+nK-1 ,:  ) ;
        
        T_simulation=GetHyperB.containBundle{bundle}.SimulateTransMFromTM2;  % simulation matrix
        Vec = T_simulation(1:3,4)*Val ;
        
        QQQ= transpose(  eye(3)*( PartHelix2' -Vec*ones(1, size(PartHelix2,1)) )   );   % parallel position
        WWW= transpose(  eye(3)*( PartScatter2' -Vec*ones(1, size(PartScatter2,1)) )   );   % parallel position
        
        lineXYZ(nn:nn+nK-1 ,:  ) = QQQ;
        scatterXYZ(nn:nn+nK-1 ,:  ) = WWW;
        
        nn=nn+nK ;
    end
    nnK=nnK+nn-1 ;
    %             fprintf( 'sti= %i ,  nnK= %i \n',   stai,nnK) ;
    
    
    pStap{stai}.XData =  lineXYZ(:,1)' ;pStap{stai}.YData =  lineXYZ(:,2)' ;pStap{stai}.ZData =  lineXYZ(:,3)' ;
    pStap_center{stai}.XData =  scatterXYZ(:,1)' ;pStap_center{stai}.YData =  scatterXYZ(:,2)' ;pStap_center{stai}.ZData =  scatterXYZ(:,3)' ;
    
end
% nnK
%%  closing strand
% CS.pCS=pCS ; CS.CShelix=CShelix ; CS.pCS_center=pCS_center ; CS.CSBaseCenterHelix=CSBaseCenterHelix ; CS.NVecCS=NVecCS ;
% ClosingCornerNotation ;


for csi=1:length(ClosingCornerNotation)   % staples
    CsAll=ClosingCornerNotation{csi}    ;
    lineXYZ=CS.CShelix{csi} ;
    scatterXYZ =  CS.CSBaseCenterHelix{csi} ;
    
    %      CsAll=GetHyperB.StapList3{csi}    ;
    nn= 1 ;
    for edgeinSCR=1:2:size(CsAll,1)
        C5Cyl=CsAll(edgeinSCR,1);
        bundle=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,1)   ;%multi-section needs be cautious
        Bundle=GetHyperB.containBundle{bundle};
        Cyl=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,2);   %C2 express
        
        BaseStart=CsAll(edgeinSCR,2);   BaseEnd=CsAll(edgeinSCR+1,2);
        if  strcmp(Bundle.Lattice, 'Square')  && Cyl~=-1   % not on overhang
            if BaseStart<BaseEnd  %go to left------------------stap
                skipP=skipPattern1;
            else
                skipP=skipPattern2;
            end
        else
            skipP=[];        %only square lattice needs skip
        end
        
        BasesArr =  setdiff(linspace(BaseStart,BaseEnd , abs(BaseStart-BaseEnd) +1  ) , skipP,'stable') ;
        nK =  length(BasesArr) ;
        
        PartHelix2 = lineXYZ(nn:nn+nK-1 ,:  ) ;
        PartScatter2 = scatterXYZ(nn:nn+nK-1 ,:  ) ;
        
        T_simulation=GetHyperB.containBundle{bundle}.SimulateTransMFromTM2;  % simulation matrix
        Vec = T_simulation(1:3,4)*Val ;
        
        QQQ= transpose(  eye(3)*( PartHelix2' -Vec*ones(1, size(PartHelix2,1)) )   );   % parallel position
        WWW= transpose(  eye(3)*( PartScatter2' -Vec*ones(1, size(PartScatter2,1)) )   );   % parallel position
        
        lineXYZ(nn:nn+nK-1 ,:  ) = QQQ;
        scatterXYZ(nn:nn+nK-1 ,:  ) = WWW;
        
        nn=nn+nK ;
    end
    
    CS.pCS{csi}.XData =  lineXYZ(:,1)' ;CS.pCS{csi}.YData =  lineXYZ(:,2)' ;CS.pCS{csi}.ZData =  lineXYZ(:,3)' ;
    CS.pCS_center{csi}.XData =  scatterXYZ(:,1)' ;CS.pCS_center{csi}.YData =  scatterXYZ(:,2)' ;CS.pCS_center{csi}.ZData =  scatterXYZ(:,3)' ;
    
end





end

