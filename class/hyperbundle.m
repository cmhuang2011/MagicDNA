classdef hyperbundle <handle
    %hybrid for SQ and HC
    % <MagicDNA (Multi-component Assembly in a Graphical Interface guided by Computation for DNA origami) is a software for designing multi-component DNA origami structures.>
%     Copyright (C) <2020>  <Chao-Min Huang, Hai-Jun Su, and Carlos E. Castro>
%     The Ohio State University 
%     
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

    %   Detailed explanation goes here
    
    properties
        
        containBundle=[];
        BundleHasCylinder=[];  % ex: [24 12] -> first bundles has 24 cylinders, 2nd bundle has 12
        GlobalCylinderIndex=[];  % [whichBundles,  localcylinderInThis bundle]
        CylinderC5AdjM=[];
        
        ForcedConnectList=[];
        FCListRelativeO=[];
        
        BasicCycleList=[];
        
        ScafRouting=[];  %index use RTable column [1 2]
        ScafGlobal=[];  %index use RTable column 3
        scafC4=[];   %index use RTable column 4
        scafC5=[];    %index use RTable column 5
        ScafdigitSQ=[]
        ScafdigitHC=[]
        ScafAllBase=[];   % 08/12/2019, multi-scaf
        AllScafXover=[];    % Nov 13 2018 , for zigzag scaf routing
        
        
        ScafXover=[];
        
        RelateTableOrr=[];
        
        RelateTable=[];
        RelateVec=[];
        OverhangRT=[];
        uniqueInplaneXY=[];
        ClosingStrand=[];     %for overhang sequencing
        ClosingCornerNotation ; 
        MonomerDsRemain = []; % already substract 1 for oxDNA indexing
        
        
        StapList=[];  %use C5
        DigitStapSQ=[];
        DigitStapHC=[];
        
        stapBP=[];
        OstapBP=[];
        stapBPinCell=[];
        StapList2=[];  %use C5
        StapList3=[];  %use C5
        SaveStapList3=[];  %use C5
        HeadOfStap=[];
        StapAllBase=[]; %Sep 26, mrDNA
        StapGraphBin=[]; % Sep 24 2019, show staple graph, for coloring chain, Update by FindStapGraph
        
        UConPlotHandle={};
        choice=[];
        TakeSeveralV3=[]; TakeSeveralV3Type=[]; % auto or manual
        LineCombination=[];  LineOption=[];
        OverRoute=2;
        ssOption=[];
        
        ScafOption=[];
        StapOption=[];
        
        Scaf_fromCadDOM = [] ; % sacffold routing generated from the algrorithm
        Scaf_fromJSON = []; % scafould routing from cadnano json file
        Stap_fromCadDOM = [] ; % staple routing generated from the algrorithm
        Stap_fromJSON = []; % staple routing from cadnano json file
        
        SavePremPair; SavePremPairManual ;
        BundleAdjM;
        %--------Multi-Scaf 08/12/2019
        SaveBlockAdjacency ; % for multi-scaffold , cell arraay(N=#ofBundle) contains undesired cylinder pairs in bundles
        SaveBlockZRange ;  % Use Edit Bundle to assgin where in a bundle without Xovers.
        OrgBlockAdj;  % [B, C1 ,B ,C2] get property once to avoid querry in loops
        SaveInternalXovers ;  % Use Edit Bundle to assgin internal Xovers as forced connections.
        
        
        CompJointInfo=[];
        UserWantOH=[];
        pSeq ; pSeqExt; ScafUnUsed; ScafUnUsedSeq; SeqsOfClosing;
        
        
        TransformHistroy=[];
        Rcenter_scatterH = [];
        ObjFigure=[];
        %------------
        CustomSkipAndInsertion = false ;  % default == false; true->use caDNAno
        DomainGraph =[];
        TMs;  % HA result
    end
    
    properties (Dependent)
        ForcedConnectListC4;   %express in Column 4
        ForcedConnectListC5;   %express in Column 5
        skipBase ;
        insertBase ; 
        
        scafC5All ;     % pseudo scaffold, put multiple scaffold together
        pSeqAll ;
        
        DuplexDomain =[];  %MagicDNA2 
        DuplexDomain_neglectNick =[];  %MagicDNA2 
    end
    
    methods
        
        function getDomainGraph(obj)
            Domains=obj.DuplexDomain;
            AdjList = zeros(10000,2) ; c=1 ;
            count=1 ;
            %-------by staple Xovers
            for k_scaf = 1: length(  obj.ScafAllBase ) 
                for j_stap = 1: length(obj.StapAllBase)
                    C5_scaf = obj.ScafAllBase{k_scaf} ;    C5_staple = obj.StapAllBase{j_stap} ;
                    [~,IndexOnScaf ]=ismember(C5_staple , C5_scaf,'rows') ;
                    IndOverhang= find(IndexOnScaf==0) ;
                    IndexOnScaf(IndOverhang) = -1:-1:-length(find(IndexOnScaf==0)) ;  % consider overhangs as continuos
                    MM =  [IndexOnScaf,[0;cumsum(diff(IndexOnScaf)~=-1)]]     ;       
                    
                    SubDomains = unique(MM(:,2)) ;
                    GlobalInd=  SubDomains+count  ;
                    NewList = [GlobalInd(1:end-1) ,GlobalInd(2:end) ] ;
                    if ~isempty(NewList) %in case of single domain in one staple
                    AdjList(c:c+size(NewList,1)-1 ,:) =  NewList ;
                    c=c+size(NewList,1);                  
                    end
                                        
                    for seg =0:max(MM(:,2))
%                        IndSCaff=  MM(:,2) ==seg  ;
%                        C5_duplex =  C5_scaf(MM(IndSCaff,1) ,:) ;
                       count=count+1 ;
                    end
                end
            end
            AdjList=AdjList(1:c-1 ,:) ; 
            %----------
            NewList = zeros(10000,2) ; c=1 ;  % along cylinder direction, including scaffold and staple
              for j_d = 1:length(Domains)
                 for j_d2 = j_d+1:length(Domains)
                               
                    if Domains{j_d}(1,1)==Domains{j_d2}(end,1) && abs(Domains{j_d}(1,2)-Domains{j_d2}(end,2))==1
                        NewList(c,:) = [j_d,j_d2 ]; c=c+1 ;
                    end
                     if Domains{j_d2}(1,1)==Domains{j_d}(end,1) && abs(Domains{j_d2}(1,2)-Domains{j_d}(end,2))==1
                        NewList(c,:) = [j_d,j_d2 ]; c=c+1 ;
                    end                                      
                 end
              end
             NewList=NewList(1:c-1 ,:) ; 
             %-------------Scaffold Xover ;
             NewList2 = zeros(10000,2) ; cN2=1 ;  
          
             [AA,~] =cellfun(@size,obj.ScafAllBase );
                GetScaf= zeros(sum(AA)  ,2) ; c=1 ;
             for sc=1:length(  obj.ScafAllBase ) 
                GetScaf(c:c+size( obj.ScafAllBase{sc},1)-1 ,:) =obj.ScafAllBase{sc} ;
                c=c+size( obj.ScafAllBase{sc},1) ;
             end
            Df= GetScaf(2:end,:) -GetScaf(1:end-1,:) ;
            Df=find(Df(:,2)==0) ;
            Df=union(Df,Df+ones(size(Df))) ;
            XoverList = GetScaf(Df,:)  ;
            XoverList(:,end+1 ) =0 ;
            for j_d = 1:length(Domains)
                [aa,~ ]=ismember(XoverList(:,1:2),Domains{j_d} , 'rows') ;
                if  sum(aa)>0
                    XoverList(aa,3) =j_d ;
                end
            end
            for k= 1:2:size(XoverList ,1)
                if XoverList(k,3)~=0 &&  XoverList(k+1,3)~=0
                    NewList2(cN2,:) = XoverList(k:k+1,3)' ; cN2=cN2+1;
                end
            end
             NewList2=NewList2(1:cN2-1 ,:) ;    
             IndF = NewList2(:,1)>  NewList2(:,2) ;
             NewList2(IndF,:) = flip(    NewList2(IndF,:) ,2) ;
             %------------------
         TotalList = [AdjList;NewList ;NewList2] ;
         TotalList=unique(TotalList ,'rows') ;
         
          obj.DomainGraph  =graph(TotalList(:,1) ,TotalList(:,2)) ;
          figure ; plot(  obj.DomainGraph ) ;
        end
        function temp_shiftScaf(obj)
            QQ = GetHyperB.Scaf_fromCadDOM{1} ;
            [a,b] = ismember([2 1 68 ; 2 1 69 ] , QQ ,'rows' ) ;
            prev = QQ(b(1)-1 ,:) ;
            QQ(b(2) ,:) = prev ;
            QQ(b(1)-1:b(1),:) = [];            
        end
        function AllDomains=get.DuplexDomain(obj)
            obj.ScafAllBase ;
            obj.StapAllBase  ;
            count=1 ;
            for k_scaf = 1: length(  obj.ScafAllBase ) 
                for j_stap = 1: length(obj.StapAllBase)
                    C5_scaf = obj.ScafAllBase{k_scaf} ;    C5_staple = obj.StapAllBase{j_stap} ;
                    [aa,IndexOnScaf ]=ismember(C5_staple , C5_scaf,'rows') ;
                    IndOverhang= find(IndexOnScaf==0) ;
                    IndexOnScaf(IndOverhang) = -1:-1:-length(find(IndexOnScaf==0)) ;  % consider overhangs as continuos
                    
                    Continue= abs(diff(IndexOnScaf)) >5 ;   %in case of insertion
                     MM =  [IndexOnScaf,[0;cumsum(Continue)] ]     ;                    
                   
                    for seg =0:max(MM(:,2))
%                        IndSCaff=  MM(:,2) ==seg  ;
%                        C5_duplex =  C5_scaf(MM(IndSCaff,1) ,:) ;
                       count=count+1 ;
                    end
                end
            end
            %---------
            AllDomains = cell(count-1 ,1 ) ; c=1 ;
            DomainLength = zeros(count-1 ,1 ) ;
            for k_scaf = 1: length(  obj.ScafAllBase ) 
                for j_stap = 1: length(obj.StapAllBase)
                    C5_scaf = obj.ScafAllBase{k_scaf} ;    C5_staple = obj.StapAllBase{j_stap} ;
                    [aa,IndexOnScaf ]=ismember(C5_staple , C5_scaf,'rows') ;
                    IndOverhang= find(IndexOnScaf==0) ;
                    IndexOnScaf(IndOverhang) = -1:-1:-length(find(IndexOnScaf==0)) ;  % consider overhangs as continuos
                    
                    Continue= abs(diff(IndexOnScaf)) >5 ;   %in case of insertion
                     MM =  [IndexOnScaf,[0;cumsum(Continue)] ]     ;                    
                    for seg =0:max(MM(:,2))
                       IndSCaff=  MM(:,2) ==seg  ;
                       C5_duplex =  C5_scaf(MM(IndSCaff,1) ,:) ;
                        AllDomains{c} =C5_duplex ; 
                        DomainLength(c)= size(C5_duplex,1) ; 
                        c=c+1 ;
                    end
                end
            end
            %--------
%             sdfdsf=3
%             Domains=[] ;
        end
        function [AllDomains_negNick, MergeTableInOri]=DuplexDomain_neglectNickFcn(obj)
           % Merge Domains created by nicks .
            %--------
            AllDomains=obj.DuplexDomain ;
            DesignSkipBase =  obj.skipBase ;
            
             stapC5 =obj.StapAllBase ;
%              stapC5Mat=cell2mat(stapC5) ;
             StapHeadTail = zeros( 2*length(stapC5),2) ;
             for k=1:length(stapC5)
                 ConsiderSkip =  setdiff(stapC5{k},DesignSkipBase,'rows','stable');
              StapHeadTail(2*k-1:2*k,:) =  ConsiderSkip([1 end],:);
             end
%              obj.skipBase
             
             DomainHeadTail= zeros( 2*length(AllDomains),2) ;
             for k=1:length(AllDomains)
              DomainHeadTail(2*k-1:2*k,:) =    AllDomains{k}([1 end],:);
             end             
             
             Adj_forDomain = repmat([-1 -1],length(stapC5),1   ) ; c=1;
             for stapHT = 1:length(stapC5)   %self-nick
                  [~, indNick_forDomain] = ismember(StapHeadTail(2*stapHT-1:2*stapHT,:) ,DomainHeadTail,'rows') ;
                  indNick_forDomain=round(indNick_forDomain/2 ) ;

%                   stapHT
%                   if ~all(indNick_forDomain);continue;end
                   MM = DomainHeadTail([2*indNick_forDomain-1 2*indNick_forDomain],:) ;
%                    sMM2 =sort(MM(:,2)) ;
                  if abs(MM(1,2)-MM(4,2))<=1 &&  MM(1,1)==MM(4,1)
                   Adj_forDomain(c,:) =indNick_forDomain(:) ; c=c+1 ;
                  end
             end
             
              for stapi = 1:length(stapC5)   % different staples 
                for stapj = stapi+1:length(stapC5) 
                   HeadTail1 = StapHeadTail([2*stapi-1, 2*stapj],:)    ;                  
                  [~, indNick_forDomain] = ismember(HeadTail1 ,DomainHeadTail,'rows') ;
                  indNick_forDomain=round(indNick_forDomain/2 ) ;
%                    if ~all(indNick_forDomain);continue;end
                   
                   MM = DomainHeadTail([2*indNick_forDomain-1 2*indNick_forDomain],:) ;
                  if abs(MM(1,2)-MM(4,2))<=2 &&  MM(1,1)==MM(4,1) %abs(MM(3,2)-MM(2,2))==1 
                   Adj_forDomain(c,:) =indNick_forDomain(:) ; c=c+1 ;
                  end                    
                    %---------
                  HeadTail2 = StapHeadTail([2*stapj-1, 2*stapi],:)    ;                  
                  [~, indNick_forDomain] = ismember(HeadTail2 ,DomainHeadTail,'rows') ;
                  indNick_forDomain=round(indNick_forDomain/2 ) ;
                   MM = DomainHeadTail([2*indNick_forDomain-1 2*indNick_forDomain],:) ;
                  if abs(MM(1,2)-MM(4,2))<=2 &&  MM(1,1)==MM(4,1)  %abs(MM(3,2)-MM(2,2))==1 
                   Adj_forDomain(c,:) =indNick_forDomain(:) ; c=c+1 ;
                  end                       
                end
              end
             
             Adj_forDomain=Adj_forDomain(1:c-1,:) ;
             Adj_forDomain=unique(Adj_forDomain,'rows','stable') ;
             
             
             AllDomains_negNick= AllDomains ;
             for k=1: size(Adj_forDomain,1)
                 AllDomains_negNick{Adj_forDomain(k,1)} = [ AllDomains{Adj_forDomain(k,2) };AllDomains{Adj_forDomain(k,1) }];
                 AllDomains_negNick{Adj_forDomain(k,2)}=[] ;                
             end
             IndRemain = cellfun(@isempty,AllDomains_negNick) ;
             AllDomains_negNick=AllDomains_negNick(~IndRemain) ;
             MergeTableInOri=Adj_forDomain ;
%         MergeTableInOri
%         size(MergeTableInOri)
%         for k =1:size(MergeTableInOri ,1)
%             AllDomains{MergeTableInOri(k,1)}
%             AllDomains{MergeTableInOri(k,2)}
%         end
%         
        
%             Domains=[] ;
        end
        function ExportCustomExtradsRemain(obj, N ,PC_mode , CustomCell)
            length( obj.ClosingCornerNotation)
            if length(CustomCell )~= length( obj.ClosingCornerNotation)
                fprintf('ERROR!! Input Custom Table dimension doesn''t match with the number of CS. Invalid!! \n')
                return ;                
            end
            
            for k =1: length(CustomCell)
                if max( CustomCell{k}) > N
                    fprintf('ERROR!! Specifying connections on non-exist instances. Invalid!! \n')
                    return ;
                end
            end
            
            StackScaf = cell2mat( obj.ScafAllBase ) ;
            StackStap = cell2mat( obj.StapAllBase ) ;
            StackAll = [ StackScaf ; StackStap ];
            
            cc= 0;
            for k_cs = 1 : length( obj.ClosingCornerNotation)
                CloseStrandAllBase = interpolateBase(obj.ClosingCornerNotation{k_cs} ) ;
                cc=cc+  size(CloseStrandAllBase,1)/2 ;
            end
            
            OverhangStaplePair= zeros(cc ,3 ) ;
            k=1;
            for k_cs = 1 : length( obj.ClosingCornerNotation)
                CloseStrandAllBase = interpolateBase(obj.ClosingCornerNotation{k_cs} ) ;               
                for li_cs = 1: size(CloseStrandAllBase,1)/2
                    [~,bb]  = ismember( CloseStrandAllBase(li_cs,:) , StackAll,'rows' ) ;
                    [~,dd]  = ismember( CloseStrandAllBase(end-li_cs+1,:) , StackAll,'rows' ) ;
                    OverhangStaplePair(k ,: ) = [bb,dd , k_cs] ; k=k+1 ;
                end
%                  hL = size(CloseStrandAllBase,1)/2;
                 RootCS = obj.ClosingStrand.RootAndExtC5(:, [1 3]) ;
                 IsHeadAsRoot = ismember(CloseStrandAllBase(1,:) , RootCS  ,'rows' );
                 IsTailAsRoot = ismember(CloseStrandAllBase(end,:) , RootCS  ,'rows' );
                 ssTL=obj.ClosingStrand.ssT_shift(k_cs);
                 
                 if xor(IsHeadAsRoot , IsTailAsRoot)  % assume as double overhang case1
                     if ~IsHeadAsRoot  % single stranded T, without mutual trap
                         OverhangStaplePair(k-ssTL:k-1 ,1) = nan;
                         OverhangStaplePair(k-li_cs:k-1 ,1)= circshift (OverhangStaplePair(k-li_cs:k-1 ,1) , ssTL);
                     else
                         OverhangStaplePair(k-ssTL:k-1 ,2) = nan;
                         OverhangStaplePair(k-li_cs:k-1 ,2)= circshift (OverhangStaplePair(k-li_cs:k-1 ,2) , ssTL);
                     end
                 else % assume as double overhang case2 
                     
                     if ~IsHeadAsRoot
                         OverhangStaplePair(k-ssTL:k-1 ,1) =nan ;
                     else
                        OverhangStaplePair(k-li_cs:k-li_cs+ssTL-1 ,1) = nan; 
                         
                     end
                     
                     
%                      if IsHeadAsRoot %Insert 'T' in both ends which are roots
%                          answer{1}(1:ssTL)='A' ;  % force these staple bases to T (CS='A') ;
%                          answer{1}(end-ssTL+1:end)='A' ;
%                      else  %Insert 'T' in the middle,  which are roots
%                          answer{1}(hL-ssTL+1:hL)='A' ;  % force these staple bases to T (CS='A') ;
%                          answer{1}(hL+1:hL+ssTL)='A' ;
%                      end
                 end
            end
            OverhangStaplePair_Ori = OverhangStaplePair ; % pair indices in monomer
            OverhangStaplePair = OverhangStaplePair(:,1:2) ;
            [u,~] = find(OverhangStaplePair ==0) ;
            OverhangStaplePair(u,:) = [];    
            
            TotalBase = size(StackAll,1) ;  % a.k.a, period over monomer
            
            MutualPairForPattern = zeros(10000,2 ) ;cM =1;
            for k =1: length(CustomCell)
               AcrossMonomer_CS =   CustomCell{ k} ; % l 
               ForThisCS = OverhangStaplePair_Ori(OverhangStaplePair_Ori(:,3)==k ,1:2) ;
               
               for j = 1 : size(AcrossMonomer_CS ,1 )
                   BtwMonomers = AcrossMonomer_CS(j, :) ;
%                    if BtwMonomers(2)<BtwMonomers(1)
%                        BtwMonomers= flip(BtwMonomers,2) ;
%                    end
                    BaseIndAcrMonomer =  repmat( TotalBase*(BtwMonomers-1) , size(ForThisCS,1) ,1  )  + ForThisCS;
%                     BaseIndAcrMonomer =  repmat( TotalBase*(BtwMonomers-1) , size(OverhangStaplePair,1) ,1  )  + OverhangStaplePair;
                   
                    
                    MutualPairForPattern(cM:cM+size(BaseIndAcrMonomer,1)-1 ,:) = BaseIndAcrMonomer ;
                    cM = cM + size(BaseIndAcrMonomer,1) ;
%                    sdfsf=2
               end
            end
            
            MutualPairForPattern=MutualPairForPattern(1:cM-1 ,:) ;
            
            
            PatterN_Extra_dsRemain = MutualPairForPattern -1 ; % oxDNA index, startiing as 0
            UU = isnan(PatterN_Extra_dsRemain) ;
            RowNan = or(UU(:,1) ,UU(:,2)) ;
            PatterN_Extra_dsRemain(RowNan ,: ) =[];

%             return ;%---------------------
            %-----------'dS_xtra.conf'
            file3_name='dS_xtra.conf'   ;
            fileID2 = fopen(file3_name,'w');           
            for iF=1:size(PatterN_Extra_dsRemain,1)
                fprintf(fileID2,'{\n' );
                fprintf(fileID2,'type = mutual_trap\n' );
                fprintf(fileID2,'particle = %u\n' ,PatterN_Extra_dsRemain(iF,1));
                fprintf(fileID2,'ref_particle  = %u\n' ,PatterN_Extra_dsRemain(iF,2));
                fprintf(fileID2,'stiff = %u \n' ,0.5 );
                fprintf(fileID2,'r0 = 1.2 \n'  );
                fprintf(fileID2,'PBC = 1 \n'  );
                fprintf(fileID2,'}\n' );
                fprintf(fileID2,'{\n' );
                fprintf(fileID2,'type = mutual_trap\n' );
                fprintf(fileID2,'particle = %u\n' ,PatterN_Extra_dsRemain(iF,2));
                fprintf(fileID2,'ref_particle  = %u\n' ,PatterN_Extra_dsRemain(iF,1));
                fprintf(fileID2,'stiff = %u \n' ,0.5 );
                fprintf(fileID2,'r0 = 1.2 \n' );
                fprintf(fileID2,'PBC = 1 \n'  );
                fprintf(fileID2,'}\n' );
            end
            fclose(fileID2);            
            %-----------
            MonoDsT =  obj.MonomerDsRemain;
            MultimerDsT = repmat(MonoDsT , N,1) ;
            N_DsT = size(MonoDsT ,1) ;
            AcrossfMonomer_DST= TotalBase*[ repelem((1:N)' ,N_DsT,1)-1 , repelem((1:N)' ,N_DsT,1)-1];
            MultimerDsT_final =MultimerDsT +  AcrossfMonomer_DST ;

            for jj= 1: 2
                if jj==1
                file4_name='dSRemain_wOH.conf'   ; Stiff_inOH=0.5 ; 
                else
                file4_name='dSRemain_wOH_p1.conf'   ;Stiff_inOH=0.05 ; 
                end
                
                fileID4 = fopen(file4_name,'w');
                for iF=1:size(PatterN_Extra_dsRemain,1)
                    fprintf(fileID4,'{\n' );
                    fprintf(fileID4,'type = mutual_trap\n' );
                    fprintf(fileID4,'particle = %u\n' ,PatterN_Extra_dsRemain(iF,1));
                    fprintf(fileID4,'ref_particle  = %u\n' ,PatterN_Extra_dsRemain(iF,2));
                    fprintf(fileID4,'stiff = %u \n' ,Stiff_inOH );
                    fprintf(fileID4,'r0 = 1.2 \n'  );
                    fprintf(fileID4,'PBC = 1 \n'  );
                    fprintf(fileID4,'}\n' );
                    fprintf(fileID4,'{\n' );
                    fprintf(fileID4,'type = mutual_trap\n' );
                    fprintf(fileID4,'particle = %u\n' ,PatterN_Extra_dsRemain(iF,2));
                    fprintf(fileID4,'ref_particle  = %u\n' ,PatterN_Extra_dsRemain(iF,1));
                    fprintf(fileID4,'stiff = %u \n' ,Stiff_inOH );
                    fprintf(fileID4,'r0 = 1.2 \n' );
                    fprintf(fileID4,'PBC = 1 \n'  );
                    fprintf(fileID4,'}\n' );
                end
                for iF=1:size(MultimerDsT_final,1)
                    fprintf(fileID4,'{\n' );
                    fprintf(fileID4,'type = mutual_trap\n' );
                    fprintf(fileID4,'particle = %u\n' ,MultimerDsT_final(iF,1));
                    fprintf(fileID4,'ref_particle  = %u\n' ,MultimerDsT_final(iF,2));
                    fprintf(fileID4,'stiff = %u \n' ,0.2 );
                    fprintf(fileID4,'r0 = 1.2 \n'  );
                    fprintf(fileID4,'PBC = 1 \n'  );
                    fprintf(fileID4,'}\n' );
                    fprintf(fileID4,'{\n' );
                    fprintf(fileID4,'type = mutual_trap\n' );
                    fprintf(fileID4,'particle = %u\n' ,MultimerDsT_final(iF,2));
                    fprintf(fileID4,'ref_particle  = %u\n' ,MultimerDsT_final(iF,1));
                    fprintf(fileID4,'stiff = %u \n' ,0.2 );
                    fprintf(fileID4,'r0 = 1.2 \n' );
                    fprintf(fileID4,'PBC = 1 \n'  );
                    fprintf(fileID4,'}\n' );
                end
                fclose(fileID4);                
            end
            
            
        end % end of ExportCustomExtradsRemain
        
        function ExportExtradsRemain(obj, N ,PC_mode)
            StackScaf = cell2mat( obj.ScafAllBase ) ;
            StackStap = cell2mat( obj.StapAllBase ) ;
            StackAll = [ StackScaf ; StackStap ];
            
            cc= 0;
            for k_cs = 1 : length( obj.ClosingCornerNotation)
                CloseStrandAllBase = interpolateBase(obj.ClosingCornerNotation{k_cs} ) ;
                cc=cc+  size(CloseStrandAllBase,1)/2 ;
            end
            
            OverhangStaplePair= zeros(cc ,2 ) ;
            k=1;
            % if use new version of single-stranded class,
            for k_cs = 1 : 0 %length( obj.ClosingCornerNotation)
                CloseStrandAllBase = interpolateBase(obj.ClosingCornerNotation{k_cs} ) ;
                
                for li_cs = 1: size(CloseStrandAllBase,1)/2
                    [~,bb]  = ismember( CloseStrandAllBase(li_cs,:) , StackAll,'rows' ) ;
                    [~,dd]  = ismember( CloseStrandAllBase(end-li_cs+1,:) , StackAll,'rows' ) ;
                    OverhangStaplePair(k ,:) = [bb,dd] ; k=k+1 ;
                end
%                  hL = size(CloseStrandAllBase,1)/2;
                 RootCS = obj.ClosingStrand.RootAndExtC5(:, [1 3]) ;
                 IsHeadAsRoot = ismember(CloseStrandAllBase(1,:) , RootCS  ,'rows' );
                 ssTL=obj.ClosingStrand.ssT_shift(k_cs);
                 if ~IsHeadAsRoot  % single stranded T, without mutual trap
                      OverhangStaplePair(k-ssTL:k-1 ,1) = 0;
                      OverhangStaplePair(k-li_cs:k-1 ,1)= circshift (OverhangStaplePair(k-li_cs:k-1 ,1) , ssTL);                     
                 else
                       OverhangStaplePair(k-ssTL:k-1 ,2) = 0;
                      OverhangStaplePair(k-li_cs:k-1 ,2)= circshift (OverhangStaplePair(k-li_cs:k-1 ,2) , ssTL);    
                 end
            end
            [u,~] = find(OverhangStaplePair ==0) ;
            OverhangStaplePair(u,:) = [];    
            
            TotalBase = size(StackAll,1) ;
            dsTable = repmat(OverhangStaplePair , N,1) ;
            NOSP = size(OverhangStaplePair ,1) ;
            
%             AcrossfMonomer= TotalBase*[ repelem((1:N)' ,NOSP,1)-1 , repelem((1:N)' ,NOSP,1)];
%             PatterN_Extra_dsRemain = dsTable + AcrossfMonomer ;
%             PatterN_Extra_dsRemain=PatterN_Extra_dsRemain-1 ;
%             
%             if strcmp( PC_mode ,'Chain')
%                 PatterN_Extra_dsRemain(end-size(OverhangStaplePair,1)+1:end ,:) =[];
%             elseif strcmp( PC_mode ,'Circular')
%                 PatterN_Extra_dsRemain(end-size(OverhangStaplePair,1)+1:end ,2)  = PatterN_Extra_dsRemain(end-size(OverhangStaplePair,1)+1:end ,2)   - TotalBase*N ;
%             else
%                 ShouldNotHappend =1231 ;
%             end
%             
%             file3_name='dS_xtra.conf'   ;
%             fileID2 = fopen(file3_name,'w');           
%             for iF=1:size(PatterN_Extra_dsRemain,1)
%                 fprintf(fileID2,'{\n' );
%                 fprintf(fileID2,'type = mutual_trap\n' );
%                 fprintf(fileID2,'particle = %u\n' ,PatterN_Extra_dsRemain(iF,1));
%                 fprintf(fileID2,'ref_particle  = %u\n' ,PatterN_Extra_dsRemain(iF,2));
%                 fprintf(fileID2,'stiff = %u \n' ,20 );
%                 fprintf(fileID2,'r0 = 1.2 \n'  );
%                 fprintf(fileID2,'PBC = 1 \n'  );
%                 fprintf(fileID2,'}\n' );
%                 fprintf(fileID2,'{\n' );
%                 fprintf(fileID2,'type = mutual_trap\n' );
%                 fprintf(fileID2,'particle = %u\n' ,PatterN_Extra_dsRemain(iF,2));
%                 fprintf(fileID2,'ref_particle  = %u\n' ,PatterN_Extra_dsRemain(iF,1));
%                 fprintf(fileID2,'stiff = %u \n' ,20 );
%                 fprintf(fileID2,'r0 = 1.2 \n' );
%                 fprintf(fileID2,'PBC = 1 \n'  );
%                 fprintf(fileID2,'}\n' );
%             end
%             fclose(fileID2);
            
            MonoDsT =  obj.MonomerDsRemain;
            MultimerDsT = repmat(MonoDsT , N,1) ;
            N_DsT = size(MonoDsT ,1) ;
            AcrossfMonomer_DST= TotalBase*[ repelem((1:N)' ,N_DsT,1)-1 , repelem((1:N)' ,N_DsT,1)-1];
            MultimerDsT_final =MultimerDsT +  AcrossfMonomer_DST ;
            
            file3_name='dSRemainPat.conf'   ;
            fileID3 = fopen(file3_name,'w');           
            for iF=1:size(MultimerDsT_final,1)
                fprintf(fileID3,'{\n' );
                fprintf(fileID3,'type = mutual_trap\n' );
                fprintf(fileID3,'particle = %u\n' ,MultimerDsT_final(iF,1));
                fprintf(fileID3,'ref_particle  = %u\n' ,MultimerDsT_final(iF,2));
                fprintf(fileID3,'stiff = %u \n' ,0.2 );
                fprintf(fileID3,'r0 = 1.2 \n'  );
                fprintf(fileID3,'PBC = 1 \n'  );
                fprintf(fileID3,'}\n' );
                fprintf(fileID3,'{\n' );
                fprintf(fileID3,'type = mutual_trap\n' );
                fprintf(fileID3,'particle = %u\n' ,MultimerDsT_final(iF,2));
                fprintf(fileID3,'ref_particle  = %u\n' ,MultimerDsT_final(iF,1));
                fprintf(fileID3,'stiff = %u \n' ,0.2 );
                fprintf(fileID3,'r0 = 1.2 \n' );
                fprintf(fileID3,'PBC = 1 \n'  );
                fprintf(fileID3,'}\n' );
            end
            fclose(fileID3);            
          
        end
        
        function pSeqAll=get.pSeqAll(obj)
            pSeqAll= obj.pSeq{1}  ;
            for k = 2: length(obj.pSeq)
                pSeqAll= [pSeqAll ,   obj.pSeq{k} ] ;
            end
        end
        
        function Output=get.scafC5All(obj)
            Output= obj.scafC5{1}  ;
            for k = 2: length(obj.scafC5)
                Output= [Output ;   obj.scafC5{k} ] ;
            end
        end
        
        function total_mem = get_mem(obj)
            props = properties(obj);
            total_mem = 0;
            for ii=1:length(props)
                curr_prop = obj.(props{ii});  %#ok<*NASGU>
                s = whos('curr_prop');
                total_mem = total_mem + s.bytes;
            end
        end
        
        function Skipbase=get.skipBase(obj)
            % Skipbase: Dependent property****
            % Depend on
            % obj.containBundle{k}.Default_skipPattern1,obj.containBundle{k}.Default_skipPattern2
            % update for MagiDNA2, July 30,2020
            if obj.CustomSkipAndInsertion ==false % default
                BaseRoute=[];
                %             for k=1: length(obj.scafC5All)
                BaseRouteOneSCaf = interpolateBase( obj.scafC5All ) ;
                BaseRoute=[ BaseRoute ;BaseRouteOneSCaf];
                MaxBase =  size(obj.ScafdigitHC{1} ,1) ;
                InherPosition= [-1,-1 ];
                for k=1 : length(obj.containBundle)
                    if strcmp(obj.containBundle{k}.type ,'SQ' )
                        skipPattern1=obj.containBundle{k}.Default_skipPattern1   ;   % skip mod2 =0
                        skipPattern2=obj.containBundle{k}.Default_skipPattern2   ; 
                        for j = 1 : length(obj.containBundle{k}.Zbase1)
                            skipPattern1=skipPattern1(and(skipPattern1>obj.containBundle{k}.Zbase1(j) , skipPattern1<obj.containBundle{k}.Zbase2(j)) ) ;
                            skipPattern2=skipPattern2(and(skipPattern2>obj.containBundle{k}.Zbase1(j) , skipPattern2<obj.containBundle{k}.Zbase2(j)) ) ;

                            BC = [k ,j ] ;
                            [~,ind] = ismember(BC , obj.RelateTable(:,1:2) , 'rows') ;
                            C4 = obj.RelateTable(ind,4) ;
                            C5 = obj.RelateTable(ind,5) ;
                            if mod(C4,2)==0
                                if ~isempty(skipPattern1)
                                    InherPosition=union(InherPosition, [C5*ones(length(skipPattern1),1) ,skipPattern1'] ,'rows') ;
                                end
                            else
                                if ~isempty(skipPattern2)
                                    InherPosition=union(InherPosition, [C5*ones(length(skipPattern2),1) ,skipPattern2'] ,'rows') ;
                                end
                            end
                        end
                    end
                end
                InherPosition = setdiff(InherPosition , [-1,-1],'rows') ;
                Skipbase=intersect(InherPosition, BaseRoute,'stable','rows') ;
            else
               Skipbase=[-1,-1] ;
               for k =1 :length(obj.containBundle )
                   if ~isempty( obj.containBundle{k}.skipPosition)
                     Skipbase=[Skipbase ;  obj.containBundle{k}.skipPosition(:,3:4) ];
                   end
               end
%                Skipbase=Skipbase(2:end ,:) ;
               Skipbase = setdiff(Skipbase , [-1,-1],'stable','rows') ;
            end            
        end
        
        function InsertBase=get.insertBase(obj)
             if obj.CustomSkipAndInsertion ==false % default
                 InsertBase=[-1,-1];
             else
               InsertBase=[-1,-1] ;
               for k =1 :length(obj.containBundle )
                   if ~isempty( obj.containBundle{k}.InsertPosition)
                     InsertBase=[InsertBase ;  obj.containBundle{k}.InsertPosition(:,3:4) ];
                   end
               end
%                InsertBase=InsertBase(2:end ,:) ;            
             end  
             InsertBase = setdiff(InsertBase , [-1,-1],'stable','rows') ;
        end
                
        function S=saveobj(obj)
            mc = ?hyperbundle ;        % get class infomation
            for k=1:length(mc.PropertyList)
                if  mc.PropertyList(k).Dependent==0   % only independent properties
                    S.(mc.PropertyList(k).Name)=obj.( mc.PropertyList(k).Name);
                end
            end
        end
        function obj=loadobj(obj,S)
            %              ListProperties=properties(obj) ;
            mc = ?hyperbundle ;        % get class infomation
            for k=1:length(mc.PropertyList)
                if  mc.PropertyList(k).Dependent==0   % only independent properties
                    if isfield(S,mc.PropertyList(k).Name)
                        obj.(mc.PropertyList(k).Name)=S.(mc.PropertyList(k).Name);
                        
                    else
                        obj.(mc.PropertyList(k).Name)= [];
                        debuggggg=22 ;
                    end
                end
            end
            
            OldData = 0 ;
            for k=1: length( obj.containBundle )
                if  obj.containBundle{k}.Tol ~=15
                    fprintf('The Tol is %i when saving this data \n' , obj.containBundle{k}.Tol);
                    OldData =1 ;
                end
                
                obj.containBundle{k}.Tol = 15 ;   % debug 02132019
            end
            
            if OldData==1
%                 obj.TakeSeveralV3=cell(size( obj.TakeSeveralV3))  ;  % assign as empty to prevent wrong FB index
                fprintf(' OldData is true............. \n' ) ;
                fprintf(' Please search forced-connection again \n' ) ;
            end
            
            %----Move old Single scaf into cell for multi-scaf, 08/12/2019
            Fields={'ScafRouting','ScafGlobal','scafC4','scafC5','ScafAllBase','Scaf_fromCadDOM'} ;
            if  ~iscell(obj.(Fields{1}))  % old formats for scaffold routing
                for k = 1 : length(Fields)
                    obj.(Fields{k})={ obj.(Fields{k})} ;
                end
            end
            %-----------
            
        end
        
        function Output=get.ForcedConnectListC4(obj)
            RTable=obj.RelateTable;
            TT=obj.ForcedConnectList;
            Output=zeros(size(TT,1),4);
            for i=1:size(TT,1)
                [~,ind]=ismember(TT(i,1:2), RTable(:,1:2),'rows');
                [~,ind2]=ismember(TT(i,4:5), RTable(:,1:2),'rows');
                Output(i,:)=[RTable(ind,4) ,TT(i,3) , RTable(ind2,4) ,TT(i,6) ];
            end
        end
        function Output=get.ForcedConnectListC5(obj)
            RTable=obj.RelateTable;
            TT=obj.ForcedConnectList;
            Output=zeros(size(TT,1),4);
            for i=1:size(TT,1)
                [~,ind]=ismember(TT(i,1:2), RTable(:,1:2),'rows');
                [~,ind2]=ismember(TT(i,4:5), RTable(:,1:2),'rows');
                Output(i,:)=[RTable(ind,5) ,TT(i,3) , RTable(ind2,5) ,TT(i,6) ];
            end
        end
        
        
        
        
        function obj=hyperbundle(type,varargin)
            if type==1
                obj.containBundle={varargin{1}};
            end
            obj=AddBundles(obj,[]);
        end
        
        function obj=AddBundles(obj,Nbun)
            if ~isempty(Nbun)
                obj.containBundle{end+1}=Nbun;
            end
            for ii=1:length( obj.containBundle)
                AA(ii)=  length(obj.containBundle{ii}.Zbase1);
            end
            obj.BundleHasCylinder=AA;
            
            Gindex=1;
            BB=zeros(sum(AA),2);
            for jj=1:length( obj.containBundle)
                BB(Gindex:Gindex+AA(jj)-1 , 1:2)=[jj*ones( AA(jj),1) ,(1:AA(jj))'];
                Gindex=Gindex+ AA(jj);
            end
            obj.GlobalCylinderIndex=BB;
        end
        
        function obj=ConvertScafSQ(obj)  %into digital form
            MCylinderIndex= length(unique(obj.RelateTable(:,4)))  ; % due to extend overhang
            DCellList=cell(MCylinderIndex,1);
            
            SeqI = obj.ScafGlobal{1} ;
            for k= 2: length(obj.ScafGlobal)
                SeqI= [ SeqI ;  obj.ScafGlobal{k} ];
            end
            SeqIAll=SeqI;
            SeqIAll(:,2)=SeqI(:,2);   %+32*ones(size(SeqI(:,2)));  %in case of reach bottom bound
            SeqAll=[SeqIAll(1:end-1,1:2) SeqIAll(2:end,1:2)];
            SeqAll(:,2)= round(SeqAll(:,2));SeqAll(:,4)= round(SeqAll(:,4));
            
            MM=max(max(SeqAll));
            MM=32*(1+ceil(MM/32)) ;   % found bug, Sep 18 2018
            for i=1:MCylinderIndex   %create shell to store information
                DCellList{i}=zeros(MM,4);
            end
            
            for scaf_j=  1:length(obj.ScafGlobal)
                SeqI=obj.ScafGlobal{scaf_j};
                SeqI(:,2)=SeqI(:,2);   %+32*ones(size(SeqI(:,2)));  %in case of reach bottom bound
                Seq=[SeqI(1:end-1,1:2) SeqI(2:end,1:2)];
                Seq(:,2)= round(Seq(:,2));Seq(:,4)= round(Seq(:,4));
                %             MCylinderIndex=size(obj.uniqueInplaneXY,1);
                
                
                %                 MM=max(max(Seq));
                %                 MM=32*(1+ceil(MM/32));   % found bug, Sep 18 2018
                %                 for i=1:MCylinderIndex   %create shell to store information
                %                     DCellList{i}=zeros(MM,4);
                %                 end
                
                C0=-1;  Z0=-1;
                for i1=1:1:size(Seq,1)     % i1=1 3 5 7 ......
                    CylinderIndex=obj.RelateTable(Seq(i1,1),5);
                    CylinderIndex2=obj.RelateTable(Seq(i1,3),5);
                    if Seq(i1,2)>Seq(i1,4) &&  CylinderIndex==CylinderIndex2  %move down
                        for j1=Seq(i1,2):-1:Seq(i1,4)
                            if (i1==size(Seq,1))  && (j1==Seq(i1,4))
                                DCellList{CylinderIndex}(j1,:)=[C0,Z0,-1,-1];
                                break
                            end
                            if    j1==Seq(i1,4)
                            else
                                C1= obj.RelateTable(Seq(i1,1),5 );
                                ZNext=j1-1;
                                DCellList{CylinderIndex}(j1,:)=[C0,Z0,C1,ZNext];
                                C0=CylinderIndex;
                                Z0=j1;
                            end
                        end
                    elseif Seq(i1,2)<Seq(i1,4) &&  CylinderIndex==CylinderIndex2   %move up
                        for j1=Seq(i1,2):1:Seq(i1,4)
                            if (i1==size(Seq,1))  && (j1==Seq(i1,4))
                                DCellList{CylinderIndex}(j1,:)=[C0,Z0,-1,-1];
                                break
                            end
                            if    j1==Seq(i1,4)
                            else
                                C1= obj.RelateTable(Seq(i1,1) ,5 );
                                ZNext=j1+1;
                                DCellList{CylinderIndex}(j1,:)=[C0,Z0,C1,ZNext];
                                C0=CylinderIndex;
                                Z0=j1;
                            end
                        end
                    else    %  Xover
                        C1= CylinderIndex2;    %  C1=  Seq(i1+1,3);
                        ZNext=  Seq(i1,4);
                        Z0=Seq(i1,2) ;
                        DCellList{CylinderIndex}(Z0,:)=[C0,Z0,C1,ZNext]    ;
                        C0=CylinderIndex;
                    end
                end
                for x=1:MCylinderIndex
                    for y=1:length(DCellList{x})
                        if sum(DCellList{x}(y,:))==0
                            DCellList{x}(y,:)=[-1 -1 -1 -1];
                        end
                    end
                end
            end
            obj.ScafdigitSQ=DCellList;
        end  % end of ConvertScafSQ
        
        function obj=ConvertScafHC(obj)  %into digital form
            MCylinderIndex= length(unique(obj.RelateTable(:,4)))  ; % due to extend overhang
            DCellList=cell(MCylinderIndex,1);
            
            SeqI = obj.ScafGlobal{1} ;
            for k= 2: length(obj.ScafGlobal)
                SeqI= [ SeqI ;  obj.ScafGlobal{k} ];
            end
            SeqIAll=SeqI;
            SeqIAll(:,2)=SeqI(:,2);   %+32*ones(size(SeqI(:,2)));  %in case of reach bottom bound
            SeqAll=[SeqIAll(1:end-1,1:2) SeqIAll(2:end,1:2)];
            SeqAll(:,2)= round(SeqAll(:,2));SeqAll(:,4)= round(SeqAll(:,4));
            
            MM=max(max(SeqAll));
            MM=21*(1+ceil(MM/21)) ;   % found bug, Sep 18 2018
            for i=1:MCylinderIndex   %create shell to store information
                DCellList{i}=zeros(MM,4);
            end
            
            
            for scaf_j=  1:length(obj.ScafGlobal)
                SeqI=obj.ScafGlobal{scaf_j};
                SeqI(:,2)=SeqI(:,2);   %+32*ones(size(SeqI(:,2)));  %in case of reach bottom bound
                Seq=[SeqI(1:end-1,1:2) SeqI(2:end,1:2)];
                Seq(:,2)= round(Seq(:,2));Seq(:,4)= round(Seq(:,4));
                
                %                 MM=max(max(Seq));
                %                 MM=21*(1+ceil(MM/21)) ;   % found bug, Sep 18 2018
                %                 for i=1:MCylinderIndex   %create shell to store information
                %                     DCellList{i}=zeros(MM,4);
                %                 end
                
                C0=-1;  Z0=-1;
                for i1=1:1:size(Seq,1)     % i1=1 3 5 7 ......
                    CylinderIndex=obj.RelateTable(Seq(i1,1),5);
                    CylinderIndex2=obj.RelateTable(Seq(i1,3),5);
                    if Seq(i1,2)>Seq(i1,4) &&  CylinderIndex==CylinderIndex2  %move down
                        for j1=Seq(i1,2):-1:Seq(i1,4)
                            if (i1==size(Seq,1))  && (j1==Seq(i1,4))
                                DCellList{CylinderIndex}(j1,:)=[C0,Z0,-1,-1];
                                break
                            end
                            if    j1==Seq(i1,4)
                            else
                                C1= obj.RelateTable(Seq(i1,1),5 );
                                ZNext=j1-1;
                                DCellList{CylinderIndex}(j1,:)=[C0,Z0,C1,ZNext];
                                C0=CylinderIndex;
                                Z0=j1;
                            end
                        end
                    elseif Seq(i1,2)<Seq(i1,4) &&  CylinderIndex==CylinderIndex2   %move up
                        for j1=Seq(i1,2):1:Seq(i1,4)
                            if (i1==size(Seq,1))  && (j1==Seq(i1,4))
                                DCellList{CylinderIndex}(j1,:)=[C0,Z0,-1,-1];
                                break
                            end
                            if    j1==Seq(i1,4)
                            else
                                C1= obj.RelateTable(Seq(i1,1) ,5 );
                                ZNext=j1+1;
                                DCellList{CylinderIndex}(j1,:)=[C0,Z0,C1,ZNext];
                                C0=CylinderIndex;
                                Z0=j1;
                            end
                        end
                    else    %  Xover
                        C1= CylinderIndex2;    %  C1=  Seq(i1+1,3);
                        ZNext=  Seq(i1,4);
                        Z0=Seq(i1,2) ;
                        DCellList{CylinderIndex}(Z0,:)=[C0,Z0,C1,ZNext]    ;
                        C0=CylinderIndex;
                    end
                end
                for x=1:MCylinderIndex
                    for y=1:length(DCellList{x})
                        if sum(DCellList{x}(y,:))==0
                            DCellList{x}(y,:)=[-1 -1 -1 -1];
                        end
                    end
                end
            end
            obj.ScafdigitHC=DCellList;
        end  % end of ConvertScafHC
        
        function obj=ExportJSON(obj)
            dat=loadjson('CadNano.json');
            NNdat=dat;
            %             jj=randi([10,99],1,1);
            %             JJ=strcat('TBund',int2str(jj));
            %             JJ=strcat(JJ,'.json');
            
            prompt = {'Enter File name:'};
            dlg_title = 'Input';
            num_lines = 1;
            defaultans = {'No'};
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            
            if strcmp(answer,'No')
                return;
            end
            for twofile2=1:2
                if twofile2==1
                    file_name=  strcat(answer{1},'SQ');
                    DCellList2=obj.ScafdigitSQ   ;
                    StappCellList2=obj.DigitStapSQ;
                else
                    file_name=  strcat(answer{1},'HC');
                    DCellList2=obj.ScafdigitHC   ;
                    StappCellList2=obj.DigitStapHC;
                    
                end
                
                file_name=strcat(file_name,'.json');
                NNdat.name=file_name;
                ColorChoice=[29184      243362     1507550     3355443     5749504     7536862     8947848    11184640    12060012    13369344    16225054];
                for cylindex=2:size(DCellList2,1)
                    NNdat.vstrands{cylindex} =NNdat.vstrands{1};       %initialized
                end
                %                 DCellList2=ScafCell;
                
                for kj=1:size(DCellList2,1)
                    for k=1:size(DCellList2{kj},1)
                        if DCellList2{kj}(k,1)~=-1
                            DCellList2{kj}(k,1)= obj.RelateVec(DCellList2{kj}(k,1));
                        end
                        if DCellList2{kj}(k,3)~=-1
                            DCellList2{kj}(k,3)= obj.RelateVec(DCellList2{kj}(k,3));
                        end
                    end
                end
                
                for kj2=1:size(DCellList2,1)
                    DCMat=DCellList2{kj2};
                    DeterCylDirectW=vertcat([-1 -1 -1 -1],DCMat);
                    DeterCylDirectW(end,:)=[];
                    DCellList2{kj2}=DeterCylDirectW;
                end
                for kj=1:size(StappCellList2,1)
                    for k=1:size(StappCellList2{kj},1)
                        if StappCellList2{kj}(k,1)~=-1
                            StappCellList2{kj}(k,1)= obj.RelateVec(StappCellList2{kj}(k,1));
                        end
                        if StappCellList2{kj}(k,3)~=-1
                            StappCellList2{kj}(k,3)= obj.RelateVec(StappCellList2{kj}(k,3));
                        end
                    end
                end
                for kj2=1:size(StappCellList2,1)
                    DCMat=StappCellList2{kj2};
                    DeterCylDirectW=vertcat([-1 -1 -1 -1],DCMat);
                    DeterCylDirectW(end,:)=[];
                    StappCellList2{kj2}=DeterCylDirectW;
                end
                MaxBase=size(DCMat,1);
                
                for cylindex=1: size(DCellList2,1)
                    C5index=find( obj.RelateTable(:,5)== cylindex);
                    C5index=C5index(1);
                    Bundle=obj.containBundle{   obj.RelateTable(C5index,1)} ;
                    
                    Cyl2=obj.RelateTable(C5index,2) ;
                    %                     Cyl2=obj.RelateTable(obj.RelateTable(:,5)==C5index,2)
                    
                    %                    if isempty(Cyl2)
                    %                        sdfsf=33;
                    %                    end
                    Cyl2=Cyl2(1);  %C2 express
%                     if strcmp(Bundle.Lattice, 'HoneyComb')  || Cyl2==-1   % not on overhang
%                         skipPattern1=[];
%                         skipPattern2=[];
%                     else
%                         skipPattern1=10:60:MaxBase;   % skip mod2 =0
%                         skipPattern2=40:60:MaxBase;
%                     end
                    
                    
                    Colorstap= find(obj.HeadOfStap(:,1)==cylindex);  %staple index in stapList3
                    ColorM=zeros(length(Colorstap),2);
                    ColorM(:,1)=obj.HeadOfStap(Colorstap,2) + 0*ones(length(Colorstap),1)   ;
                    
                    %                     ColorM(:,2)= ColorChoice(randi([1 length(ColorChoice)],size(Colorstap)))   ;
                    t_json=findobj(gcf,'Tag','t_json') ;
                    for color_j=1 :length(Colorstap)
                        Hex =t_json.Data{Colorstap(color_j),5} ;
                        ColorM(color_j,2)= hex2dec(Hex) ;
                    end
                    
                    
                    
                    
                    ColorM2=sortrows(ColorM,1);
                    ColorM3=[ColorM2;-999,-888];
                    NNdat.vstrands{cylindex}.stap_colors ={ColorM3};
                    %                    [cylindex, obj.RelateVec(cylindex)];
                    NNdat.vstrands{cylindex}.num =obj.RelateVec(cylindex);
                    NNdat.vstrands{cylindex}.scafLoop =cell(0,1);
                    NNdat.vstrands{cylindex}.stap=StappCellList2{cylindex};
                    NNdat.vstrands{cylindex}.scaf=DCellList2{cylindex};     %
                    
                    VerVecSkip=transpose(find(sum(DCellList2{cylindex},2)~=-4));
                    if mod( NNdat.vstrands{cylindex}.num,2)==0
                        
%                         skipPattern1 =  Bundle.Default_skipPattern1 +1 ;
                        skipPattern1 = obj.skipBase(obj.skipBase(:,1)==C5index ,2  )' +1 ; %Matlab and Pyhton indexing
                        SkipIndex=intersect(VerVecSkip,skipPattern1);
                    else
                        skipPattern2 = obj.skipBase(obj.skipBase(:,1)==C5index ,2  )' +1 ; 
                        SkipIndex=intersect(VerVecSkip,skipPattern2);
                    end
                    skipVec=zeros(1,size(DCellList2{1},1));
                    skipVec(SkipIndex)=-1;
                    NNdat.vstrands{cylindex}.skip=round(skipVec);
                    %------add insert July 30 2020
                    VerVecInsert=transpose(find(sum(DCellList2{cylindex},2)~=-4));
 
                        insertPattern = obj.insertBase(obj.insertBase(:,1)==C5index ,2  )' +1 ; 
                        InsertIndex=intersect(VerVecSkip,insertPattern);
                    
                    insertVec=zeros(1,size(DCellList2{1},1));
                    insertVec(InsertIndex)=1;
                    NNdat.vstrands{cylindex}.loop=round(insertVec);                    
                    %----------
                    %                    szofstapVec=size(skipVec)
                    NNdat.vstrands{cylindex}.stapLoop=cell(0,1);
%                     loopVec=zeros(1,size(DCellList2{1},1)); %loopVec(100)=-1;
%                     NNdat.vstrands{cylindex}.loop=loopVec;
                    QQ=find(obj.RelateTable(:,5)==cylindex); QQ=QQ(1);
                    NNdat.vstrands{cylindex}.col=round(obj.RelateTable(QQ,6));
                    NNdat.vstrands{cylindex}.row=round(obj.RelateTable(QQ,7));
                end
                
                %-------adding spacing cylinders
                AddspacingCyl=1;
                if AddspacingCyl==0
                    NNdat2=NNdat;
                elseif AddspacingCyl==1
                    NNdat2=NNdat;
                    nBundle=max(obj.RelateTable(:,1));
                    if ismember(-1 , obj.RelateTable(:,2) ) % with overhange two more spacing
                        extraCyl= 2*(nBundle) ;
                    else
                        extraCyl= 2*(nBundle-1) ;
                    end
                    ExtReTableC5= max(obj.RelateTable(:,5))+1: max(obj.RelateTable(:,5))+extraCyl  ;
                    RestofExt=zeros(3,extraCyl);
                    for k=1:length(ExtReTableC5)
                        RestofExt(2,k)=0;
                        RestofExt(3,k)=k-1 ;
%                         RestofExt(3,k)=0;
%                         RestofExt(2,k)=k-1 ;
                        RestofExt(1,k) =RestofExt(3,k) +max(obj.RelateTable(:,4))+1;
                    end
                    ExtTab=[RestofExt(1,:);ExtReTableC5;RestofExt(2:3,:)]' ;
                    NNdat2.name=NNdat.name;
                    NNdat2.vstrands{1}= NNdat.vstrands{1} ; Cyli=2; extj=1; ToOH_cylinders = 0 ;
                    for k2=2:length(NNdat.vstrands)
                        if ( obj.RelateTable(k2,1)~= obj.RelateTable(k2-1,1) && obj.RelateTable(k2-1,2)>0 && obj.RelateTable(k2,2)>0  )   ||  (obj.RelateTable(k2,2)==-1 && ToOH_cylinders==0)   %cross bundle, extra
                            for tw=1:2
                                %                         NNdat2.vstrands{Cyli}.stap_colors={[-999,-888]};
                                NNdat2.vstrands{Cyli}.stap_colors=cell(0,1);
                                NNdat2.vstrands{Cyli}.num=ExtTab(extj,1);       %ExtTab(extj,2);
                                NNdat2.vstrands{Cyli}.scafLoop=cell(0,1);
                                NNdat2.vstrands{Cyli}.stap = -1*ones(MaxBase,4) ;
                                NNdat2.vstrands{Cyli}.skip= zeros(1,MaxBase);
                                NNdat2.vstrands{Cyli}.scaf = -1*ones(MaxBase,4) ;
                                NNdat2.vstrands{Cyli}.stapLoop=cell(0,1);
                                NNdat2.vstrands{Cyli}.loop =zeros(1,MaxBase);
                                NNdat2.vstrands{Cyli}.col=ExtTab(extj,4);  
                                NNdat2.vstrands{Cyli}.row =ExtTab(extj,3);
%                                 NNdat2.vstrands{Cyli}.col=ExtTab(extj,3);  
%                                 NNdat2.vstrands{Cyli}.row =ExtTab(extj,4);

                                Cyli=Cyli+1;
                                extj=extj+1;
                            end
                            if obj.RelateTable(k2,2)==-1
                                ToOH_cylinders=1 ;
                            end
                        end
                        NNdat2.vstrands{Cyli}= NNdat.vstrands{k2} ;
                        Cyli=Cyli+1;
                    end
                end
                %--------
                TTtext=savejson('Title',NNdat2,'ArrayIndent',0,'Compact',1 );
                
                %                 TTtext=savejson('Title',NNdat,'ArrayIndent',0,'Compact',1 );
                TTtext(1:10)=[];
                TTtext(end-1:end)=[];
                IOfSC2=strfind(TTtext, ',[-999,-888]');  %for color json export
                %                 Cop=IOfSC2;
                %                 for removedd=1:length(IOfSC2)
                %                    Exxtraindex= Cop(1);
                %                      TTtext(Exxtraindex:Exxtraindex+11)=[];
                %                      Cop=strfind(TTtext, ',[-999,-888]');
                %                 end
                
                
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
                
                [status, msg, msgID] = mkdir('cadnano_jsonfile') ;  % Feb 13 2020, create the folder saving the JSON file, preventing  this folder won't be create in package process, Feb 13 2020
                JsonFolder=[ pwd filesep 'cadnano_jsonfile' filesep];
                fileID = fopen([JsonFolder file_name],'w');
                fprintf(fileID,TTtext);
                fclose(fileID);
                fprintf('print json file %s  \n',file_name) ;
                %                 ReadTest=loadjson([JsonFolder file_name])
            end
            
            
        end %end of ExportJSON
        
        function obj=findRT(obj)
            RTable=zeros(size(obj.GlobalCylinderIndex,1), 7)  ;
            RTable(:,1:2)=obj.GlobalCylinderIndex;
            % [ [1 2]->[bundle, localindex], [3]->Gindex, [4]->Candnano ]
            nC=size(RTable,1);
            for ii=1:nC
                RTable(ii,3)=ii;
            end
            
            
            
            UniXY=cell(length(obj.containBundle),1);
            for iB=1:length(UniXY)
                QQ  =unique(obj.containBundle{iB}.CylInplanePosition,'rows');
                UniXY{iB}=[iB*ones(size(QQ,1),1) QQ];
            end
            AllUniXY=[];
            for iB2=1:length(UniXY)
                AllUniXY=[ AllUniXY; UniXY{iB2}  ];
            end
            obj.uniqueInplaneXY=AllUniXY;
            
            for col5=1:nC
                bi=RTable(col5,1);  % bundle index
                bundle=obj.containBundle{RTable(col5,1)} ;
                XY=bundle.CylInplanePosition(RTable(col5,2),:)  ;
                [~,searchR]= ismember( [bi XY], AllUniXY,'rows');
                RTable(col5,5)=  searchR;
            end
            
            CountEven=0;  CountOdd=1;
            for col4=1:nC
                bundle=obj.containBundle{RTable(col4,1)} ;
                CylisAGroup=ismember(RTable(col4,2) ,bundle.AGroup);
                CylGoUp= ~xor(CylisAGroup , bundle.AGroupGoUp) ;
                switch CylGoUp
                    case 1
                        if  sum(RTable(:,5)==RTable(col4,5))==1
                            RTable(col4,4)=  CountEven ;
                            CountEven=CountEven+2;
                        else
                            WW=find(RTable(:,5)==RTable(col4,5));
                            if  col4==min(WW)
                                RTable(WW,4)=  CountEven;
                                CountEven=CountEven+2;
                            end
                        end
                    case 0
                        if  sum(RTable(:,5)==RTable(col4,5))==1
                            RTable(col4,4)=  CountOdd ;
                            CountOdd=CountOdd+2;
                        else
                            WW=find(RTable(:,5)==RTable(col4,5));
                            if  col4==min(WW)
                                RTable(WW,4)=  CountOdd;
                                CountOdd=CountOdd+2;
                            end
                        end
                end
            end
            
            %             for correctC4=1:nC   % correct Column 4
            %                 ind=find(RTable(:,5)==RTable(correctC4,5));
            %                 ind=ind(1);
            %                  RTable(correctC4,4)=RTable(ind,4);
            %             end
            
            [ HClattice ] = findHClattice( 1 ,[40 40]) ;
            %             ConstXDel_HC= HClattice.HCcenter(8,1)-HClattice.HCcenter(7,1);
            ConstXDel_HC=sqrt(3) ;
            
            shiftamongDifBundle=[0 0];
            for Bi=1:length(obj.containBundle)
                updateShift=shiftamongDifBundle;
                bundle=obj.containBundle{Bi} ;
                %                BundleminXY=min(bundle.CylInplanePosition)/(bundle.CylRadius*2);
                %                BundlemaxXY=max(bundle.CylInplanePosition)/(bundle.CylRadius*2);
                InthisBundle= find(RTable(:,1)==Bi);
                RefCyl=RTable(InthisBundle(1),2);
                if strcmp(bundle.Lattice, 'Square')
                    REfXY=bundle.CylInplanePosition(RefCyl,:)/(bundle.CylRadius*2)+ shiftamongDifBundle ;  %--------only for SQ
                    BundleminXY=min(bundle.CylInplanePosition)/(bundle.CylRadius*2);
                else
                    NX=  round(bundle.CylInplanePosition(RefCyl,1)/ConstXDel_HC)   + shiftamongDifBundle(1);   % only X
                    NY=  round( bundle.CylInplanePosition(RefCyl,2)/3)  +  shiftamongDifBundle(2)+1 ;   % only Y
                    REfXY=[NX,NY];
                    
                    mm1=round(min(bundle.CylInplanePosition(:,1))/ConstXDel_HC)+1  ;
                    mm2=round( min(bundle.CylInplanePosition(:,2))/3)+1 ;
                    BundleminXY=[mm1,mm2] ;
                    %                 sdfsf=2344
                end
                
                CylisAGroup=ismember(RefCyl,bundle.AGroup);
                CylGoUp= ~xor(CylisAGroup , bundle.AGroupGoUp);
                
                %                if xor(~xor(CylGoUp==1 ,mod(sum(REfXY),2)==1) , strcmp(bundle.Lattice, 'HoneyComb') )
                %                   comp=1;
                %                else
                %                  comp=0;
                %                end
                
                if xor(~xor(CylGoUp==1 ,mod(sum(REfXY),2)==1) , 0 )
                    comp=1;
                else
                    comp=0;
                end
                %
                if mod(sum(BundleminXY),2)==0
                    tar=[2,2]; CompXY=BundleminXY-tar;
                else
                    tar=[2,3];CompXY=BundleminXY-tar;
                end
                
                
                RTable(InthisBundle(1), 6:7)=REfXY+[comp,0]-CompXY ;
                
                for cylj= 2:length( InthisBundle)
                    thisCyl=RTable(InthisBundle(cylj),2);
                    if strcmp(bundle.Lattice, 'Square')
                        thisXY=bundle.CylInplanePosition(thisCyl,:)/(bundle.CylRadius*2)+updateShift ;
                    else
                        NX2=  round(bundle.CylInplanePosition(thisCyl,1)/ConstXDel_HC)   ;   % only X
                        NY2=  round( bundle.CylInplanePosition(thisCyl,2)/3)+1  ;   % only Y
                        thisXY=[NX2,NY2]+updateShift  ;
                    end
                    RTable(InthisBundle(cylj), 6:7)=thisXY+[comp,0]-CompXY;
                    QQ=thisXY+[comp,0]-CompXY;
                    shiftamongDifBundle=max(shiftamongDifBundle,QQ) ;
                end
            end
            obj.RelateTable=RTable;
            obj.RelateTableOrr = RTable ; %-----------for extend overhangs
            %             obj.RelateVec=RTable(:,5);
            
            
            
            [~,index]=sortrows(RTable,5);
            Vec=RTable(:,4);
            TTRelateVec=Vec(index);
            WQER =union(TTRelateVec,[],'stable');
            obj.RelateVec=WQER';
            
            %                 [~,index]=sortrows(obj.RTable',3);
            %                 Vec=obj.RTable(2,:);
            %                 TTRelateVec=Vec(index);
            %                 obj.RelateVec=union(TTRelateVec,[],'stable');
            
            %             C3Expree=obj.ScafGlobal;
            %             C4Express=zeros(size(C3Expree));
            %             C5Express=zeros(size(C3Expree));
            %             for k=1:size(C5Express,1)
            %                 cylC3=C3Expree(k,1);
            %                 C4Express(k,1:2)= [ RTable(cylC3,4)  C3Expree(k,2)] ;
            %                 C5Express(k,1:2)= [ RTable(cylC3,5)  C3Expree(k,2)] ;
            %             end
            %
            %             obj.scafC4=C4Express;
            %             obj.scafC5=C5Express;
            
            C5AdjM=zeros(max(obj.RelateTable(:,5)) , max(obj.RelateTable(:,5)) );  %ex 40 cylinders
            
            for morecyli= 1: size( obj.RelateTable,1) %ex 1:52
                for morecylj= morecyli+1: size( obj.RelateTable,1)
                    bundlei= obj.RelateTable(morecyli,1);
                    bundlej= obj.RelateTable(morecylj,1);
                    if bundlei==bundlej
                        bundlepart=obj.containBundle{bundlei} ;
                        if  bundlepart.CylAdjMat( obj.RelateTable(morecyli,2),  obj.RelateTable(morecylj,2))==1
                            C5AdjM(obj.RelateTable(morecyli,5), obj.RelateTable(morecylj,5))=1 ;
                            C5AdjM(obj.RelateTable(morecylj,5), obj.RelateTable(morecyli,5))=1   ;
                        end
                    end
                end
            end
            obj.CylinderC5AdjM=   C5AdjM;
            %             obj.ForcedConnectListC4;
        end   % end of fcn findRT
        
        function  get_ScafAllBase(obj) % get ScafAllBase in cells
            ScafForScatterAll = [];
            ScafForScatterIndividual = cell( length( obj.scafC5), 1)  ;  %scaffold 2D surf graphic handles in cell, later inserting spacing ~= C5
            
            for scaf_j = 1 : length( obj.scafC5)
                CellMat= obj.scafC5{scaf_j} ;
                BaseRoute = interpolateBase( CellMat ) ;
                ScafForScatter= setdiff(BaseRoute , obj.skipBase,'rows' ,'stable') ;  %C5 notation
                if scaf_j==1
                    ScafForScatterAll=   ScafForScatter ;
                else
                    ScafForScatterAll=   [ScafForScatterAll; ScafForScatter] ;
                end   
                ScafForScatterIndividual{scaf_j} = ScafForScatter ;
            end
            obj.ScafAllBase=ScafForScatterIndividual ;
        end
        
        function StapToScaf_Corr_NBase_FollowObjFunc=InspectRouting_export(obj,   varargin)
%             varargin : t_json,pStapleH,  plotH ,jsonSlider2
            if ~isempty( obj.StapList3)  % evaluate cropping scaffold with staples
                obj.get_ScafAllBase ;   % get ScafAllBase in cells
                NOriGPairBCB = obj.ScafRouting; NBaseOri=0;
                for k=1: length(NOriGPairBCB)
                    BaseRoutOri = interpolateBase_ThreeCol( NOriGPairBCB{k} );
                    nn= size(BaseRoutOri,1) ;
                    NBaseOri=NBaseOri+nn ;
                end
                
                fprintf('Checking sacffold and staple mapping. \n')
                StapToScaf_Corr= zeros( length(obj.StapList3) ,  length(obj.ScafRouting)  ) ;
                StapToScaf_Corr_NBase= zeros( length(obj.StapList3) ,  length(obj.ScafRouting)  ) ;
                
                %-------collect and stack scaffold for efficiency
                AllScafBaseByBase =  zeros(NBaseOri , 4   ) ; cc=1;  % without considering skip [B C B scaf]
                for scafj= 1:  length(obj.ScafRouting)
                    BaseRout2 = interpolateBase_ThreeCol( obj.ScafRouting{scafj} );
                    AllScafBaseByBase(cc:cc+size(BaseRout2,1)-1 ,:) =[BaseRout2, scafj*ones(size(BaseRout2,1),1)] ; cc=cc+size(BaseRout2,1) ;
                end
                
                for stapi=  1:  length(obj.StapList3)
                    C5rep_Stap = obj.StapList3{stapi} ;
                    %        CornerC4 =ScafCornerC4{k} ;
                    CornerBCB= [zeros(size(C5rep_Stap,1),1) ,C5rep_Stap ];
                    [~,b]=ismember(C5rep_Stap(:,1),obj.RelateTable(:,5) ) ;
                    CornerBCB(:,1:2) = obj.RelateTable(b,1:2) ;
                    % NewScaf =CornerBCB;
                    StapBaseBCB =  interpolateBase_ThreeCol( CornerBCB );
                    [~,BB] = ismember(StapBaseBCB , AllScafBaseByBase(:,1:3) ,'rows'   ) ; BB=setdiff(BB,0);
                    BelongToBundles = unique( AllScafBaseByBase(BB,4) )' ;
                    StapToScaf_Corr(stapi, BelongToBundles) =1;
                    
                    [aa,bb]=hist(AllScafBaseByBase(BB,4),unique(AllScafBaseByBase(BB,4))) ;
                    StapToScaf_Corr_NBase(stapi, bb) =aa ;
                end
                StapToScaf_Corr_NBase=StapToScaf_Corr_NBase(1: length(obj.StapList3), 1:length(obj.ScafRouting)  ) ;
                [u,v] = find(StapToScaf_Corr_NBase~=0) ;
                QQ= StapToScaf_Corr_NBase~=0 ;
                Q=sum(QQ,2) ; %Q(randi(size(Q,1) ,5,1))=3 
                
%                  CircularH = InspectCircularMapping(obj ) ;
%                 fH2 = figure(562) ;clf  ;fH2.WindowButtonDownFcn=''; fH2.WindowButtonMotionFcn='';
%                 fH2.Name='Scaffold and Staple mapping (Without considering skip in SQ bundles' ; set(fH2,'NumberTitle','off');

              
%                 subplot(1,4,1) ;
%                 ax=gca; CCorder =ax.ColorOrder;
%                 for k = 1: length(Q)
%                     v = [0 0; 1 0; 1 1; 0 1]+ [zeros(4,1), (k-1)*ones(4,1)]  ;
%                     f = [1 2 3 4];                   
%                     pH= patch('Faces',f,'Vertices',v,'FaceColor',CCorder( mod( Q(k)-1,size(CCorder,1))+1 ,:) ) ;
%                     pH.EdgeColor=0.1*ones(1,3) ; pH.LineWidth=0.1 ;
%                     CircularH{k}.Color =pH.FaceColor ;
%                 end
%                 axis ij ;
%                 %                 fffH =figure(561) ; clf; %h = image(u,v,StapToScaf_Corr_NBase(u,v)) ;
%                 subplot(1,4,[3 4]) ;
%                 h = image(StapToScaf_Corr_NBase) ;
%                 
%                 
% %                 axx =gca; colormap copper;
%                 %                 h.ButtonDownFcn=@(src,evn)showMapping(src,evn,StapToScaf_Corr_NBase);
% %                 fH2.WindowButtonDownFcn  =@(src,evn)showMapping(src,evn,StapToScaf_Corr_NBase,axx,ax ,Q ,CircularH , varargin);
%                                 fH2.WindowButtonMotionFcn =@(src,evn)showMapping(src,evn,StapToScaf_Corr_NBase,axx,ax ,Q ,CircularH , varargin  );
%                 colorbar;
%                 xlabel('Scaffold index');     ylabel('Staple index');
%                 ax.YLim=axx.YLim; ax.XTick=[]; axes(ax); ax.Color='none';
%                 color = get(fH2,'Color');
%                 set(gca,'XColor',[0 0 0],'YColor',[0 0 0 ],'TickDir','out')
%                 axes(axx);  axx.FontSize =14;
%                 
                UU_Corr = unique(StapToScaf_Corr,'rows') ;
                
%                 %        RowsSingle = find(sum(UU_Corr,2)==1 ) ;
%                 for k=1: size(UU_Corr, 1)
%                     [QQ,~] = ismember(StapToScaf_Corr,UU_Corr(k,:),'rows'   ) ;
%                     fprintf( ' Scaf mapping case [%s] has %i staples.   \n',  num2str(UU_Corr(k,:)),sum(QQ) ) ;
%                 end
%                 fprintf( '\n' ) ;
                
                CheckStapleOverMultipleScaf = sum(UU_Corr, 2) ;
                for k= 1 : max(CheckStapleOverMultipleScaf)
                    casesMultiScaf  = CheckStapleOverMultipleScaf==k ;
                    [QQ2,~] = ismember(StapToScaf_Corr,UU_Corr(casesMultiScaf,:),'rows'   ) ;
%                     fprintf( ' Staple across (%s) scaffolds =  %i staples.   \n',  num2str(k),sum(QQ2) ) ;
                end
                fprintf( '\n' ) ; TolScaf =size(StapToScaf_Corr ,2) ; StapToScaf_Corr_NBase_FollowObjFunc=zeros(1,size(StapToScaf_Corr ,2)) ;
                for scafi = 1 : size(StapToScaf_Corr ,2)
                    Ind_ThisScaf = StapToScaf_Corr(:,scafi)==1 ;
                    Ind_ToOther =  StapToScaf_Corr(:,setdiff([1:TolScaf],  scafi    ))==1 ;
                    Ind_ToOther = sum(Ind_ToOther,2)>0   ;
                    ConnectStaple = and( Ind_ThisScaf ,Ind_ToOther) ;
                    fprintf( 'Scaf %s has %i staples that connect to other scaffolds\n',  num2str(scafi),sum(ConnectStaple) ) ;
                    StapToScaf_Corr_NBase_FollowObjFunc(scafi) =sum(ConnectStaple) ;
                    
                end
%                 StapToScaf_Corr_NBase_FollowObjFunc; % ----only numerical
%                 output for this function, remove unnecessary
%                 visualization
%                 StapToScaf_Corr_NBase; colorbar;
%                 fffH.UserData.StapToScaf_Corr_NBase=StapToScaf_Corr_NBase ;
                
%                 hlink =linkprop([ax axx],{ 'YLim' }); % The axes should stay aligned
%                 ax.UserData.hlink=hlink ;
                
            end
           
        end  % end of InspectRouting_export
        
        function InspectRouting(obj,   varargin)
%             varargin : t_json,pStapleH,  plotH ,jsonSlider2
            if ~isempty( obj.StapList3)  % evaluate cropping scaffold with staples
                obj.get_ScafAllBase ;   % get ScafAllBase in cells
                NOriGPairBCB = obj.ScafRouting; NBaseOri=0;
                for k=1: length(NOriGPairBCB)
                    BaseRoutOri = interpolateBase_ThreeCol( NOriGPairBCB{k} );
                    nn= size(BaseRoutOri,1) ;
                    NBaseOri=NBaseOri+nn ;
                end
                
                fprintf('Checking sacffold and staple mapping. \n')
                StapToScaf_Corr= zeros( length(obj.StapList3) ,  length(obj.ScafRouting)  ) ;
                StapToScaf_Corr_NBase= zeros( length(obj.StapList3) ,  length(obj.ScafRouting)  ) ;
                
                %-------collect and stack scaffold for efficiency
                AllScafBaseByBase =  zeros(NBaseOri , 4   ) ; cc=1;  % without considering skip [B C B scaf]
                for scafj= 1:  length(obj.ScafRouting)
                    BaseRout2 = interpolateBase_ThreeCol( obj.ScafRouting{scafj} );
                    AllScafBaseByBase(cc:cc+size(BaseRout2,1)-1 ,:) =[BaseRout2, scafj*ones(size(BaseRout2,1),1)] ; cc=cc+size(BaseRout2,1) ;
                end
                
                for stapi=  1:  length(obj.StapList3)
                    C5rep_Stap = obj.StapList3{stapi} ;
                    %        CornerC4 =ScafCornerC4{k} ;
                    CornerBCB= [zeros(size(C5rep_Stap,1),1) ,C5rep_Stap ];
                    [~,b]=ismember(C5rep_Stap(:,1),obj.RelateTable(:,5) ) ;
                    CornerBCB(:,1:2) = obj.RelateTable(b,1:2) ;
                    % NewScaf =CornerBCB;
                    StapBaseBCB =  interpolateBase_ThreeCol( CornerBCB );
                    [~,BB] = ismember(StapBaseBCB , AllScafBaseByBase(:,1:3) ,'rows'   ) ; BB=setdiff(BB,0);
                    BelongToBundles = unique( AllScafBaseByBase(BB,4) )' ;
                    StapToScaf_Corr(stapi, BelongToBundles) =1;
                    
                    [aa,bb]=hist(AllScafBaseByBase(BB,4),unique(AllScafBaseByBase(BB,4))) ;
                    StapToScaf_Corr_NBase(stapi, bb) =aa ;
                end
                StapToScaf_Corr_NBase=StapToScaf_Corr_NBase(1: length(obj.StapList3), 1:length(obj.ScafRouting)  ) ;
                [u,v] = find(StapToScaf_Corr_NBase~=0) ;
                QQ= StapToScaf_Corr_NBase~=0 ;
                Q=sum(QQ,2) ; %Q(randi(size(Q,1) ,5,1))=3 
                
                 CircularH = InspectCircularMapping(obj ) ;
                fH2 = figure(562) ;clf  ;fH2.WindowButtonDownFcn=''; fH2.WindowButtonMotionFcn='';
                fH2.Name='Scaffold and Staple mapping (Without considering skip in SQ bundles' ; set(fH2,'NumberTitle','off');

              
                subplot(1,4,1) ;
                ax=gca; CCorder =ax.ColorOrder;
                for k = 1: length(Q)
                    v = [0 0; 1 0; 1 1; 0 1]+ [zeros(4,1), (k-1)*ones(4,1)]  ;
                    f = [1 2 3 4];                   
                    pH= patch('Faces',f,'Vertices',v,'FaceColor',CCorder( mod( Q(k)-1,size(CCorder,1))+1 ,:) ) ;
                    pH.EdgeColor=0.1*ones(1,3) ; pH.LineWidth=0.1 ;
                    CircularH{k}.Color =pH.FaceColor ;
                end
                axis ij ;
                %                 fffH =figure(561) ; clf; %h = image(u,v,StapToScaf_Corr_NBase(u,v)) ;
                subplot(1,4,[3 4]) ;
                h = image(StapToScaf_Corr_NBase) ;
                
                
                axx =gca; colormap copper;
                %                 h.ButtonDownFcn=@(src,evn)showMapping(src,evn,StapToScaf_Corr_NBase);
%                 fH2.WindowButtonDownFcn  =@(src,evn)showMapping(src,evn,StapToScaf_Corr_NBase,axx,ax ,Q ,CircularH , varargin);
                                fH2.WindowButtonMotionFcn =@(src,evn)showMapping(src,evn,StapToScaf_Corr_NBase,axx,ax ,Q ,CircularH , varargin  );
                colorbar;
                xlabel('Scaffold index');     ylabel('Staple index');
                ax.YLim=axx.YLim; ax.XTick=[]; axes(ax); ax.Color='none';
                color = get(fH2,'Color');
                set(gca,'XColor',[0 0 0],'YColor',[0 0 0 ],'TickDir','out')
                axes(axx);  axx.FontSize =14;
                UU_Corr = unique(StapToScaf_Corr,'rows') ;
                %        RowsSingle = find(sum(UU_Corr,2)==1 ) ;
                for k=1: size(UU_Corr, 1)
                    [QQ,~] = ismember(StapToScaf_Corr,UU_Corr(k,:),'rows'   ) ;
                    fprintf( ' Scaf mapping case [%s] has %i staples.   \n',  num2str(UU_Corr(k,:)),sum(QQ) ) ;
                end
                fprintf( '\n' ) ;
                
                CheckStapleOverMultipleScaf = sum(UU_Corr, 2) ;
                for k= 1 : max(CheckStapleOverMultipleScaf)
                    casesMultiScaf  = CheckStapleOverMultipleScaf==k ;
                    [QQ2,~] = ismember(StapToScaf_Corr,UU_Corr(casesMultiScaf,:),'rows'   ) ;
                    fprintf( ' Staple across (%s) scaffolds =  %i staples.   \n',  num2str(k),sum(QQ2) ) ;
                end
                fprintf( '\n' ) ; TolScaf =size(StapToScaf_Corr ,2) ;
                for scafi = 1 : size(StapToScaf_Corr ,2)
                    Ind_ThisScaf = StapToScaf_Corr(:,scafi)==1 ;
                    Ind_ToOther =  StapToScaf_Corr(:,setdiff([1:TolScaf],  scafi    ))==1 ;
                    Ind_ToOther = sum(Ind_ToOther,2)>0   ;
                    ConnectStaple = and( Ind_ThisScaf ,Ind_ToOther) ;
                    fprintf( 'Scaf %s has %i staples that connect to other scaffolds\n',  num2str(scafi),sum(ConnectStaple) ) ;
                end
                StapToScaf_Corr_NBase; colorbar;
                fffH.UserData.StapToScaf_Corr_NBase=StapToScaf_Corr_NBase ;
                fH2.UserData.StapToScaf_Corr=StapToScaf_Corr ;
                hlink =linkprop([ax axx],{ 'YLim' }); % The axes should stay aligned
                ax.UserData.hlink=hlink ;
                
            end
           
        end  % end of InspectRouting
        
        
        function  varargout= InspectCircularMapping(obj )
            % only be trigger by InspectRouting to use to UI
            %             ScafBase =obj.ScafAllBase   ;
            [A,B]=cellfun(@size,obj.ScafAllBase) ;
            StapBase =obj.StapAllBase   ;
            %             skipBase= GetHyperB.skipBase ;
            
            %             sdsf=3
            ScafForScatterAll=zeros(sum(A) ,2 ); c=1;
            for scaf_j=1: length(obj.scafC5)
                CellMat= obj.scafC5{scaf_j} ;
                BaseRoute = interpolateBase( CellMat ) ;
                ScafForScatter= setdiff(BaseRoute , obj.skipBase,'rows' ,'stable') ;  %C5 notation
                %                 if scaf_j==1
                %                     ScafForScatterAll=   ScafForScatter ;
                %                 else
                %                     ScafForScatterAll=   [ScafForScatterAll; ScafForScatter] ;
                %                 end
                ScafForScatterAll(c:c+size(ScafForScatter,1)-1 ,:) =ScafForScatter ; c=c+size(ScafForScatter,1) ;
            end
            
            
            [A2,B2]=cellfun(@size,StapBase) ;
            MergeStapForMapping= zeros(sum(A2) ,2 ) ; c=1 ;
            for k=1: length(StapBase)
                MergeStapForMapping(c:c+A2(k)-1 ,:) = StapBase{k} ; c=c+A2(k) ;
            end
            
            [Q,  stapToScaf ] =ismember(MergeStapForMapping ,  ScafForScatterAll ,'rows') ;
            cumA =cumsum(A)  ;
            FindkCir =(stapToScaf-cumA')>0 ;
            LocationCases = sortrows(unique( FindkCir  ,'rows')  ) ;
            [~,BB ]= ismember(FindkCir,LocationCases ,'rows' ) ;
            
            cumA_v2 = [0;cumA(1:end-1)] ;
            FindPhase = stapToScaf-cumA_v2' ;
            FindPhase(FindPhase<0) = 100*max(max(FindPhase)) ;
            if length(A) ==1  % single scaffold case
            BaseIndOnIndv=  FindPhase ;  
            else
            BaseIndOnIndv = min(FindPhase')' ;
            end
            Ref = A(BB) ;
            Phase = (BaseIndOnIndv./Ref)*360 ;
            
            kCir_phase = zeros(length(stapToScaf) ,2 ) ;
            kCir_phase(:,1) =BB ;
            kCir_phase(:,2) =Phase ;
            
            f565=figure(565);clf ;  hold on ;
            f565.Name='Scaffold and Staple circular mapping' ; set(f565,'NumberTitle','off');
            
            r1=0.25 ; theta =[0:360]-90 ;
            for k = 1 : size(A)    % gradient scaffold
                x=r1.*cosd(theta) + k ;
                y=r1.*sind(theta) + 0.5*mod(k,2)  ;
                z= zeros(size(x));col = (1:length(x))*1000;
                surfH=surface([x;x],[y;y],[z;z],[col;col], 'facecol','no', 'edgecol','interp', 'linew',2 ,  'EdgeColor','flat');
                text(k,0.5*mod(k,2)  ,num2str(strcat('Scaf-', ' ',num2str(k)  )) ) ;
            end
            axis equal ;axis off; axis ij;
            A2 ;c= 1 ;
            r2 =0.2 ;  plotH= cell(length(StapBase) ,1) ;
            for stpj = 1 : length(StapBase)
                Sub_kCir_phase = kCir_phase(c:c+A2(stpj)-1, :) ; c=c+A2(stpj) ;
                
                x= r2.*cosd(Sub_kCir_phase(:,2)-90) + Sub_kCir_phase(:,1) ;
                y= r2.*sind(Sub_kCir_phase(:,2)-90) + 0.5*mod(Sub_kCir_phase(:,1),2);
                
                plotH{stpj}=plot(x,y) ;
            end
            ax= gca;
            ax.Position = [0.02 0.02 0.96 0.96] ;
            varargout{1}= plotH ;
            %               dsdfs=3
        end   % end of InspectCircularMapping
        
        function [BundlesIndex, TF_CrossBundles]=ShowStapleInBundles(obj)
            Stap3 = obj.StapList3 ; % C5 Rep
            
            BundlesIndex=zeros(size(Stap3)) ;
            TF_CrossBundles=zeros(size(Stap3)) ;
            for k =1: length(Stap3)
                C5 = Stap3{k}(:,1) ;
                [~,b] =ismember(C5 , obj.RelateTable(:,5) )    ;
                BInd =  obj.RelateTable(b,1) ;
                if  length(unique(BInd)) ==1
                    BundlesIndex(k) =unique(BInd) ;
                    
                else
                    QQ= unique(BInd) ;
                    BundlesIndex(k)= QQ(1) ; % temporiry put the first one.
                    %                                         BundlesIndex(k)= nan ; % temporiry put nan. later can use cell to store.
                    
                    TF_CrossBundles(k)=1 ;
                end
            end
            
            %             sdsf=3
        end
        
        function obj=ConvertScafG(obj)
            obj.ScafGlobal=cell( length(obj.ScafRouting),1);
            obj.scafC4=cell( length(obj.ScafRouting),1);
            obj.scafC5=cell( length(obj.ScafRouting),1);
            for scaf_j = 1 : length(obj.ScafRouting)
                TT =obj.ScafRouting{scaf_j} ;
                
                ScafG =zeros(size(TT,1),2) ;
                for i=1:size(TT,1)
                    [~,ind]=ismember(TT(i,1:2),  obj.GlobalCylinderIndex,'rows');
                    ScafG(i,1:2)=[ind,TT(i,3)];
                end
                %                 obj.ScafGlobal=ScafG;
                RTable=obj.RelateTable;
                C3Expree=ScafG;
                C4Express=zeros(size(C3Expree));
                C5Express=zeros(size(C3Expree));
                for k=1:size(C5Express,1)
                    cylC3=C3Expree(k,1);
                    C4Express(k,1:2)= [ RTable(cylC3,4)  C3Expree(k,2)] ;
                    C5Express(k,1:2)= [ RTable(cylC3,5)  C3Expree(k,2)] ;
                end
                
                
                obj.ScafGlobal{scaf_j}=ScafG;
                obj.scafC4{scaf_j}=C4Express;
                obj.scafC5{scaf_j}=C5Express;
            end
            
        end
        
        function findCycleList_maxScafXover(obj,~,mainHandles,type)
            %              AllScafXover;
            
            figure(5);clf;
            subplot(2,3,1);axCell=cell(6,1);
            for k=1:6
                subplot(2,3,k);
                axCell{k}=gca;
            end
            
            
            AssignScafXovers=[];
            %             PrePairInfo=mainHandles.PremPair;
            if type==2
                PrePairInfo=mainHandles;
            else
                PrePairInfo=mainHandles.PremPair;
            end
            %             obj.OverRoute=obj.ssOption.OverRoute;
            nPrePaired=0;
            for nPPCount=1:length(obj.containBundle)
                nPrePaired=nPrePaired+ size(PrePairInfo.CellPairList{nPPCount},1);
            end
            
            OriGPairBCB=cell(nPrePaired,1);  %No over route
            count=1;
            OBJOR=obj.OverRoute;
            for bunI=1:length(PrePairInfo.CellPairList)   % Use prePair to genearte initial cylces
                for cycJ=1:size( PrePairInfo.CellPairList{bunI} ,1)
                    TwoCyl=PrePairInfo.CellPairList{bunI}(cycJ,:);   %cylinder index
                    TT=[bunI, TwoCyl(1), obj.containBundle{bunI}.Zbase1(TwoCyl(1))- OBJOR ;...
                        bunI, TwoCyl(1), obj.containBundle{bunI}.Zbase2(TwoCyl(1))+ OBJOR ;...
                        bunI, TwoCyl(2), obj.containBundle{bunI}.Zbase2(TwoCyl(2))+ OBJOR ;...
                        bunI, TwoCyl(2), obj.containBundle{bunI}.Zbase1(TwoCyl(2))- OBJOR ];
                    %check cycle direction--
                    FirstEdgeGoUp = ~xor( ismember(TT(1,2), obj.containBundle{TT(1,1)}.AGroup),  obj.containBundle{TT(1,1)}.AGroupGoUp==1);
                    ReallyGoUpInTT = TT(2,3)>TT(1,3);
                    if xor(FirstEdgeGoUp,ReallyGoUpInTT)
                        TT=flip(TT);
                    end
                    %------
                    OriGPairBCB{ count}=TT;
                    count=count+1 ;
                end
            end
            %             figure(5566);clf;hold on; ax2=gca ; axis off;
            %                         obj.drawBCBdata(OriGPairBCB);
            obj.drawBCBdata(OriGPairBCB,axCell{1});
            
            ConnectList=obj.ForcedConnectList;
            % random shuffle to increase variety
            randOrder = randperm(size(ConnectList,1)/2 ) ;
            randOrder2=reshape([2*randOrder'-1, 2*randOrder']' ,size(ConnectList,1),1) ;
            ConnectList=ConnectList( randOrder2 ,:) ;
            
            
            %             ClistEasy=obj.CListsortCListNotFollow53(ConnectList);
            NOriGPairBCB=OriGPairBCB;
            %-----------------
            %--------------
            for AddXoveri=1:size(ConnectList,1)/2
                Xover=ConnectList(2*AddXoveri-1:2*AddXoveri,:);
                [NOriGPairBCB,UnUsedX, CycleInvolved] =obj.AddXover22cycles(NOriGPairBCB,Xover);
                if ~isempty(UnUsedX) && length(CycleInvolved)==1
                    DD=NOriGPairBCB;
                    NOriGPairBCB=obj.AddXover2OneCycle(NOriGPairBCB,Xover);
                end
                
                %                 obj.drawBCBdata(NOriGPairBCB,ax2);
            end
            %                          obj.drawBCBdata(NOriGPairBCB)
            obj.drawBCBdata(NOriGPairBCB,axCell{2});
            %                          return%--------------------------
            
            % tic
            %             GetHyperB.getScafXoverAsStap  ;
            
            AllScafXover = obj.AllScafXover ;
            %             ConnectList=obj.ForcedConnectList;
            % random shuffle to increase variety
            %             randOrder = randperm(size(AllScafXover,1)/2 ) ;
            %             randOrder2=reshape([2*randOrder'-1, 2*randOrder']' ,size(AllScafXover,1),1) ;
            %            AllScafXover=AllScafXover( randOrder2 ,:) ;
            
            ForceCons =[ConnectList(:,1:3)  ;ConnectList(:,4:6)  ];
            ttol =6;
            for AddXoveri=1:size(AllScafXover,1)/2
                %                   if  AddXoveri==355
                %                   AddXoveri
                %                   end
                Xover=AllScafXover(2*AddXoveri-1:2*AddXoveri,:);
                %---------
                %                 ConnectList ;
                
                [af,Inds3p] =ismember( ForceCons(:,1:2),Xover(1,1:2),'rows');
                %                 if  af~=0
                if min(abs(ForceCons(af,3)-Xover(1,3) ) )< ttol   % ignore some of scaf Xovers
                    Xover
                    
                    continue ;
                end
                %                 end
                [af,Inds3p] = ismember( ForceCons(:,1:2),Xover(1,4:5),'rows') ;
                %                 if  af~=0
                if min(abs(ForceCons(af,3)-Xover(1,6) ) )< ttol   % ignore some of scaf Xovers
                    Xover
                    continue ;
                end
                %                 end
                %-----------
                [NOriGPairBCB,UnUsedX, CycleInvolved] =obj.AddXover22cycles(NOriGPairBCB,Xover);
                if ~isempty(UnUsedX) && length(CycleInvolved)==1
                    DD=NOriGPairBCB;
                    NOriGPairBCB=obj.AddXover2OneCycle(NOriGPairBCB,Xover);
                end
                %                                 obj.drawBCBdata(NOriGPairBCB)     ;
            end
            % toc
            %             obj.drawBCBdata(NOriGPairBCB)
            fprintf('Above Xovers are ignored due to closure to forced connections.  \n')
            obj.drawBCBdata(NOriGPairBCB,axCell{3});
            %             SaveAxes ;
            
            Adj2=zeros(length(NOriGPairBCB) ,length(NOriGPairBCB) );
            for icycle=1:length(NOriGPairBCB)   % version 2 to find AdjM of graph
                for jcycle=icycle+1:length(NOriGPairBCB)
                    %                      for kscOver= 1:2:size(AllScafXover,1)
                    Cyclei = NOriGPairBCB{icycle} ;
                    Cyclej = NOriGPairBCB{jcycle} ;
                    %                     [ia,ib]=ismember( Cyclei(:,1:2) , Cyclej(:,1:2) ,'rows')  ;
                    [~,ib]=ismember( Cyclei, [AllScafXover(:,1:3);AllScafXover(:,4:6)]  ,'rows') ; ib(ib>size(AllScafXover,1)) =ib(ib>size(AllScafXover,1))-size(AllScafXover,1) ;
                    [~,ib2]=ismember( Cyclej, [AllScafXover(:,1:3);AllScafXover(:,4:6)]  ,'rows') ;ib2(ib2>size(AllScafXover,1)) =ib2(ib2>size(AllScafXover,1))-size(AllScafXover,1) ;
                    Sib=sort(round(ib/2) ) ;
                    Sib2=sort(round(ib2/2) ) ;
                    uniAll = unique([Sib;Sib2])  ;uniAll=setdiff(uniAll,0) ;
                    %                     [N,m] = histcounts(Sib,uniAll)      ;
                    
                    %                     [N2,m2] = histcounts(Sib2,uniAll)     ;
                    %                 if sum(Sib==0)==length(Sib)
                    %                     N=Sib
                    if ~isempty(uniAll)
                        N=countCM(Sib,uniAll) ;
                        N2=countCM(Sib2,uniAll) ;
                        XoverBetween= and(N2==2,  N==2) ;
                        if sum(XoverBetween)> 0
                            Adj2( icycle,jcycle)=1;
                            Adj2( jcycle,icycle)=1;
                        end
                    end
                end
            end
            
            
            
            
            gAdjM2=graph(Adj2);
            %              obj.drawBCBdata(NOriGPairBCB,axCell{2});
            
            varietyoption=2 ;
            switch varietyoption
                case 1
                    [T,predM] = minspantree(gAdjM2,'Method', 'dense');  % old way, variety is low
                case 2
                    [T,predM] = minspantree(gAdjM2,'Method', 'dense', 'Root' , randi(length(NOriGPairBCB) ));
            end
            %               figure(45);clf;
            axes(axCell{4});
            p=plot(gAdjM2) ;     highlight(p,T,'EdgeColor','r','LineWidth',1.5)
            %              figure(5566);clf; p=plot(gAdjM2) ; highlight(p,T,'EdgeColor','r','LineWidth',1.5);
            
            fprintf('spantree =  %s \n', num2str(predM) )
            if  sum(isnan(predM))>0
                % some loops can't not integrate after adding forced connection. May be due to too short for bundles or limited scaf xover.
                figure(452);clf; p=plot(gAdjM2) ;  title(' cycle graph, some of them are isolated  ')
                fprintf('some loops can''t not integrate after adding forced connection \n')
                %                 error('Error. check forced connection, lengths of bundle, or limited scaf xover ')
            end
            
            RemoveXovers =zeros(2*(length(predM)-1-sum(isnan(predM)) ) ,6) ;ck=1;
            for k=1: length(predM)
                if predM(k)~=0 && ~isnan( predM(k))
                    [k,predM(k)]  ;
                    Cyclei = NOriGPairBCB{k} ;
                    Cyclej = NOriGPairBCB{predM(k)} ;
                    [~,ib]=ismember( Cyclei, [obj.AllScafXover(:,1:3);obj.AllScafXover(:,4:6)]  ,'rows') ;
                    ib(ib>size(obj.AllScafXover,1)) =ib(ib>size(obj.AllScafXover,1))-size(obj.AllScafXover,1) ;
                    [~,ib2]=ismember( Cyclej, [obj.AllScafXover(:,1:3);obj.AllScafXover(:,4:6)]  ,'rows') ;
                    ib2(ib2>size(obj.AllScafXover,1)) =ib2(ib2>size(obj.AllScafXover,1))-size(obj.AllScafXover,1) ;
                    
                    Sib=sort(round(ib/2) ) ;
                    Sib2=sort(round(ib2/2) ) ;
                    uniAll = unique([Sib;Sib2])  ;uniAll=setdiff(uniAll,0) ;
                    %                     [N,~] = hist(Sib,uniAll)      ;
                    %                     [N2,~] = hist(Sib2,uniAll)     ;
                    N=countCM(Sib,uniAll) ;
                    N2=countCM(Sib2,uniAll) ;
                    
                    XoverBetween= and(N2==2,  N==2) ;
                    XoverInd =uniAll(XoverBetween) ;XoverInd=XoverInd(1) ;
                    
                    XoverBridgeTwoCycle = obj.AllScafXover(2*XoverInd-1:2*XoverInd,:) ;
                    RemoveXovers(ck:ck+1,:) =XoverBridgeTwoCycle ;ck=ck+2;
                    %                     NOriGPairBCB= AddBridgeXover2TwoCycle(obj,NOriGPairBCB,XoverBridgeTwoCycle) ;
                    
                end
            end
            
            %      RemoveXovers
            for k=1:2:size(RemoveXovers,1)
                [k , size(NOriGPairBCB) ] ;
                Xover= RemoveXovers(k:k+1,:) ;
                
                NOriGPairBCB= AddBridgeXover2TwoCycle(obj,NOriGPairBCB,Xover,RemoveXovers) ;
                %                   obj.drawBCBdata(NOriGPairBCB)     ;
                if k==1
                    obj.drawBCBdata(NOriGPairBCB,axCell{5});
                end
            end
            obj.drawBCBdata(NOriGPairBCB,axCell{6});
            
            
            
            if length(NOriGPairBCB)~=1  % use old algorithm, using vertical Xovers
                fprintf('Due to Forced connection, the code found some isoloated scaf circuits without using vertical Xover. Use the old algorithm to integrate scaf cycles  \n')
                Adj2old=zeros(length(NOriGPairBCB) ,length(NOriGPairBCB) );
                for icycle=1:length(NOriGPairBCB)   % version 2 to find AdjM of graph
                    for jcycle=icycle+1:length(NOriGPairBCB)
                        XX=Given2CylceFindXoverList(obj,NOriGPairBCB{icycle},NOriGPairBCB{jcycle},AssignScafXovers, ConnectList  ) ;
                        if ~isempty(XX) %&& ~isempty(YY)
                            Adj2old( icycle,jcycle)=1;
                            Adj2old( jcycle,icycle)=1;
                        end
                    end
                end
                gAdjM2old=graph(Adj2old);
                
                varietyoption=1 ;
                switch varietyoption
                    case 1
                        [~,predMOld] = minspantree(gAdjM2old,'Method', 'sparse');  % old way, variety is low
                    case 2
                        [~,predMOld] = minspantree(gAdjM2old,'Method', 'sparse', 'Root' , randi(length(NOriGPairBCB) ));
                end
                fprintf('spantree(Old) =  %s \n', num2str(predMOld) )  ;
                
                if  sum(isnan(predMOld))>0
                    % some loops can't not integrate after adding forced connection. May be due to too short for bundles or limited scaf xover.
                    figure(4456); p=plot(predMOld) ;  title(' cycle graph, some of them are isolated even using vertical scaf Xovers   ')
                    fprintf('some loops can''t not integrate after adding forced connection \n')
                    error('Error. check forced connection, lengths of bundle, or limited scaf xover ')
                end
                CombiningXover=zeros( 2*(length(predMOld)-1),6);
                for k=2: length(predMOld)
                    [k,predMOld(k)];
                    OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{k},NOriGPairBCB{predMOld(k)},[] ) ;
                    if isempty(OneXover)
                        OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{predMOld(k)},NOriGPairBCB{k}) ;
                    end
                    CombiningXover(2*k-3:2*k-2,:)= OneXover;
                end
                for kk=1:size(CombiningXover,1)/2
                    AppXover=CombiningXover(2*kk-1:2*kk,:);
                    [NOriGPairBCB,UnUsedX, CycleInvolved] = obj.AddXover22cycles(NOriGPairBCB,AppXover);
                end
                
            end
            
            
            
            
            obj.drawBCBdata(NOriGPairBCB)     ;
            NewscafR2=  NOriGPairBCB{1} ;
            %----------
            scafR=NOriGPairBCB{1};
            
            Scaf0=scafR;
            % fixing overrouting at ends------
            ResortCL=[ConnectList(:,1:3) ;ConnectList(:,4:6)] ;
            NofEnds=size(obj.RelateTable,1);
            EndsOnBottom=zeros(NofEnds,3);   nE1=1;
            EndsOnTop=zeros(NofEnds,3);   nE2=1;
            for Bi=1:length( obj.containBundle)
                NCylinBi=length(obj.containBundle{Bi}.Zbase1);
                AA=[Bi*ones(NCylinBi,1),(1:NCylinBi)',obj.containBundle{Bi}.Zbase1'-obj.OverRoute];  %detect
                EndsOnBottom(nE1:NCylinBi+nE1-1,:)=AA;
                nE1=nE1+NCylinBi;
                
                BB=[Bi*ones(NCylinBi,1),(1:NCylinBi)',obj.containBundle{Bi}.Zbase2'+obj.OverRoute];
                EndsOnTop(nE2:NCylinBi+nE2-1,:)=BB;
                nE2=nE2+NCylinBi;
            end
            
            [isBottom,~]=  ismember(ResortCL,EndsOnBottom,'rows');
            TargetNodeofBottom=ResortCL(isBottom,:);
            [~,indComp1]= ismember(TargetNodeofBottom,scafR,'rows');
            
            [isTop2,~]= ismember(scafR,EndsOnTop,'rows');
            [isBottom2,~]= ismember(scafR,EndsOnBottom,'rows');
            
            [~,Bx] = ismember(scafR(indComp1,1:3),obj.ssOption.individual(:,1:3) ,'rows')    ;
            if isempty(obj.ssOption.individual)
                CompV1=0;
            else
                CompV1= obj.ssOption.individual(Bx,4);
            end
            %-- - - - - -
            [isTop,~]=  ismember(ResortCL,EndsOnTop,'rows')  ;
            TargetNodeofTop=ResortCL(isTop,:);
            [~,indComp2]= ismember(TargetNodeofTop,scafR,'rows');
            
            [~,Bx] = ismember(scafR(indComp2,1:3),obj.ssOption.individual(:,1:3) ,'rows');
            %                CompV2= obj.ssOption.individual(Bx,4);
            if isempty(obj.ssOption.individual)
                CompV2=0;
            else
                CompV2= obj.ssOption.individual(Bx,4);
            end
            %---------------------
            isTop2A=setdiff(find(isTop2),indComp2);  isBottom2A=setdiff(find(isBottom2),indComp1);
            scafR(isTop2A,3)= scafR(isTop2A,3)-obj.OverRoute+obj.ssOption.OverRoute-obj.ssOption.BundleShiftTwoSide( scafR(isTop2A,1), 2)   ;   %compensate normal overrouting scaf
            scafR(isBottom2A,3)= scafR(isBottom2A,3)+obj.OverRoute-obj.ssOption.OverRoute+obj.ssOption.BundleShiftTwoSide( scafR(isBottom2A,1),1)   ;   %compensate
            % bug: Sep 20, flip two sides' shifting
            scafR(indComp2,3)= scafR(indComp2,3)+CompV2   ;   %compensate  isTop forceconnection
            scafR(indComp1,3)= scafR(indComp1,3)-CompV1 ;   %compensate  isBottom belong forceconnection
            %--------------
            scafR2= scafR;
            
            searchDebugBunCyl = [4,2];
            for updateSS=1:length(obj.ssOption.ForceConnUpdate)
                OldsPM2= obj.ssOption.ForceConnUpdate{updateSS}.BCBPM2 ;
                NewEnds=obj.ssOption.ForceConnUpdate{updateSS}.BCBCur ;
                [~,iiid1]= ismember(Scaf0,OldsPM2(1:3),'rows') ; fd1=find(iiid1);
                scafR2( fd1,:)=NewEnds(1:3)+scafR(fd1,:)-Scaf0(fd1,:)   ;
                [~,iiid2]= ismember(Scaf0,OldsPM2(4:6),'rows') ; fd2=find(iiid2);
                scafR2(fd2,:)=NewEnds(4:6)+scafR(fd2,:)-Scaf0(fd2,:)    ;
            end
            scafR3= scafR2;
            scafR= scafR2;
            %---------
            if strcmp(obj.ScafOption.prevStack,'scafloop')
                %                     TargetPolyBundles= [9,10 , 11];
                TargetPolyBundles= obj.ScafOption.PolyBundle ;
                for Vedge = 2:2 :size(scafR2 ,1)-1
                    %                     [Vedge,size(scafR2,1)]
                    MM =   scafR2(Vedge:Vedge+1 ,:) ;
                    if ismember(MM(1,1) , TargetPolyBundles)
                        Bundle= obj.containBundle{MM(1,1)} ;
                        if MM(1,1) ==MM(2,1) && MM(1,3) ==MM(2,3) && (MM(1,3)<Bundle.Zbase1(MM(1,2)) || MM(1,3)>Bundle.Zbase2(MM(1,2)) ) % locate out of z1z2
                            
                            CylMaster=MM(1,2) ;
                            CylInSlave=MM(2,2) ;
                            CinMoveUp=  ~xor( ismember(CylMaster,Bundle.AGroup),Bundle.AGroupGoUp);
                            rnd= [MM(1,3)-10,MM(1,3),MM(1,3)+10]   ;XOinZ=zeros(1,3) ;
                            if strcmp(Bundle.Lattice, 'Square')            %change for hybrid structure
                                XOinZ(1)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(1),CinMoveUp);
                                XOinZ(2)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(2),CinMoveUp);
                                XOinZ(3)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(3),CinMoveUp);
                            else
                                %                                 XOinZ= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd,CinMoveUp);
                                XOinZ(1)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(1),CinMoveUp);
                                XOinZ(2)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(2),CinMoveUp);
                                XOinZ(3)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(3),CinMoveUp);
                            end
                            
                            QQ = abs(XOinZ-MM(1,3));
                            choose= XOinZ(QQ==min(QQ)) ;
                            if choose<Bundle.Zbase1(MM(1,2))
                                choose=choose+1;
                            end
                            scafR2(Vedge:Vedge+1 ,3) =choose ;
                        end
                    end
                end
                
                %-----------------
                % hard code
                %                 BundleNoScafLoop = 10:17 ;    fprintf('BundleNoScafLoop = %s \n', num2str(BundleNoScafLoop)  ) ;
                BundleNoScafLoop=obj.ScafOption.BundleNoScafLoop ;
                for k =1:size(scafR2,1)
                    if ismember(scafR2(k,1) ,BundleNoScafLoop)
                        BB=scafR2(k,1) ;
                        if   scafR2(k,3) < obj.containBundle{BB}.Zbase1(scafR2(k,2))   % over-route on Z1
                            scafR2(k,3)=obj.containBundle{BB}.Zbase1(scafR2(k,2))  ;
                        elseif scafR2(k,3) > obj.containBundle{BB}.Zbase2(scafR2(k,2))   % over-route on Z1
                            scafR2(k,3)=obj.containBundle{BB}.Zbase2(scafR2(k,2))   ;
                        end
                    end
                end
                %---------
            elseif  strcmp(obj.ScafOption.prevStack,'polyT')
                BundleNoScafLoop=obj.ScafOption.BundleNoScafLoop  ;
                polyTBundle= setdiff(1:length(obj.containBundle) , BundleNoScafLoop) ;
                for Vedge = 2:2 :size(scafR2 ,1)-1
                    MM =   scafR2(Vedge:Vedge+1 ,:) ;
                    if ismember(MM(1,1) , polyTBundle)
                        Bundle= obj.containBundle{MM(1,1)} ;
                        if MM(1,1) ==MM(2,1) && MM(1,3) ==MM(2,3) && (MM(1,3)<Bundle.Zbase1(MM(1,2)) || MM(1,3)>Bundle.Zbase2(MM(1,2)) ) % locate out of z1z2, exclude interal Xovers
                            
                            CylMaster=MM(1,2) ;
                            CylInSlave=MM(2,2) ;
                            CinMoveUp=  ~xor( ismember(CylMaster,Bundle.AGroup),Bundle.AGroupGoUp);
                            if MM(1,3)>Bundle.Zbase2(MM(1,2)) % Z2 side
                                rnd= [MM(1,3)-20,MM(1,3)-10,MM(1,3)]   ;
                            else
                                rnd= [MM(1,3),MM(1,3)+10,MM(1,3)+20]   ;
                            end
                            
                            XOinZ=zeros(1,3) ;
                            if strcmp(Bundle.Lattice, 'Square')            %change for hybrid structure
                                XOinZ(1)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(1),CinMoveUp);
                                XOinZ(2)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(2),CinMoveUp);
                                XOinZ(3)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(3),CinMoveUp);
                            else
                                %                                 XOinZ= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd,CinMoveUp);
                                XOinZ(1)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(1),CinMoveUp);
                                XOinZ(2)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(2),CinMoveUp);
                                XOinZ(3)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(3),CinMoveUp);
                            end
                            
                            if MM(1,3)>Bundle.Zbase2(MM(1,2)) % Z2 side
                                QQ = XOinZ-Bundle.Zbase2(MM(1,2));
                                QQ(QQ>0)=200;
                                QQ=abs(QQ) ;
                            else   % Z1 side
                                QQ = XOinZ-Bundle.Zbase1(MM(1,2));
                                QQ(QQ<0)=200;
                                QQ=abs(QQ) ;
                            end
                            
                            %                                 QQ = abs(XOinZ-MM(1,3));
                            
                            choose= XOinZ(QQ==min(QQ)) ;
                            if MM(1,3)<Bundle.Zbase1(MM(1,2))
                                choose=choose+1;
                            end
                            scafR2(Vedge:Vedge+1 ,3) =choose ;
                        end
                    end
                end
                BundleNoScafLoop=obj.ScafOption.BundleNoScafLoop ;
                for k =1:size(scafR2,1)
                    if ismember(scafR2(k,1) ,BundleNoScafLoop)
                        BB=scafR2(k,1) ;
                        if   scafR2(k,3) < obj.containBundle{BB}.Zbase1(scafR2(k,2))   % over-route on Z1
                            scafR2(k,3)=obj.containBundle{BB}.Zbase1(scafR2(k,2))  ;
                        elseif scafR2(k,3) > obj.containBundle{BB}.Zbase2(scafR2(k,2))   % over-route on Z1
                            scafR2(k,3)=obj.containBundle{BB}.Zbase2(scafR2(k,2))   ;
                        end
                    end
                end
            end
            
            %--------
            obj.ScafRouting=scafR2   ;
            
            
            
            
        end  % end of function findCycleList_maxScafXover
        
        function ScafXovers=getXoverinScaf(obj,scafR)
            
            AllXovers= zeros(2000,4 ) ;cc=1 ;
            
            Corners = sortrows(scafR{1}) ;
            Cylinders = unique( Corners(:,1:2) ,'rows' ) ;
            for cyli = 1 :size(Cylinders,1)
                [aa,bb] = ismember( Corners(:,1:2),Cylinders(cyli,1:2)  ,'rows') ;
                Base = Corners(aa,3) ;
                
                dBase = diff(Base) ==1 ;
                Xovers = [Base(dBase),Base(dBase)+1] ;
                n =size(Xovers,1) ;
                AllXovers(cc:cc+n-1 ,:) =[ ones(n,1)*Cylinders(cyli,1:2),Xovers];
                cc=cc+n;
                
            end
            AllXovers=AllXovers(1:cc-1,:) ;
            
            [~,Ind1] =  ismember( AllXovers(:,1:3)  ,scafR{1},'rows') ;
            [~,Ind2] =  ismember( AllXovers(:,[1 2 4])  ,scafR{1},'rows') ;
            
            AllXovers= [AllXovers,Ind1,Ind2 ] ;
            %             AllXovers=[ AllXovers , Ind1+Ind2]
            
            sumV = Ind1+Ind2 ;
            %             notCase = ismember( )
            
            ScafXovers = zeros( ceil(size(AllXovers,1)/2) ,12 ) ; count =1 ;
            for srchi = 1: size(AllXovers,1)
                if mod(AllXovers(srchi ,5) ,2)==0
                    P1 = AllXovers(srchi , [1 2 3])  ;
                    P2 = AllXovers(srchi , [1 2 4])  ;
                    
                    ind = find( AllXovers(:,5)== AllXovers(srchi ,5)+1) ;
                    if ~isempty(ind) %&& abs(OriIndP1-OriIndP2)==1
                        P3 = AllXovers(ind , [1 2 4])  ;
                        P4 = AllXovers(ind , [1 2 3])  ;
                        
                        [~,OriIndP1 ]= ismember(P1, scafR{1},'rows') ;
                        [~,OriIndP2 ]= ismember(P2, scafR{1},'rows') ;
                        [~,OriIndP3 ]= ismember(P3, scafR{1},'rows') ;
                        [~,OriIndP4 ]= ismember(P4, scafR{1},'rows') ;
                        %                    [OriIndP1,OriIndP2,OriIndP3,OriIndP4]
                        
                        ScafXovers(count,: )=[P1,P2,P3,P4] ; count=count+1;
                    end
                end
            end
            ScafXovers=ScafXovers(1:count-1, :) ;
            
            
            
            %             sdf=4
        end
        
        
        function NewNOriGPairBCB=AddBridgeXover2TwoCycle(obj,NOriGPairBCB,XoverBridgeTwoCycle,RemoveXovers)
            NewNOriGPairBCB=NOriGPairBCB ;
            InvolveCycle = [];
            FourArr =[XoverBridgeTwoCycle(1,1:3) ; XoverBridgeTwoCycle(1,4:6 ) ;XoverBridgeTwoCycle(2,1:3)    ; XoverBridgeTwoCycle(2,4:6 ) ];
            
            shiftii=0;  shiftjj=0;
            cw=1;
            while 1   % while loop to prevend
                sbreak=0;
                for icycle=1:length(NOriGPairBCB)   % version 2 to find AdjM of graph
                    for jcycle=icycle+1:length(NOriGPairBCB)
                        %                      for kscOver= 1:2:size(AllScafXover,1)
                        Cyclei = circshift(NOriGPairBCB{icycle} , shiftii) ;
                        Cyclej = circshift(NOriGPairBCB{jcycle} ,shiftjj);
                        
                        XoverBridgeTwoCycle;
                        [~, cycci] = ismember(FourArr,Cyclei  ,'rows') ;
                        [~, cyccj] = ismember(FourArr,Cyclej  ,'rows') ;
                        Indii= find(cycci) ;
                        Indjj= find(cyccj) ;
                        
                        if isempty(setdiff(1:4,unique([Indii;Indjj]) ) ) && length(union(1:4,unique([Indii;Indjj]) ) )==4
                            InvolveCycle=[icycle,jcycle] ;sbreak=1 ;
                            break    ;
                        end
                        
                    end
                    if sbreak==1
                        break;
                    end
                end
                
                %                 Cyclei = NOriGPairBCB{InvolveCycle(1)} ;
                %                 Cyclej = NOriGPairBCB{InvolveCycle(2)} ;
                
                %             if min(Indii)==1
                %             NewNOriGPairBCB{InvolveCycle(1)} =[ Cyclei(Indii(2)-1 , :) ; Cyclej(Indjj(2)+1:end ,:) ; Cyclej(1:Indjj(1)-1 ,:) ; Cyclei(Indii(1)+1, :)] ;
                
                %             elseif min(Indjj)~=1
                %                  min(Indii)~=1
                %                  min(Indjj)~=1
                %                  if isempty(Indii)
                %                      sf=3
                %                  end
                
                if  ~ismember(1,cycci)   &&  ~ismember(1,cyccj)
                    break;
                else
                    if rand()>0.5
                        shiftii=shiftii+2;
                    else
                        shiftjj=shiftjj+2;
                        
                    end
                end
                cw=cw+1 ;
            end
            Indii=setdiff(cycci,0) ;   Indjj=setdiff(cyccj,0) ;
            
            NewNOriGPairBCB{InvolveCycle(1)} =[ Cyclei(1:Indii(1)-1 , :) ; Cyclej(Indjj(2)+1:end ,:) ; Cyclej(1:Indjj(1)-1 ,:) ; Cyclei(Indii(2)+1:end, :)] ;
            NewNOriGPairBCB(InvolveCycle(2))=[];
        end
        
        function Output=IntegrateScaffold(obj ,NOriGPairBCB)
            %  mainly modify from findCycleList, Use for multi-scaffold maximizing the numbers of tangle staples.
            NOriGPairBCBInitial =NOriGPairBCB ;
            Adj2=zeros(length(NOriGPairBCB) ,length(NOriGPairBCB) );
            for icycle=1:length(NOriGPairBCB)   % version 2 to find AdjM of graph
                for jcycle=icycle+1:length(NOriGPairBCB)
                    
                    %                     XX=Given2CylceFindXoverList(obj,NOriGPairBCB{icycle},NOriGPairBCB{jcycle}) ;
                    %                     YY=Given2CylceFindXoverList(obj,NOriGPairBCB{jcycle},NOriGPairBCB{icycle}) ;
                    XX=Given2CylceFindXoverList(obj,NOriGPairBCB{icycle},NOriGPairBCB{jcycle},[]  ) ;
                    
                    if ~isempty(XX) %&& ~isempty(YY)
                        Adj2( icycle,jcycle)=1;
                        Adj2( jcycle,icycle)=1;
                    end
                end
            end
            gAdjM2=graph(Adj2);
            [T2,predM2] = minspantree(gAdjM2,'Type','forest' )  ;
            %            figure(235); hh= plot(gAdjM2) ;highlight(hh,T2,'EdgeColor','r','LineWidth',1.5) ;
            cc=1 ;
            CombiningXover=zeros( 2*(length(predM2)-sum(predM2==0)),6);
            for k=1: length(predM2)
                if predM2(k) ~=0
                    %                         [k,predM(k)];
                    OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{k},NOriGPairBCB{predM2(k)},[] ) ;
                    if isempty(OneXover)
                        OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{predM2(k)},NOriGPairBCB{k}) ;
                    end
                    CombiningXover(cc:cc+1,:)= OneXover;
                    cc=cc+2 ;
                end
            end
            
            for kk=1:size(CombiningXover,1)/2
                AppXover=CombiningXover(2*kk-1:2*kk,:);
                [NOriGPairBCB,UnUsedX, CycleInvolved] = obj.AddXover22cycles(NOriGPairBCB,AppXover);
            end
            SaveNOriGPairBCB=  NOriGPairBCB ;
            
            tol =6 ;
            [ ~,whichCyl ] = CheckScafRXover(NOriGPairBCB, tol  ,[]) ;  % or add more bundle
            [~,iy1]=ismember( whichCyl, CombiningXover(:,1:3),'rows' ) ;
            [~,iy2]=ismember( whichCyl, CombiningXover(:,4:6),'rows' ) ;
            iyy= union(iy1(iy1~=0) ,iy2(iy2~=0) ) ;
            iyyEven= zeros( length(iyy) ,2) ;
            for ci= 1:length(iyy)
                if mod(iyy(ci),2)==1
                    iyyEven(ci,:) = [iyy(ci) ,iyy(ci)+1];
                else
                    iyyEven(ci,:) = [iyy(ci)-1 ,iyy(ci)];
                end
            end
            iyyEven= unique(iyyEven,'rows') ;
            
            RemoveXovers =CombiningXover(iyyEven' ,:) ;
            
            
            %                RemoveXovers =CombiningXover(iyyEven' ,:) ;
            
            nNotGood =size(whichCyl ,1)  ;   % if initial solution is too bad, just back to outer loop to get a better one
            if  nNotGood< 500
                %                     fprintf( ' Seems to be able to find a good scaf \n')
                for tworound= 1:2
                    if tworound==1
                        SaveAppliedXoverInFirstRound = zeros(200, 6) ; kx =1;
                    else
                        fprintf('start second round \n') ;
                        if nNotGood==0 || isempty(RemoveXovers) ; break;end         % in case that first round is good
                        
                        RemoveXovers= SaveAppliedXoverInFirstRound(1:kx-1 , :) ; % second round try to remove newly added Xovers
                        IndRemain=   find( or( ismember( RemoveXovers(:,1:3), whichCyl2,'rows'),ismember( RemoveXovers(:,4:6), whichCyl2,'rows')) ) ;  % only reserve xover has problem
                        IndRemain= union(IndRemain(mod(IndRemain,2)==1), IndRemain(mod(IndRemain,2)==0)-1 )  ; IndRemain=union(IndRemain, IndRemain+1) ;
                        RemoveXovers=RemoveXovers( IndRemain,:) ;
                        %                             RemoveXovers
                    end
                    
                    for riXover=  1:2:size(RemoveXovers,1)
                        SS=NOriGPairBCB ;
                        Xover = RemoveXovers(riXover:riXover+1 ,:) ;
                        TwoCycles=    removeScafXover_general(obj,NOriGPairBCB,Xover)   ;
                        %                               removeScafXover_general(GetHyperB,NOriGPairBCB,Xover)
                        [~, XoverList]=obj.Given2CylceFindXoverList(TwoCycles{1},TwoCycles{2},  []) ; % avoid in while loop. Query one time for a list, 09/11/2019
                        
                        if isempty(XoverList) || length(TwoCycles)>2
                            EEemptyHappened =1 ;
                            continue;
                        end
                        cw= 1;
                        while 1
                            OneXover = XoverList(  randi(size(XoverList,1),1) , :)  ;
                            OneXover= [OneXover(1,1:6) ;OneXover(1,7:12) ];
                            [NOriGPairBCBOut,UnUsedX, CycleInvolved] = obj.AddXover22cycles(TwoCycles,OneXover)     ;
                            [ ~,whichCyl2 ] = CheckScafRXover(NOriGPairBCBOut{1}, tol,[]) ;
                            cw=cw+1;
                            if cw>30
                                whichCyl2;
%                                 fprintf('reach max loop\n');
                                break;   %
                            end
                            
                            if size(whichCyl2,1) < nNotGood
                                nNotGood= size(whichCyl2,1) ;
                                NOriGPairBCB= NOriGPairBCBOut;
                                if length(NOriGPairBCB) ==2
                                    sdsf=3
                                end
                                
                                SaveAppliedXoverInFirstRound(kx:kx+1,: )= OneXover ; kx=kx+2 ;
                                break
                            end
                            if nNotGood==0
                                break
                            end
                        end
%                         fprintf('nCornerOnBadXover = %i\n',   size(whichCyl2,1))
                    end
                end  % for tworound= 1:2
            end
            [ Goodcc,whichCylcc ] = CheckScafRXover(NOriGPairBCB{1}, tol,[] ) ; % whichCylcc
            
            Output =NOriGPairBCBOut ;
            if length(Output)==2
                sdfsf=3;
            end
            
            %            obj.ScafRouting=NOriGPairBCB ;
        end
        
        
        function [obj,Goodcc]=findCycleList(obj,~,mainHandles,type)
            
            %             profile on
            
            obj.GetOrgBlockAdj ;
            %             figure(6);clf;
            %             subplot(2,3,1);axCell=cell(6,1);
            %             for k=1:6
            %                 subplot(2,3,k);
            %                 axCell{k}=gca;
            %             end
            
            tol=obj.ScafOption.minDist_btwXover;
            AssignScafXovers=[];
            %             PrePairInfo=mainHandles.PremPair;
            if type==2
                PrePairInfo=mainHandles;
            else
                PrePairInfo=mainHandles.PremPair;
            end
            %             sdff=23424
            
            %             obj.OverRoute=obj.ssOption.OverRoute;
            nPrePaired=0;
            for nPPCount=1:length(obj.containBundle)
                nPrePaired=nPrePaired+ size(PrePairInfo.CellPairList{nPPCount},1);
            end
            
            
            OriGPairBCB=cell(nPrePaired,1);  %No over route
            count=1;
            OBJOR=obj.OverRoute;
            for bunI=1:length(PrePairInfo.CellPairList)   % Use prePair to genearte initial cylces
                for cycJ=1:size( PrePairInfo.CellPairList{bunI} ,1)
                    TwoCyl=PrePairInfo.CellPairList{bunI}(cycJ,:);   %cylinder index
                    TT=[bunI, TwoCyl(1), obj.containBundle{bunI}.Zbase1(TwoCyl(1))- OBJOR ;...
                        bunI, TwoCyl(1), obj.containBundle{bunI}.Zbase2(TwoCyl(1))+ OBJOR ;...
                        bunI, TwoCyl(2), obj.containBundle{bunI}.Zbase2(TwoCyl(2))+ OBJOR ;...
                        bunI, TwoCyl(2), obj.containBundle{bunI}.Zbase1(TwoCyl(2))- OBJOR ];
                    %check cycle direction--
                    FirstEdgeGoUp = ~xor( ismember(TT(1,2), obj.containBundle{TT(1,1)}.AGroup),  obj.containBundle{TT(1,1)}.AGroupGoUp==1);
                    ReallyGoUpInTT = TT(2,3)>TT(1,3);
                    if xor(FirstEdgeGoUp,ReallyGoUpInTT)
                        TT=flip(TT);
                    end
                    %------
                    OriGPairBCB{ count}=TT;
                    count=count+1 ;
                end
            end
            
            %             figure(5566);clf;hold on; ax2=gca ; axis off;
            %              obj.drawBCBdata(OriGPairBCB,axCell{1});
            
            %             obj.drawBCBdata(OriGPairBCB,ax2);
            ConnectList=obj.ForcedConnectList;
            if ~isempty(obj.SaveInternalXovers)  %-------Muli-Scaf, isolate regions
                AddInterXovers = obj.SaveInternalXovers ;
                AddInterXovers2 = zeros(size(AddInterXovers,1)*2  ,size(AddInterXovers,2)/2) ;
                for k = 1:size(AddInterXovers,1)
                    AddInterXovers2(2*k-1,:) =AddInterXovers(k, 1:6) ;
                    AddInterXovers2(2*k,:) =AddInterXovers(k, 7:12) ;
                end
                AddInterXovers2= [AddInterXovers2(:,4:6) ,AddInterXovers2(:,1:3)] ;
                
                ConnectList=[ConnectList;  AddInterXovers2] ;
            end
            
            %----------------------
            ss_json=findall(gcf,'Tag','ss_json') ;
            sH=findobj(ss_json,'Tag','sH') ;
            if ~isempty(sH)
                if isfield(sH.UserData,'ExtraForcedScafXover')  % && size(obj.ScafRouting ,1)==size(sH.Parent.UserData.PlottedScafR ,1 )
                    if size(obj.ScafRouting{1} ,1)==size(sH.Parent.UserData.PlottedScafR{1} ,1)
                        if sum(sum(obj.ScafRouting{1}==sH.Parent.UserData.PlottedScafR{1}) )== numel(sH.Parent.UserData.PlottedScafR{1})
                            % check the routing showing in cadnano tab is the
                            % current one.
                            nExtra = size(sH.UserData.ExtraForcedScafXover,1)/2 ;
                            str=strcat('ExtraForcedXover (',num2str(nExtra),')' ) ;
                            answerExtra = questdlg('Found Extra Overs in cadnano tab. Do you want to use it ?',str, ...
                                'Yes','No' ,'Yes');
                            if strcmp(answerExtra,'Yes')
                                %                         ConnectList= [ConnectList ; sH.UserData.ExtraForcedScafXover] ;
                                %                         ConnectList
                                ExtraXover =sH.UserData.ExtraForcedScafXover;
                                OriginalScafRouting  = obj.ScafRouting ;
                                
                                for AddXoveri=1:size(ExtraXover,1)/2  % add extra xover from cadnano tab
                                    Xover=ExtraXover(2*AddXoveri-1:2*AddXoveri,:);
                                    [OriginalScafRouting,UnUsedX, CycleInvolved] =obj.AddXover22cycles(OriginalScafRouting,Xover);
                                    if ~isempty(UnUsedX) && length(CycleInvolved)==1
                                        DD=OriginalScafRouting;
                                        OriginalScafRouting=obj.AddXover2OneCycle(OriginalScafRouting,Xover);
                                    end
                                end
                                Adj2=zeros(length(OriginalScafRouting) ,length(OriginalScafRouting) );
                                for icycle=1:length(OriginalScafRouting)   % version 2 to find AdjM of graph
                                    for jcycle=icycle+1:length(OriginalScafRouting)
                                        XX=Given2CylceFindXoverList(obj,OriginalScafRouting{icycle},OriginalScafRouting{jcycle},AssignScafXovers, ExtraXover  ) ;
                                        if ~isempty(XX) %&& ~isempty(YY)
                                            Adj2( icycle,jcycle)=1;
                                            Adj2( jcycle,icycle)=1;
                                        end
                                    end
                                end
                                gAdjM2=graph(Adj2);
                                [T,predM] = minspantree(gAdjM2,'Method', 'sparse');  % old way, variety is low
                                if  sum(isnan(predM))>0
                                    % some loops can't not integrate after adding forced connection. May be due to too short for bundles or limited scaf xover.
                                    figure; p=plot(gAdjM2) ;  title(' cycle graph, some of them are isolated  ')
                                    fprintf('some loops can''t not integrate after adding forced connection \n')
                                    error('Error. check forced connection, lengths of bundle, or limited scaf xover ')
                                end
                                CombiningXover=zeros( 2*(length(predM)-1),6);
                                for k=2: length(predM)
                                    [k,predM(k)];
                                    OneXover=obj.Given2CylceFindXoverList(OriginalScafRouting{k},OriginalScafRouting{predM(k)},[] ) ;
                                    if isempty(OneXover)
                                        OneXover=obj.Given2CylceFindXoverList(OriginalScafRouting{predM(k)},OriginalScafRouting{k}) ;
                                    end
                                    CombiningXover(2*k-3:2*k-2,:)= OneXover;
                                end
                                for kk=1:size(CombiningXover,1)/2
                                    AppXover=CombiningXover(2*kk-1:2*kk,:);
                                    [OriginalScafRouting,UnUsedX, CycleInvolved] = obj.AddXover22cycles(OriginalScafRouting,AppXover);
                                    
                                end
                                obj.ScafRouting= OriginalScafRouting ;
                                [ Goodcc,whichCylcc ] = CheckScafRXover(OriginalScafRouting{1}, tol,obj.ScafOption.BundleNoScafLoop ) ;
                                %                         sdsdf=3
                                
                                return ;
                            end
                        end
                    end
                end
            end
            %--------------
            %              fw=msgbox('Using hard code for pairing !! ');
            %             waitfor(fw);
            %              warning('Using hard code in Scaf algorithm')
            %             for k=4:6   % hard code
            %                ConnectList(end+1,:) = [k,5,74,k,3,74];
            %                ConnectList(end+1,:) = [k,3,75,k,5,75];
            %                ConnectList(end+1,:) = [k,4,78,k,6,78];
            %                ConnectList(end+1,:) = [k,6,79,k,4,79];
            %             end
            %---------------
            
            
            % random shuffle to increase variety
            randOrder = randperm(size(ConnectList,1)/2 ) ;
            randOrder2=reshape([2*randOrder'-1, 2*randOrder']' ,size(ConnectList,1),1) ;
            ConnectList=ConnectList( randOrder2 ,:) ;
            
            
            %             ClistEasy=obj.CListsortCListNotFollow53(ConnectList);
            NOriGPairBCB=OriGPairBCB;
            %-----------------
            %--------------
            for AddXoveri=1:size(ConnectList,1)/2
                Xover=ConnectList(2*AddXoveri-1:2*AddXoveri,:);
                [NOriGPairBCB,UnUsedX, CycleInvolved] =obj.AddXover22cycles(NOriGPairBCB,Xover);
                if ~isempty(UnUsedX) && length(CycleInvolved)==1
                    DD=NOriGPairBCB;
                    NOriGPairBCB=obj.AddXover2OneCycle(NOriGPairBCB,Xover);
                end
                
                %                 obj.drawBCBdata(NOriGPairBCB,ax2);
            end
            %               obj.drawBCBdata(NOriGPairBCB,axCell{2});
            %             obj.drawBCBdata(NOriGPairBCB,ax2)
            %                          obj.drawBCBdata(NOriGPairBCB)
            %                          return%--------------------------
            
            Adj2=zeros(length(NOriGPairBCB) ,length(NOriGPairBCB) );
            for icycle=1:length(NOriGPairBCB)   % version 2 to find AdjM of graph
                for jcycle=icycle+1:length(NOriGPairBCB)
                    
                    %                     XX=Given2CylceFindXoverList(obj,NOriGPairBCB{icycle},NOriGPairBCB{jcycle}) ;
                    %                     YY=Given2CylceFindXoverList(obj,NOriGPairBCB{jcycle},NOriGPairBCB{icycle}) ;
                    XX=Given2CylceFindXoverList(obj,NOriGPairBCB{icycle},NOriGPairBCB{jcycle},AssignScafXovers, ConnectList  ) ;
                    
                    if ~isempty(XX) %&& ~isempty(YY)
                        Adj2( icycle,jcycle)=1;
                        Adj2( jcycle,icycle)=1;
                    end
                end
            end
            gAdjM2=graph(Adj2);
            
            varietyoption=1 ;
            %             switch varietyoption
            %                 case 1
            %                     [T,predM] = minspantree(gAdjM2,'Method', 'sparse');  % old way, variety is low
            %                 case 2
            %                     [T,predM] = minspantree(gAdjM2,'Method', 'sparse', 'Root' , randi(length(NOriGPairBCB) ));
            %             end
            
            %-------------- multiple scaffolds, 08/08/2019
            %             if  sum(isnan(predM))>0
            [T2,predM2] = minspantree(gAdjM2,'Type','forest' ) ; predM2
            
            figure(235); hh= plot(gAdjM2) ;highlight(hh,T2,'EdgeColor','r','LineWidth',1.5) ;
            
            cc=1 ;
            CombiningXover=zeros( 2*(length(predM2)-sum(predM2==0)),6);
            for k=1: length(predM2)
                
                if predM2(k) ~=0
                    %                         [k,predM(k)];
                    OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{k},NOriGPairBCB{predM2(k)},[] ) ;
                    if isempty(OneXover)
                        OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{predM2(k)},NOriGPairBCB{k}) ;
                    end
                    CombiningXover(cc:cc+1,:)= OneXover;
                    cc=cc+2 ;
                end
            end
            
            for kk=1:size(CombiningXover,1)/2
                AppXover=CombiningXover(2*kk-1:2*kk,:);
                [NOriGPairBCB,UnUsedX, CycleInvolved] = obj.AddXover22cycles(NOriGPairBCB,AppXover);
            end
            obj.ScafRouting=NOriGPairBCB ;
            
%             figure(234);clf; obj.drawBCBdata(NOriGPairBCB,gca) ;
            %                 dfgdg=4
            %                 figure; obj.drawBCBdata(NOriGPairBCB,gca) ;
            
            %                 sdfsf=3
            %                 figure(556) ;clf ;hold on ;
            %                 obj.plotScafR_cylindermodelMulti ;
            %                 return
            %             end
            
            
            
            %             axes(axCell{4});
            %               p=plot(gAdjM2) ;     highlight(p,T,'EdgeColor','r','LineWidth',1.5)
            %              predM
            %              obj.drawBCBdata(OriGPairBCB,axCell{1});
            %                        figure; p=plot(gAdjM2) ;    % highlight(p,T,'EdgeColor','r','LineWidth',1.5)
            
            %             fprintf('spantree =  %s \n', num2str(predM) ) ;
            
            %             if  sum(isnan(predM))>0
            %                 % some loops can't not integrate after adding forced connection. May be due to too short for bundles or limited scaf xover.
            %                 figure; p=plot(gAdjM2) ;  title(' cycle graph, some of them are isolated  ')
            %                 fprintf('some loops can''t not integrate after adding forced connection \n')
            %                 error('Error. check forced connection, lengths of bundle, or limited scaf xover ')
            %             end
            
            
            %             switch varietyoption
            %                 case 1
            %                     CombiningXover=zeros( 2*(length(predM)-1),6);
            %                     for k=2: length(predM)
            %                         [k,predM(k)];
            %                         OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{k},NOriGPairBCB{predM(k)},[] ) ;
            %                         if isempty(OneXover)
            %                             OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{predM(k)},NOriGPairBCB{k}) ;
            %                         end
            %                         CombiningXover(2*k-3:2*k-2,:)= OneXover;
            %                     end
            %                 case 2
            %                     AssignScafXovers = [] ;  % for breakable bundle feature, hasn't integrated
            %                     CombiningXover=zeros( 2*(length(predM)-1),6); cK=2;   %Bund2Xover=1;
            %                     for k=1: length(predM)
            %                         if predM(k)~=0
            %                             CurrAllScafxover=[ AssignScafXovers;CombiningXover];
            %                             nw=1;
            %                             while 1
            %                                 OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{k},NOriGPairBCB{predM(k)} , CurrAllScafxover,ConnectList ) ;
            %                                 if isempty(OneXover)
            %                                     OneXover=obj.Given2CylceFindXoverList(NOriGPairBCB{predM(k)},NOriGPairBCB{k},[]   ) ;
            %                                 end
            %                                 Range=[131,134];
            %                                 if OneXover(1,1)~=2 || nw>2
            %                                     break
            %                                 end
            %                                 nw=nw+1;
            %                             end
            %                             CombiningXover(2*cK-3:2*cK-2,:)= OneXover;
            %                             cK=cK+1;
            %                         end
            %                     end
            %             end
            
            
            
            
            %             for kk=1:size(CombiningXover,1)/2
            %                 AppXover=CombiningXover(2*kk-1:2*kk,:);
            %                 [NOriGPairBCB,UnUsedX, CycleInvolved] = obj.AddXover22cycles(NOriGPairBCB,AppXover);
            %                 %              obj.drawBCBdata(NOriGPairBCB);
            %             end
            
            %           obj.drawBCBdata(NOriGPairBCB,ax2);
            %             tol=10;
            
            
            Goodcc=0;
            %               obj.drawBCBdata(NOriGPairBCB,axCell{5});
            
            %             if length(NOriGPairBCB)==1
            for scaf_j = 1 : length(NOriGPairBCB)
                scafR=NOriGPairBCB{scaf_j};
                
                Scaf0=scafR;
                % fixing overrouting at ends------
                ResortCL=[ConnectList(:,1:3) ;ConnectList(:,4:6)] ;
                NofEnds=size(obj.RelateTable,1);
                EndsOnBottom=zeros(NofEnds,3);   nE1=1;
                EndsOnTop=zeros(NofEnds,3);   nE2=1;
                for Bi=1:length( obj.containBundle)
                    NCylinBi=length(obj.containBundle{Bi}.Zbase1);
                    AA=[Bi*ones(NCylinBi,1),(1:NCylinBi)',obj.containBundle{Bi}.Zbase1'-obj.OverRoute];  %detect
                    EndsOnBottom(nE1:NCylinBi+nE1-1,:)=AA;
                    nE1=nE1+NCylinBi;
                    
                    BB=[Bi*ones(NCylinBi,1),(1:NCylinBi)',obj.containBundle{Bi}.Zbase2'+obj.OverRoute];
                    EndsOnTop(nE2:NCylinBi+nE2-1,:)=BB;
                    nE2=nE2+NCylinBi;
                end
                
                [isBottom,~]=  ismember(ResortCL,EndsOnBottom,'rows');
                TargetNodeofBottom=ResortCL(isBottom,:);
                [TFindComp1,indComp1]= ismember(TargetNodeofBottom,scafR,'rows');
                indComp1=indComp1(TFindComp1) ;
                
                [isTop2,~]= ismember(scafR,EndsOnTop,'rows');
                [isBottom2,~]= ismember(scafR,EndsOnBottom,'rows');
                
                [~,Bx] = ismember(scafR(indComp1,1:3),obj.ssOption.individual(:,1:3) ,'rows')    ;
                if isempty(obj.ssOption.individual)
                    CompV1=0;
                else
                    CompV1= obj.ssOption.individual(Bx,4);
                end
                %-- - - - - -
                [isTop,~]=  ismember(ResortCL,EndsOnTop,'rows')  ;
                TargetNodeofTop=ResortCL(isTop,:);
                [TFindComp2,indComp2]= ismember(TargetNodeofTop,scafR,'rows');
                indComp2=indComp2(TFindComp2) ;
                [~,Bx] = ismember(scafR(indComp2,1:3),obj.ssOption.individual(:,1:3) ,'rows');
                %                CompV2= obj.ssOption.individual(Bx,4);
                if isempty(obj.ssOption.individual)
                    CompV2=0;
                else
                    Ind0 = find(Bx==0) ;  % temp to give some value, later set to zeros
                    Bx(Ind0) =1 ;
                    CompV2= obj.ssOption.individual(Bx,4);
                    CompV2(Ind0) =0 ;
                end
                %---------------------
                isTop2A=setdiff(find(isTop2),indComp2);  isBottom2A=setdiff(find(isBottom2),indComp1);
                %               scafR(indComp2,3)= scafR(indComp2,3)-obj.OverRoute+CompV2   ;   %compensate  isTop forceconnection
                %               scafR(indComp1,3)= scafR(indComp1,3)+obj.OverRoute-CompV1 ;   %compensate  isBottom belong forceconnection
                
                scafR(isTop2A,3)= scafR(isTop2A,3)-obj.OverRoute+obj.ssOption.OverRoute-obj.ssOption.BundleShiftTwoSide( scafR(isTop2A,1), 2)   ;   %compensate normal overrouting scaf
                scafR(isBottom2A,3)= scafR(isBottom2A,3)+obj.OverRoute-obj.ssOption.OverRoute+obj.ssOption.BundleShiftTwoSide( scafR(isBottom2A,1),1)   ;   %compensate
                % bug: Sep 20, flip two sides' shifting
                
                scafR(indComp2,3)= scafR(indComp2,3)+CompV2   ;   %compensate  isTop forceconnection
                scafR(indComp1,3)= scafR(indComp1,3)-CompV1 ;   %compensate  isBottom belong forceconnection
                %              scafR(isTop2A,3)= scafR(isTop2A,3)+obj.ssOption.OverRoute   ;   %compensate
                %              scafR(isBottom2A,3)= scafR(isBottom2A,3)-obj.ssOption.OverRoute   ;   %compensate
                %--------------
                scafR2= scafR;
                for updateSS=1:length(obj.ssOption.ForceConnUpdate)
                    %               PM2A=  obj.ssOption.ForceConnUpdate{updateSS}.BCBPM2(3)-obj.ssOption.ForceConnUpdate{updateSS}.BCB0(3);
                    %               PM2B=  obj.ssOption.ForceConnUpdate{updateSS}.BCBPM2(6)-obj.ssOption.ForceConnUpdate{updateSS}.BCB0(6);
                    %               PM2A=0;PM2B=0;
                    OldsPM2= obj.ssOption.ForceConnUpdate{updateSS}.BCBPM2 ;
                    NewEnds=obj.ssOption.ForceConnUpdate{updateSS}.BCBCur ;
                    
                    [~,iiid1]= ismember(Scaf0,OldsPM2(1:3),'rows') ; fd1=find(iiid1);
                    scafR2( fd1,:)=NewEnds(1:3)+scafR(fd1,:)-Scaf0(fd1,:)   ;
                    %              scafR2( fd1,:)=NewEnds(1:3);
                    [~,iiid2]= ismember(Scaf0,OldsPM2(4:6),'rows') ; fd2=find(iiid2);
                    scafR2(fd2,:)=NewEnds(4:6)+scafR(fd2,:)-Scaf0(fd2,:)    ;
                    %               scafR2(fd2,:)=NewEnds(4:6);
                end
                scafR3= scafR2;
                scafR= scafR2;
                %------------
                % %----------shift some scaf loop to "exact position" so that
                %  sticky ends can be used to polymerize
                % hard code
                %                 fprintf('sticky end option, hard code in obj.findCycleList \n' )
                if strcmp(obj.ScafOption.prevStack,'scafloop')
                    TargetPolyBundles= obj.ScafOption.PolyBundle ;
                    for Vedge = 2:2 :size(scafR2 ,1)-1
                        %                     [Vedge,size(scafR2,1)]
                        MM =   scafR2(Vedge:Vedge+1 ,:) ;
                        if ismember(MM(1,1) , TargetPolyBundles)
                            Bundle= obj.containBundle{MM(1,1)} ;
                            if MM(1,1) ==MM(2,1) && MM(1,3) ==MM(2,3) && (MM(1,3)<Bundle.Zbase1(MM(1,2)) || MM(1,3)>Bundle.Zbase2(MM(1,2)) ) % locate out of z1z2
                                
                                CylMaster=MM(1,2) ;
                                CylInSlave=MM(2,2) ;
                                CinMoveUp=  ~xor( ismember(CylMaster,Bundle.AGroup),Bundle.AGroupGoUp);
                                rnd= [MM(1,3)-10,MM(1,3),MM(1,3)+10]   ;XOinZ=zeros(1,3) ;
                                if strcmp(Bundle.Lattice, 'Square')            %change for hybrid structure
                                    XOinZ(1)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(1),CinMoveUp);
                                    XOinZ(2)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(2),CinMoveUp);
                                    XOinZ(3)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(3),CinMoveUp);
                                else
                                    %                                 XOinZ= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd,CinMoveUp);
                                    XOinZ(1)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(1),CinMoveUp);
                                    XOinZ(2)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(2),CinMoveUp);
                                    XOinZ(3)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(3),CinMoveUp);
                                end
                                
                                QQ = abs(XOinZ-MM(1,3));
                                choose= XOinZ(QQ==min(QQ)) ;
                                if choose<Bundle.Zbase1(MM(1,2))
                                    choose=choose+1;
                                end
                                scafR2(Vedge:Vedge+1 ,3) =choose ;
                            end
                        end
                    end
                    %-----------------
                    %
                    %                 BundleNoScafLoop = 10:17 ;    fprintf('BundleNoScafLoop = %s \n', num2str(BundleNoScafLoop)  ) ;
                    BundleNoScafLoop=obj.ScafOption.BundleNoScafLoop ;
                    for k =1:size(scafR2,1)
                        if ismember(scafR2(k,1) ,BundleNoScafLoop)
                            BB=scafR2(k,1) ;
                            if   scafR2(k,3) < obj.containBundle{BB}.Zbase1(scafR2(k,2))   % over-route on Z1
                                scafR2(k,3)=obj.containBundle{BB}.Zbase1(scafR2(k,2))  ;
                            elseif scafR2(k,3) > obj.containBundle{BB}.Zbase2(scafR2(k,2))   % over-route on Z1
                                scafR2(k,3)=obj.containBundle{BB}.Zbase2(scafR2(k,2))   ;
                            end
                        end
                    end
                    %---------
                elseif  strcmp(obj.ScafOption.prevStack,'polyT')
                    BundleNoScafLoop=obj.ScafOption.BundleNoScafLoop  ;
                    polyTBundle= setdiff(1:length(obj.containBundle) , BundleNoScafLoop) ;
                    for Vedge = 2:2 :size(scafR2 ,1)-1
                        MM =   scafR2(Vedge:Vedge+1 ,:) ;
                        if ismember(MM(1,1) , polyTBundle)
                            Bundle= obj.containBundle{MM(1,1)} ;
                            if MM(1,1) ==MM(2,1) && MM(1,3) ==MM(2,3) && (MM(1,3)<Bundle.Zbase1(MM(1,2)) || MM(1,3)>Bundle.Zbase2(MM(1,2)) ) % locate out of z1z2, exclude interal Xovers
                                
                                CylMaster=MM(1,2) ;
                                CylInSlave=MM(2,2) ;
                                CinMoveUp=  ~xor( ismember(CylMaster,Bundle.AGroup),Bundle.AGroupGoUp);
                                if MM(1,3)>Bundle.Zbase2(MM(1,2)) % Z2 side
                                    rnd= [MM(1,3)-20,MM(1,3)-10,MM(1,3)]   ;
                                else
                                    rnd= [MM(1,3),MM(1,3)+10,MM(1,3)+20]   ;
                                end
                                
                                XOinZ=zeros(1,3) ;
                                if strcmp(Bundle.Lattice, 'Square')            %change for hybrid structure
                                    XOinZ(1)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(1),CinMoveUp);
                                    XOinZ(2)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(2),CinMoveUp);
                                    XOinZ(3)= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(3),CinMoveUp);
                                else
                                    %                                 XOinZ= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd,CinMoveUp);
                                    XOinZ(1)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(1),CinMoveUp);
                                    XOinZ(2)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(2),CinMoveUp);
                                    XOinZ(3)= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd(3),CinMoveUp);
                                end
                                
                                if MM(1,3)>Bundle.Zbase2(MM(1,2)) % Z2 side
                                    QQ = XOinZ-Bundle.Zbase2(MM(1,2));
                                    QQ(QQ>0)=200;
                                    QQ=abs(QQ) ;
                                else   % Z1 side
                                    QQ = XOinZ-Bundle.Zbase1(MM(1,2));
                                    QQ(QQ<0)=200;
                                    QQ=abs(QQ) ;
                                end
                                choose= XOinZ(QQ==min(QQ)) ;
                                if MM(1,3)<Bundle.Zbase1(MM(1,2))
                                    choose=choose+1;
                                end
                                scafR2(Vedge:Vedge+1 ,3) =choose ;
                            end
                        end
                    end
                    BundleNoScafLoop=obj.ScafOption.BundleNoScafLoop ;
                    for k =1:size(scafR2,1)
                        if ismember(scafR2(k,1) ,BundleNoScafLoop)
                            BB=scafR2(k,1) ;
                            if   scafR2(k,3) < obj.containBundle{BB}.Zbase1(scafR2(k,2))   % over-route on Z1
                                scafR2(k,3)=obj.containBundle{BB}.Zbase1(scafR2(k,2))  ;
                            elseif scafR2(k,3) > obj.containBundle{BB}.Zbase2(scafR2(k,2))   % over-route on Z1
                                scafR2(k,3)=obj.containBundle{BB}.Zbase2(scafR2(k,2))   ;
                            end
                        end
                    end
                end
                
                %circular shift to longest edge
                edgwsL = abs(scafR2(2:2:end ,3) -scafR2(1:2:end ,3));
                Movestep=  find(edgwsL==max(edgwsL)) ;    Movestep=Movestep(1) ;
                scafR2=  circshift(scafR2,-2*Movestep+2) ;
                %                  [scafR,scafR2] ;
                %-----
                Vec=round((scafR2(2,:)-scafR2(1,:))/2);        Veci=Vec/norm(Vec);
                
                NewscafR2=[ scafR2(1,:);scafR2(1,:)+Vec;scafR2(1,:)+Vec+Veci; scafR2(2:end,:)];
                NewscafR2=circshift(NewscafR2,-2);
                %---------------------------------------------------------++++++++++++
                NOriGPairBCB{scaf_j}=NewscafR2 ;
                ScafR_temp =NOriGPairBCB{scaf_j} ;
                [ ~,whichCyl ] = CheckScafRXover(NOriGPairBCB{scaf_j}, tol  ,obj.ScafOption.BundleNoScafLoop) ;  % or add more bundle
                [~,iy1]=ismember( whichCyl, CombiningXover(:,1:3),'rows' ) ;
                [~,iy2]=ismember( whichCyl, CombiningXover(:,4:6),'rows' ) ;
                iyy= union(iy1(iy1~=0) ,iy2(iy2~=0) ) ;
                iyyEven= zeros( length(iyy) ,2) ;
                for ci= 1:length(iyy)
                    if mod(iyy(ci),2)==1
                        iyyEven(ci,:) = [iyy(ci) ,iyy(ci)+1];
                    else
                        iyyEven(ci,:) = [iyy(ci)-1 ,iyy(ci)];
                    end
                end
                iyyEven= unique(iyyEven,'rows') ;
                
                RemoveXovers =CombiningXover(iyyEven' ,:) ;
                
                nNotGood =size(whichCyl ,1)  ;   % if initial solution is too bad, just back to outer loop to get a better one
                if  nNotGood< 50
                    fprintf( ' Seems to be able to find a good scaf \n')
                    %                 [nNotGood ,  size(RemoveXovers,1)]
                    for tworound= 1:2
                        if tworound==1
                            SaveAppliedXoverInFirstRound = zeros(200, 6) ; kx =1;
                        else
                            fprintf('start second round \n') ;
                            if nNotGood==0 || isempty(RemoveXovers) ; break;end         % in case that first round is good
                            
                            RemoveXovers= SaveAppliedXoverInFirstRound(1:kx-1 , :) ; % second round try to remove newly added Xovers
                            IndRemain=   find( or( ismember( RemoveXovers(:,1:3), whichCyl2,'rows'),ismember( RemoveXovers(:,4:6), whichCyl2,'rows')) ) ;  % only reserve xover has problem
                            IndRemain= union(IndRemain(mod(IndRemain,2)==1), IndRemain(mod(IndRemain,2)==0)-1 )  ; IndRemain=union(IndRemain, IndRemain+1) ;
                            RemoveXovers=RemoveXovers( IndRemain,:) ;
                            %                             RemoveXovers
                        end
                        
                        for riXover=  1:2:size(RemoveXovers,1)
                            Xover = RemoveXovers(riXover:riXover+1 ,:) ;
                            TwoCycles=    removeScafXover_general(obj,NOriGPairBCB(scaf_j),Xover)  ;
                            [~, XoverList]=obj.Given2CylceFindXoverList(TwoCycles{1},TwoCycles{2},  []) ; % avoid in while loop. Query one time for a list, 09/11/2019
                            
                            if isempty(XoverList)
                                EEemptyHappened =1
                                continue;
                            end
                            cw= 1;
                            while 1
                                %                                                 OneXover=obj.Given2CylceFindXoverList(TwoCycles{1},TwoCycles{2},  CombiningXover,ConnectList   )
                                %                                 OneXover=obj.Given2CylceFindXoverList(TwoCycles{1},TwoCycles{2},  []  ) ;
                                OneXover = XoverList(  randi(size(XoverList,1),1) , :)  ;
                                OneXover= [OneXover(1,1:6) ;OneXover(1,7:12) ];
                                [NOriGPairBCBOut,~, ~] = obj.AddXover22cycles(TwoCycles,OneXover)     ;
                                [ ~,whichCyl2 ] = CheckScafRXover(NOriGPairBCBOut{1}, tol,obj.ScafOption.BundleNoScafLoop) ;
                                cw=cw+1;
                                if cw>30
                                    whichCyl2;
                                    fprintf('reach max loop\n');
                                    break;   %
                                end
                                if size(whichCyl2,1) < nNotGood
                                    nNotGood= size(whichCyl2,1) ;
                                    NOriGPairBCB{scaf_j}= NOriGPairBCBOut{1};
                                    %                             CombiningXover=[CombiningXover ; OneXover];
                                    SaveAppliedXoverInFirstRound(kx:kx+1,: )= OneXover ; kx=kx+2 ;
                                    break
                                end
                                if nNotGood==0
                                    break
                                end
                            end
                            %                     cw ;
                            fprintf('nCornerOnBadXover = %i\n',   size(whichCyl2,1))
                        end
                    end  % for tworound= 1:2
                    
                end
                [ Goodcc,whichCylcc ] = CheckScafRXover(NOriGPairBCB{1}, tol,obj.ScafOption.BundleNoScafLoop ) ; whichCylcc
                fprintf(' In findCyclist,   Goodcc = %i \n '  ,Goodcc)
                NewscafR2=  NOriGPairBCB{scaf_j} ;
                %----------
                obj.ScafRouting{scaf_j}=[NewscafR2  ] ;
            end
            %           obj.getbetterTMatrix
            %             profile viewer
            
        end  %end of findCycleList
        
        function Hubs=GetHubBundleInd(obj)
            IsHubs=zeros(size(obj.containBundle)) ;
            for k=1:length(IsHubs)
                if strcmp(obj.containBundle{k}.type , 'Hub')
                    IsHubs(k)=1 ;
                end
            end
            Hubs=find(IsHubs) ;
        end
        
        
        function TwoCycles=removeScafXover(obj,NOriGPairBCB,Xover)
            OneCycle = NOriGPairBCB{1} ;
            
            %             NOriGPairBCB
            %             Xover
            [~,ib] =ismember([Xover(:,1:3) ;Xover(:,4:6)] , NOriGPairBCB{1} ,'rows') ;
            inds =sort(ib) ;  %these four corners need to be removed.
            if inds(1)~=1
                SecondCycle= circshift(OneCycle(inds(2)+1 : inds(3)-1 ,:)  ,1);
                TwoCycles= { [OneCycle(1:inds(1)-1,:) ;OneCycle(inds(end)+1:end,: )] , SecondCycle } ;
            else  % happend to be the start
                SecondCycle= circshift(OneCycle(inds(3)+1 : end-1 ,:)  ,1);
                TwoCycles= { circshift(OneCycle(2:inds(2)-1,:),1)  , SecondCycle } ;
            end
        end % end of  function removeScafXover
        
        function OutputPairBCB=removeScafXover_general(obj,InputPairBCB,Xover)
            
            ListXoverPoints = [Xover(:,1:3) ; Xover(:,4:6) ];
            InvolvedCycles = [];
            for k =1: length(InputPairBCB)
                [AA,~] =ismember(ListXoverPoints ,InputPairBCB{k}, 'rows'  ) ;
                %             [AA,BB] =ismember( InputPairBCB{k},ListXoverPoints, 'rows'  )
                if sum(AA)>0
                    InvolvedCycles=union(InvolvedCycles , k) ;
                end
            end
            OutputPairBCB=InputPairBCB(setdiff(1:length(InputPairBCB),InvolvedCycles)  );
            if length(InvolvedCycles) ==2  % Xover on two cylcle, result should reduce by one.
                %             NewCylcle = [];
                [~,BB1] =ismember(ListXoverPoints ,InputPairBCB{InvolvedCycles(1)}, 'rows'  ) ;
                [~,BB2] =ismember(ListXoverPoints ,InputPairBCB{InvolvedCycles(2)}, 'rows'  ) ;
                BB1=BB1(BB1~=0) ;
                BB2=BB2(BB2~=0) ;
                NewCylcle= [ InputPairBCB{InvolvedCycles(1)}(1:BB1(1)-1,:) ;  InputPairBCB{InvolvedCycles(2)}(BB2(2)+1:end,:) ;...
                    InputPairBCB{InvolvedCycles(2)}(1:BB2(1)-1,:) ; InputPairBCB{InvolvedCycles(1)}(BB1(2)+1:end,:) ;];
                
                nw=0;
                while ~(NewCylcle(1,1)==NewCylcle(2,1) && NewCylcle(1,2)==NewCylcle(2,2))
                    NewCylcle =circshift(NewCylcle,1) ; nw=nw+1 ;
                    if nw>1000;break;end
                end
                
                
                OutputPairBCB{end+1} =NewCylcle;
            elseif length(InvolvedCycles) ==1  % Xover on one cylcle, result should increase by one.
                [~,BB] =ismember(ListXoverPoints ,InputPairBCB{InvolvedCycles}, 'rows'  ) ;
                
                nw=0;
                while max(BB)== size(InputPairBCB{InvolvedCycles} ,1)
                    InputPairBCB{InvolvedCycles} =circshift(InputPairBCB{InvolvedCycles},2) ;
                    [~,BB] =ismember(ListXoverPoints ,InputPairBCB{InvolvedCycles}, 'rows'  ) ;
                    nw=nw+1 ;
                    if nw>1000;break;end
                end
                
                BB=sort(BB) ;
                NewCylcle1 = [InputPairBCB{InvolvedCycles}(1:BB(1)-1,:) ; InputPairBCB{InvolvedCycles}(BB(4)+1:end,:)] ;
                NewCylcle2 = [InputPairBCB{InvolvedCycles}(BB(2)+1:BB(3)-1,:) ] ;
                
                nw=0;
                while ~(NewCylcle1(1,1)==NewCylcle1(2,1) && NewCylcle1(1,2)==NewCylcle1(2,2))
                    NewCylcle1 =circshift(NewCylcle1,1) ; nw=nw+1 ;
                    if nw>1000;break;end
                end
                
                nw=0;
                while ~(NewCylcle2(1,1)==NewCylcle2(2,1) && NewCylcle2(1,2)==NewCylcle2(2,2))
                    NewCylcle2 =circshift(NewCylcle2,1) ; nw=nw+1 ;
                    if nw>1000;break;end
                end
                
                OutputPairBCB{end+1} =NewCylcle1;
                OutputPairBCB{end+1} =NewCylcle2;
                %                 case2 =1;
            else
                %                 dssdf=3
            end
        end
        
        function getbetterTMatrix(obj)
            
            scafR=obj.ScafRouting;
            NodePosition=zeros(size(scafR,1),3);
            for k=1:size(scafR,1)
                
                Acc=findHelixQ(obj.containBundle{scafR(k,1)},scafR(k,2),scafR(k,3));
                
                TT=obj.containBundle{scafR(k,1)}.SimulateTransMFromTM2;
                
                NodePosition(k,:)= (TT(1:3,1:3)*Acc{1}'+TT(1:3,4) )'  ;
            end
            OriNodePosition=NodePosition;
            
            SSDist=zeros(length(2:2:size(scafR,1)),1);
            for updateSS=2:2:size(scafR,1)-1
                SSDist(updateSS/2) = norm(NodePosition(updateSS,:)-NodePosition(updateSS+1,:));
            end
            %             OriObjValue=sum(SSDist>1)+sum(SSDist<0.6) ; % number of neiboring node out of bound
            %             OriObjValue=  10*std(SSDist)+1*sum(SSDist>1)+2*sum(SSDist<0.5);
            OriObjValue=         0.1*std(SSDist)+5*sum(SSDist>1.5)+20*sum(SSDist<0.5);
            nBun=max(scafR(:,1));
            
            
            nloopMax=10000; nloop=1;
            SaveObjVal=100000*ones(nloopMax,2) ; SaveObjVal(1,1)=OriObjValue; SaveObjVal(1,2)= 1*sum(SSDist>1.5)+ 2*sum(SSDist<0.5) ;
            
            CurrentNodePosition=NodePosition;
            BestNodePosit=NodePosition;
            tic
            while nloop<nloopMax
                
                NodePosition=     CurrentNodePosition;
                for TFi=1: nBun
                    NNN=sum(scafR(:,1)==TFi);
                    center=mean(NodePosition(scafR(:,1)==TFi,:));
                    NodePosition(scafR(:,1)==TFi,:)=NodePosition(scafR(:,1)==TFi,:)-ones(NNN,1)*center  ;
                    RMat= randomRotM(obj,0.05);
                    tVec=0.5*(rand(3,1)-0.5);
                    
                    NodePosition(scafR(:,1)==TFi,:)=(RMat*NodePosition(scafR(:,1)==TFi,:)'+ tVec*ones(1,  NNN )   )';
                    NodePosition(scafR(:,1)==TFi,:)=NodePosition(scafR(:,1)==TFi,:)+ones(NNN,1)*center;
                end
                
                SSDist=zeros(length(2:2:size(scafR,1)),1);
                for updateSS2=2:2:size(scafR,1)-1
                    SSDist(updateSS2/2) = norm(NodePosition(updateSS2,:)-NodePosition(updateSS2+1,:));
                end
                %                 ObjValue=sum(SSDist>1)+sum(SSDist<0.6) ; % number of neiboring node out of bound
                ObjValue=         0.1*std(SSDist)+5*sum(SSDist>1.5)+20*sum(SSDist<0.5);
                if ObjValue<min(SaveObjVal(:,1))
                    CurrentNodePosition=  NodePosition;
                    BestNodePosit=NodePosition;
                    
                end
                
                SaveObjVal(nloop+1,1)=ObjValue;
                SaveObjVal(nloop+1,2)= 1*sum(SSDist>1.5)+ 2*sum(SSDist<0.5) ;
                nloop=nloop+1  ;
            end
            toc
            figure(1234);clf;subplot(2,1,1)
            plot(1:nloop,SaveObjVal(1:nloop,1))
            subplot(2,1,2)
            plot(1:nloop,SaveObjVal(1:nloop,2))
            
            BestNodePosit;
            %----harvest Add new TF
            for TFi=1: nBun
                IndInvole=scafR(:,1)==TFi;
                Ori=OriNodePosition(IndInvole,:);
                Best= BestNodePosit(IndInvole,:);
                
                [regParams,~,Acc2]=absor(Ori',Best');
                %                  obj.containBundle{TFi}.SimulateTransMFromTM2= regParams.M*obj.containBundle{TFi}.SimulateTransMFromTM2;
            end
            figure(120); clf; hold on;
            SaveGHelix=cell(1,length(obj.containBundle));
            for Bundlei=1:length(obj.containBundle)
                QQWithPM10Bases =obj.containBundle{Bundlei}.HelixXYZG;
                SaveGHelix{Bundlei}= QQWithPM10Bases;
            end
            MaxBase=max(scafR(:,3));
            skipPattern1=10:60:MaxBase;   % skip mod2 =0   %test 3/2
            skipPattern2=40:60:MaxBase;
            XYZforBundles=cell(max(scafR(:,1)),1) ;
            OriNodePosition ; scafR  ;
            
            
            
            for edgeinSCR=1:2:size(scafR,1)-1
                %     edgeinSCR
                bundle=scafR(edgeinSCR,1);  Cyl=scafR(edgeinSCR+1,2);
                BaseStart=scafR(edgeinSCR,3);   BaseEnd=scafR(edgeinSCR+1,3);
                RelativeBS=BaseStart-obj.containBundle{bundle}.Zbase1(Cyl)+11;
                RelativeBE=BaseEnd-obj.containBundle{bundle}.Zbase1(Cyl)+11;
                QQ=linspace(RelativeBS,RelativeBE,abs(RelativeBE-RelativeBS)+1);
                
                if BaseStart>BaseEnd  %go to left------------------
                    QQ=setdiff(QQ, skipPattern2-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable');
                else
                    QQ=setdiff(QQ, skipPattern1-obj.containBundle{bundle}.Zbase1(Cyl)+11,'stable');
                end
                PartHelix=SaveGHelix{bundle}{Cyl}(QQ,:);
                TAll=obj.containBundle{bundle}.SimulateTransMFromTM2;
                QQQ=(TAll(1:3,1:3)*PartHelix'+TAll(1:3,4)*ones(1, size(PartHelix,1)  ))';
                scatter3(QQQ(:,1) ,QQQ(:,2), QQQ(:,3),'.')
                if isempty(  XYZforBundles{bundle})
                    XYZforBundles{bundle}=QQQ;
                else
                    XYZforBundles{bundle}= [ XYZforBundles{bundle} ;QQQ] ;
                end
            end
            
            
            
            
            
            %             scafR;
            %             BestNodePosit;
            
            figure(121); clf; hold on;
            nBase=0;
            alphSh=cell(  max(scafR(:,1)),1 ) ;
            for iB=1:length(alphSh);   nBase=nBase+   size(XYZforBundles{iB},1)  ;  end
            
            AllP=zeros(nBase,3); jAP=1;
            BundleCenters=zeros(length(alphSh),3);
            
            alphaBoundary=cell(  max(scafR(:,1)),1 ) ;
            for k=1:length(alphSh)
                points=  XYZforBundles{k};
                alphSh{k} = alphaShape(points(:,1),points(:,2),points(:,3),'HoleThreshold',2);
                %                  alphSh{k} = alphaShape(sortP(:,1),sortP(:,2),sortP(:,3),'HoleThreshold',2);
                alphSh{k}.Alpha=2 ;
                Bundlecenters=mean(points);
                plot( alphSh{k} );
                text(Bundlecenters(1),Bundlecenters(2),Bundlecenters(3)+5 ,num2str(k));
                
                AllP(jAP:size(points,1)+jAP-1,:)=points; jAP=jAP+size(points,1);
                BundleCenters(k,:)=mean(points);
                
                [~, xyz] = boundaryFacets( alphSh{k});
                testalpha=alphaShape(xyz(:,1),xyz(:,2),xyz(:,3),'HoleThreshold',2);
                testalpha.Alpha=2 ;
                
                alphaBoundary{k}=testalpha ;
            end
            
            SSaveIndex=zeros(size(OriNodePosition,1),3);
            for k=1:size(OriNodePosition,1)
                bund=scafR(k,1) ;
                [A,B]= ismember( OriNodePosition(k,:),  alphaBoundary{bund}.Points ,'rows');
                if A==0 && B==0
                    SSaveIndex(k,:)=[bund,1,1] ;      %internal Xover refers to same point so that calculating distance =0
                else
                    SSaveIndex(k,:)=[bund,A,B] ;
                end
            end
            
            
            
            GCenter=mean(AllP);
            scatter3(GCenter(1) ,GCenter(2),GCenter(3),'O');
            
            MoveUnitVector=BundleCenters - ones(size(BundleCenters,1) ,1)*GCenter ;
            for iU=1:size(MoveUnitVector,1)
                MoveUnitVector(iU,:)=    MoveUnitVector(iU,:) /norm(MoveUnitVector(iU,:));
            end
            
            stepsize=1;
            
            alphaBoundary;
            XYZforBundles2Sorface=cell(  max(scafR(:,1)),1 ) ;
            InitailInside=0;
            for checkinter=1:length(alphSh)
                otherBundle=setdiff(1:length(alphSh), checkinter) ;
                ZZZ=[];
                for kk=1:length(otherBundle)
                    if kk==1
                        %                           alphaBoundary{otherBundle(kk)}.Points
                        ZZZ= alphaBoundary{otherBundle(kk)}.Points;
                    else
                        ZZZ=[ZZZ; alphaBoundary{otherBundle(kk)}.Points];
                    end
                end
                
                tf = inShape( alphaBoundary{checkinter},ZZZ(:,1),ZZZ(:,2),ZZZ(:,3));
                XYZforBundles2Sorface{checkinter}=alphaBoundary{checkinter}.Points ;
                InitailInside=InitailInside+sum(tf) ;
            end
            %             CurrentXYZ=XYZforBundles;
            CurrentXYZ=XYZforBundles2Sorface;
            
            nloopMax=200; nloop=1;
            
            SSDist=zeros(length(2:2:size(scafR,1))-1,1);
            for updateSS2=2:2:size(scafR,1)-1
                PA=   CurrentXYZ{SSaveIndex(updateSS2,1)}(SSaveIndex(updateSS2,3) ,:) ;
                PB=   CurrentXYZ{SSaveIndex(updateSS2+1,1)}(SSaveIndex(updateSS2+1,3) ,:)   ;
                SSDist(updateSS2/2) = norm(PA-PB);
            end
            ObjValue=std( SSDist)   ;%+5*sum(SSDist>1.5);%  +20*sum(SSDist<0.5);
            SaveObjVal=100000*ones(nloopMax,2) ; SaveObjVal(1,1)=ObjValue; SaveObjVal(1,2)= 1*sum(SSDist>1.5)+ 2*sum(SSDist<0.5) ;
            tic
            
            nnInNewXYZ=zeros(1,length(alphSh));
            for ii=1:length(alphSh); nnInNewXYZ(ii)=size(CurrentXYZ{ii},1);end;
            
            while  nloop<nloopMax
                %                   tic
                NewXYZ=cell(1,length(alphSh));
                for update=1:length(alphSh)
                    PPXYZ= CurrentXYZ{update};
                    MovingBcent=mean(PPXYZ);
                    RMat= randomRotM(obj,sqrt(ObjValue));   % angle amplitutte=2 deg
                    tVec=1*(rand(1,3)-0.5);      % displacement amplitutte=1
                    NewXYZ{update}= (PPXYZ- ones(size(PPXYZ,1) ,1)*MovingBcent)* RMat + ones(size(PPXYZ,1) ,1)*tVec+ ones(size(PPXYZ,1) ,1)*MovingBcent;%   - 0.2*ones(size(PPXYZ,1) ,1)* MoveUnitVector(update,:)      ;
                    %                   XYZforBundles{update} = XYZforBundles{update}+stepsize*ones(size(XYZforBundles{update},1) ,1)* MoveUnitVector(update,:) ;
                    %                   alphSh{update}.Points = XYZforBundles{update} ;
                end
                %               toc
                %               tic
                %
                BundleInterfece=zeros(1,length(alphSh));
                %             NumberInside=0;
                alphaBoundary2=cell(1,length(alphSh));
                for checkinter=1:length(alphSh)
                    otherBundle=setdiff(1:length(alphSh), checkinter) ;
                    ZZZ=zeros(sum(nnInNewXYZ(otherBundle)),3); nZ=1;
                    %                  ZZZ=[];
                    %                  toc
                    %                  tic
                    for kk=1:length(otherBundle)
                        ZZZ(nZ:size(NewXYZ{otherBundle(kk)},1)+nZ-1,:)=NewXYZ{otherBundle(kk)};
                        nZ=nZ+size(NewXYZ{otherBundle(kk)},1);
                        
                    end
                    
                    %                 alphaBoundary{checkinter}.Points= NewXYZ{checkinter};
                    alphaBoundary2{checkinter}=alphaShape(NewXYZ{checkinter}(:,1),NewXYZ{checkinter}(:,2),NewXYZ{checkinter}(:,3),'HoleThreshold',2);
                    alphaBoundary2{checkinter}.Alpha=10 ;
                    tf = inShape( alphaBoundary2{checkinter},ZZZ(:,1),ZZZ(:,2),ZZZ(:,3));
                    % %                 NumberInside=NumberInside+sum(tf) ;
                    BundleInterfece(checkinter) =sum( tf)>0;
                    %                   toc
                    %                  tic
                end
                %               toc
                %               tic
                %
                SSDist=zeros(length(2:2:size(scafR,1))-1,1);
                for updateSS2=2:2:size(scafR,1)-1
                    PA=   NewXYZ{SSaveIndex(updateSS2,1)}(SSaveIndex(updateSS2,3) ,:) ;
                    PB=   NewXYZ{SSaveIndex(updateSS2+1,1)}(SSaveIndex(updateSS2+1,3) ,:)   ;
                    SSDist(updateSS2/2) = norm(PA-PB);
                end
                ObjValue=std( SSDist)  ;     %  0.1*std(SSDist)  ;% +5*sum(SSDist>1.5);%    +20*sum(SSDist<0.5);
                if ObjValue<min(SaveObjVal(:,1))  && sum(BundleInterfece)==0
                    CurrentXYZ=  NewXYZ;
                    %                 else
                    %                     for reI=1:length(alphSh)
                    %                       alphaBoundary{reI}.Points= CurrentXYZ{reI};
                    %                     end
                end
                %               [NumberInside, InitailInside ,nloop,ObjValue];
                SaveObjVal(nloop+1,1)=ObjValue;
                SaveObjVal(nloop+1,2)= 1*sum(SSDist>1.5)+ 2*sum(SSDist<0.5) ;
                nloop=nloop+1  ;
                %                 toc
                
                
            end   % end of while loop
            toc
            min(SaveObjVal(:,1))
            SSDist;
            figure(121); clf; hold on;
            for ploti=1:length(alphSh)
                plot( alphaBoundary2{ploti} );
            end
            sdfsff=234;
            
            
        end
        
        function R=randomRotM(obj,angRan)
            %angRan  in degree
            randAng= 2*angRan*rand(1,3)-[angRan,angRan,angRan];
            Rx=[1,0,0 ;0 ,cosd( randAng(1)), -sind(randAng(1)); 0, sind(randAng(1)), cosd( randAng(1))];
            Ry=[cosd( randAng(2)), 0, sind(randAng(2)); 0 , 1 , 0; -sind(randAng(2)),0 ,cosd( randAng(2))];
            Rz=[ cosd( randAng(3)),-sind(randAng(3)),0 ;sind(randAng(3)),  cosd( randAng(3)),0;0,0,1];
            R=Rx*Ry*Rz;
        end
        
        
        function OneXover=Given2CylceFindXoverList2(obj,Cycle1,Cycle2)
            %convert to C5 express
            NC1=zeros(size(Cycle1,1),2);
            NC2=zeros(size(Cycle2,1),2);
            for iNC1=1:size(NC1,1)
                [~,ind]=ismember(Cycle1(iNC1,1:2) ,  obj.RelateTable(:,1:2),'rows');
                NC1(iNC1,:)=[ obj.RelateTable(ind,5), Cycle1(iNC1,3)];
            end
            for iNC2=1:size(NC2,1)
                [~,ind]=ismember(Cycle2(iNC2,1:2) ,  obj.RelateTable(:,1:2),'rows');
                NC2(iNC2,:)=[ obj.RelateTable(ind,5), Cycle2(iNC2,3)];
            end
            %------------F-Convert   % to C5 express
            List=[];  OneXover=[];
            for edgeInNC1= 1:2:size(NC1,1)
                for edgeInNC2= 1:2:size(NC2,1)
                    C5cyl1=NC1(edgeInNC1,1);   C5cyl2=NC2(edgeInNC2,1);
                    ZRange1=sort([NC1(edgeInNC1:edgeInNC1+1,2)]) ;
                    ZRange2=sort([NC2(edgeInNC2:edgeInNC2+1,2)]) ;
                    if obj.CylinderC5AdjM(C5cyl1,C5cyl2)==1 && Cycle1(edgeInNC1,1)==Cycle2(edgeInNC2,1)
                        Bundle= obj.containBundle{Cycle1(edgeInNC1,1)};
                        ChooseBundle=Cycle1(edgeInNC1,1);
                        
                        if ChooseBundle~= Cycle2(edgeInNC2,1);
                            sdfsf=1244;
                        end
                        
                        CylMaster=Cycle1(edgeInNC1,2) ;
                        CylInSlave=Cycle2(edgeInNC2,2) ;
                        Intersec=intersect(ZRange1(1):ZRange1(2),ZRange2(1):ZRange2(2));
                        CinMoveUp=  ~xor( ismember(CylMaster,Bundle.AGroup),Bundle.AGroupGoUp);
                        %                        bound=5 ;
                        for eee=1:length(Intersec)
                            %                        if  ~isempty(Intersec)  &&  Intersec(end)-Intersec(1)>2*bound
                            %                        rnd=randi([Intersec(1)+bound ,Intersec(end)-bound]);
                            rnd=Intersec(eee);
                            XOinZ= FindZ2( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd,CinMoveUp,[],[] );
                            if XOinZ<max(ZRange2) && XOinZ> min(ZRange2)
                                XOinZ2=  XOinZ -1+2*CinMoveUp;
                                if CinMoveUp==0   %  debug result
                                    XOinZ2=XOinZ2+1;
                                    XOinZ=XOinZ+1;
                                end
                                OneXover=[ChooseBundle,CylMaster,XOinZ,ChooseBundle,CylInSlave,  XOinZ,...
                                    ChooseBundle,CylInSlave,XOinZ2,ChooseBundle,CylMaster,  XOinZ2]  ;
                                List=[ List; OneXover];
                                %                         return
                            end
                        end
                    end
                end
            end
            
            AA=unique(List,'rows');
            GoodAA=zeros(size(AA,1),1);
            for filteri=1:length(GoodAA)
                subAA= AA(filteri,:);
                Bundle=obj.containBundle{AA(filteri,1)};
                TwoCyls=unique(subAA([2, 5 ,8, 11]));
                Limits= [Bundle.Zbase1(TwoCyls); Bundle.Zbase2(TwoCyls)];
                Zavg=mean(subAA([3, 6 ,9, 12]));
                if Zavg- max(Limits(1,:)) <8 || -Zavg + min(Limits(2,:)) <8
                    GoodAA(filteri)=1;
                end
            end
            
            AA(GoodAA==1,:)=[];
            if ~isempty(AA)
                BB=AA(randi(size(AA,1)),:);
                OneXover=  [ BB(1:6); BB(7:12) ];
                TestTwoCycle=[{Cycle1}; {Cycle2}];
                [~, UnUsedXover,~]=AddXover22cycles(obj,TestTwoCycle,OneXover);
                nW=1;
                while   ~isempty(UnUsedXover)    % check this Xover is valid
                    BB=AA(randi(size(AA,1)),:);
                    OneXover=  [ BB(1:6); BB(7:12) ];
                    [~, UnUsedXover2,~]=AddXover22cycles(obj,TestTwoCycle,OneXover);
                    nW=nW+1;
                    if nW>=100
                        OneXover=[];
                        return
                    end
                end
            end
            %           OneXover=OneXover2;
        end  %end of function Given2CylceFindXoverList2
        
        function [OneXover2,varargout]=Given2CylceFindXoverList(obj,varargin)
            %              [s,varargout]  OneXover2
            %convert to C5 express
            Cycle1= varargin{1};
            Cycle2= varargin{2}  ;
            %            ~isempty(varargin{3})
            %             nargin;
            if  length(varargin)>2
                if ~isempty(varargin{3})
                    CombiningXover= varargin{3}  ;
                    CombiningXoverAA= unique([CombiningXover(1:2:end,:),CombiningXover(2:2:end,:)],'rows') ;
                    CombiningXoverBB= unique([CombiningXoverAA(:,1:3);CombiningXoverAA(:,4:6);CombiningXoverAA(:,7:9);CombiningXoverAA(:,10:12)],'rows') ;
                    CombiningXoverBB=setdiff(CombiningXoverBB,[0 0 0],'rows' );
                    
                    ForceConn=varargin{4}  ;
                    CombiningXoverBB=[CombiningXoverBB;ForceConn(:,1:3);ForceConn(:,4:6)];
                    
                else
                    CombiningXover=[];
                    CombiningXoverAA=[];
                    CombiningXoverBB=[];
                    ForceConn=[];
                end
            end
            
            NC1=zeros(size(Cycle1,1),2);
            NC2=zeros(size(Cycle2,1),2);
            for iNC1=1:size(NC1,1)
                [~,ind]=ismember(Cycle1(iNC1,1:2) ,  obj.RelateTable(:,1:2),'rows');
                NC1(iNC1,:)=[ obj.RelateTable(ind,5), Cycle1(iNC1,3)];
            end
            for iNC2=1:size(NC2,1)
                [~,ind]=ismember(Cycle2(iNC2,1:2) ,  obj.RelateTable(:,1:2),'rows');
                NC2(iNC2,:)=[ obj.RelateTable(ind,5), Cycle2(iNC2,3)];
            end
            %------------F-Convert   % to C5 express
            %             List=[];
            OneXover=[];
            
            List=zeros(900000, 12 ) ; ccL =1 ;
            for edgeInNC1= 1:2:size(NC1,1)-1
                for edgeInNC2= 1:2:size(NC2,1)-1
                    C5cyl1=NC1(edgeInNC1,1)   ;
                    C5cyl2=NC2(edgeInNC2,1);
                    %                     if C5cyl2==34 && C5cyl1==111
                    %                         sdfsf=3
                    %                     end
                    
                    ZRange1=sort([NC1(edgeInNC1:edgeInNC1+1,2)]) ;
                    %                                        [size(NC2), edgeInNC2]
                    
                    ZRange2=sort([NC2(edgeInNC2:edgeInNC2+1,2)]) ;
                    if obj.CylinderC5AdjM(C5cyl1,C5cyl2)==1 && Cycle1(edgeInNC1,1)==Cycle2(edgeInNC2,1)  %in same bundle and adjacent cylinders
                        Bundle= obj.containBundle{Cycle1(edgeInNC1,1)};
                        ChooseBundle=Cycle1(edgeInNC1,1);
                        
                        %                        if ChooseBundle~= Cycle2(edgeInNC2,1)
                        %                        sdfsf=1244;
                        %                        end
                        
                        CylMaster=Cycle1(edgeInNC1,2) ;
                        CylInSlave=Cycle2(edgeInNC2,2) ;
                        Intersec=intersect(ZRange1(1):ZRange1(2),ZRange2(1):ZRange2(2));
                        CinMoveUp=  ~xor( ismember(CylMaster,Bundle.AGroup),Bundle.AGroupGoUp);
                        %                        bound=5 ;
                        for eee=1:length(Intersec)
                            %                        if  ~isempty(Intersec)  &&  Intersec(end)-Intersec(1)>2*bound
                            %                        rnd=randi([Intersec(1)+bound ,Intersec(end)-bound]);
                            rnd=Intersec(eee);
                            if strcmp(Bundle.Lattice, 'Square')            %change for hybrid structure
                                XOinZ= FindZ2SQ( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd,CinMoveUp);
                            else
                                XOinZ= FindZ2HC( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd,CinMoveUp);
                            end
                            
                            if XOinZ<max(Intersec) && XOinZ> min(Intersec)
                                XOinZ2=  XOinZ -1+2*CinMoveUp;
                                if CinMoveUp==0   %  debug result
                                    XOinZ2=XOinZ2+1;
                                    XOinZ=XOinZ+1;
                                end
                                OneXover=[ChooseBundle,CylMaster,XOinZ,ChooseBundle,CylInSlave,  XOinZ,...
                                    ChooseBundle,CylInSlave,XOinZ2,ChooseBundle,CylMaster,  XOinZ2]  ;
                                %                                 List=[ List; OneXover];
                                List(ccL,:) =OneXover ;
                                ccL=ccL+1 ;
                                %                                 %                         return
                                %                                 if sum(OneXover==[20     3   124    20     4   124    20     4   123    20     3   123])==12
                                %                                     sdf=3
                                %                                     edgeInNC1
                                %                                     edgeInNC2
                                %                                 end
                            end
                        end
                    end
                end
            end
            List=List(1:ccL-1 ,:) ;
            
            tol=2  ; %%scaffold xover clearance from edges
%             tol=8  ;
            AA=unique(List,'rows');
            BadAA=zeros(size(AA,1),1);
%             IDontWantXoverOn = []  ; % hard code, exclude crossover on these bundles.
            IDontWantXoverOn = obj.ScafOption.NoScafXover ;
            for filteri=1:length(BadAA)
                subAA= AA(filteri,:);
                Bundle=obj.containBundle{AA(filteri,1)};
                TwoCyls=unique(subAA([2, 5 ,8, 11]));
                Limits= [Bundle.Zbase1(TwoCyls); Bundle.Zbase2(TwoCyls)];
                Zavg=mean(subAA([3, 6 ,9, 12]));
                if Zavg- max(Limits(1,:)) <tol || -Zavg + min(Limits(2,:)) <tol   %xover close to ends
                    BadAA(filteri)=1;
                    continue
                end
                
                CheckA=[subAA(1:2), (subAA(6)+subAA(9))/2] ;
                CheckB=[subAA(4:5), (subAA(6)+subAA(9))/2] ;
                
                for f2=1:size(CombiningXoverBB,1)          % potential xover compare with existing Xovers
                    if size(unique([CheckA(1:2) ;  CombiningXoverBB(f2,1:2)],'rows'),1)==1 && abs(CheckA(3)-CombiningXoverBB(f2,3))<tol
                        BadAA(filteri)=1; break;
                    elseif size(unique([CheckB(1:2) ;  CombiningXoverBB(f2,1:2)],'rows'),1)==1 && abs(CheckB(3)-CombiningXoverBB(f2,3))<tol
                        BadAA(filteri)=1; break;
                    end
                end
                
                if ismember(subAA(1) , IDontWantXoverOn)
                    BadAA(filteri)=1;
                end
                
            end
            
            
            
            OneXover2=[];
            AA(BadAA==1,:)=[];
            
            %+----- take some Xovers out for multi-scaf algorithm, block some adjc.
            if ~isempty(AA)
                
                BadAA2 = ismember( AA(:,[1,2,4,5]),obj.OrgBlockAdj ,'rows'   ) ;
                AA(BadAA2==1,:)=[];
            end
            if ~isempty(AA)
                
                BadAA3 = ismember( AA(:,[4,5,1,2]),obj.OrgBlockAdj ,'rows'   ) ;
                AA(BadAA3==1,:)=[];
            end
            %-----------  take some Xovers out for multi-scaf algorithm, block some Z range in each bundle.
            if ~isempty(obj.SaveBlockZRange) && ~isempty(AA)
                
                AllBlockZRange = unique( vertcat(obj.SaveBlockZRange{:}) ,'rows') ;  % [bundle, Zleft, Zright  ]
                AAO =AA ;
                IndBB=zeros(size(AA,1) ,1) ;
                for ka = 1 : length(IndBB)
                    Zmean  =mean(AA(ka,[3,6,9,12]) ,2) ;
                    IndBBAll  = and( ismember(AllBlockZRange(:,1),AA(ka,1) ) ,   and(Zmean> AllBlockZRange(:,2) ,Zmean< AllBlockZRange(:,3) )    )  ;
                    IndBB(ka) =  sum(IndBBAll)>0 ;
                end
                
                AA(IndBB==1,:)=[];
                
                
                %             sdfsdf=3
            end
            %---------
            
            if ~isempty(AA)
                varargout{1} =AA ;
                BB=AA(randi(size(AA,1)),:);
                OneXover2=  [ BB(1:6); BB(7:12) ];
                TestTwoCycle=[{Cycle1}; {Cycle2}];
                [~, UnUsedXover,~]=AddXover22cycles(obj,TestTwoCycle,OneXover2);
                nW=1;
                while   ~isempty(UnUsedXover)    % check this Xover is valid
                    BB=AA(randi(size(AA,1)),:);
                    OneXover2=  [ BB(1:6); BB(7:12) ];
                    [~, UnUsedXover,~]=AddXover22cycles(obj,TestTwoCycle,OneXover2);
                    nW=nW+1;
                    if nW>=100
                        OneXover2=[];
                        return
                    end
                end
            else
                varargout{1} =[];
            end
            if  isempty(OneXover2)
                fsdf=3 ;
            end
            
            
        end  %end of function Given2CylceFindXoverList
        
        function GetOrgBlockAdj(obj)
            OrgBlockAdjQQ= [-1,-1,-1,-1] ;
            
            if ~isempty(obj.SaveBlockAdjacency)
                for k = 1 : length( obj.containBundle)
                    if          ~isempty(obj.SaveBlockAdjacency{k})
                        MM =  obj.SaveBlockAdjacency{k} ;
                        OrgBlockAdjQQ=[OrgBlockAdjQQ ; [k*ones(size(MM,1) ,1) , MM(:,1)  ,k*ones(size(MM,1) ,1) , MM(:,2)   ] ] ;
                    end
                end
            end
            OrgBlockAdjQQ= setdiff(OrgBlockAdjQQ ,[-1,-1,-1,-1],'rows') ;
            
            obj.OrgBlockAdj = OrgBlockAdjQQ ;
            %             OrgBlockAdj
        end
        
        %         function OneXover=Given2CylceFindXoverListCOPY(obj,Cycle1,Cycle2)
        %             %determin in which bundle
        %             ShoedBundle=[Cycle1(:,1);Cycle2(:,1)];
        %             [a,b]=hist(ShoedBundle,unique(ShoedBundle));
        %             ChooseBundle=b(a==max(a));ChooseBundle=ChooseBundle(1);
        %             Bundle=obj.containBundle{ChooseBundle};
        %             NC1= Cycle1(Cycle1(:,1)==ChooseBundle,2:3);  %master
        %             NC2= Cycle2(Cycle2(:,1)==ChooseBundle,2:3);  %slave
        %
        %             NC2AllCyl=unique(NC2(:,1));
        %             OneXover=[]; List=[];
        %             for edgeiInNC1=1:2:size(NC1,1)
        %                 CylMaster= NC1(edgeiInNC1);
        %                 TheyConnect= Bundle.CylAdjMat(CylMaster,NC2AllCyl) ;
        %                 KK=find(TheyConnect==1);
        %                 if sum(TheyConnect)>0    %find neighboring cylinders
        %                     for slaveCyl=1:length(KK)
        %                         CylInSlave=NC2AllCyl(KK(slaveCyl)) ;
        %                         ZrangeInCyl=sort([NC1(edgeiInNC1,2) , NC1(edgeiInNC1+1,2)]);
        %                         ZrangeInCylSlave= sort(NC2(NC2(:,1)==CylInSlave,2))' ;
        %
        %                         Intersec=intersect(ZrangeInCyl(1):ZrangeInCyl(2),ZrangeInCylSlave(1):ZrangeInCylSlave(2));
        %                         Intersec=sort(Intersec);
        %                         if ~isempty(Intersec)
        %                             for repeat=1:5
        %                                 bound=5;
        %                                 Intersec;
        %                                 rnd=randi([Intersec(1)+bound ,Intersec(end)-bound]);
        %                                  CinMoveUp=  ~xor( ismember(CylMaster,Bundle.AGroup),Bundle.AGroupGoUp);
        %                                 XOinZ= FindZ2( CylMaster,CylInSlave,Bundle.CylInplanePosition,rnd,CinMoveUp,[],[] );
        %                                 if XOinZ<max(ZrangeInCylSlave) && XOinZ> min(ZrangeInCylSlave)
        %                                  XOinZ2=  XOinZ -1+2*CinMoveUp;
        %                                  if CinMoveUp==0   %  debug result
        %                                    XOinZ2=XOinZ2+1;
        %                                    XOinZ=XOinZ+1;
        %                                  end
        %                                 OneXover=[ChooseBundle,CylMaster,XOinZ,ChooseBundle,CylInSlave,  XOinZ,...
        %                                           ChooseBundle,CylInSlave,XOinZ2,ChooseBundle,CylMaster,  XOinZ2]  ;
        %                                 List=[ List; OneXover];
        %         %                         return
        %                                 end
        %                             end
        %                         end
        %                     end
        %                 end
        %             end
        %
        %           AA=unique(List,'rows');
        %           if ~isempty(AA)
        %           BB=AA(randi(size(AA,1)),:);
        %           OneXover=  [ BB(1:6); BB(7:12) ];
        %           end
        %         end  %end of function Given2CylceFindXoverList
        
        
        function NOriGPairBCB=AddXover2OneCycle(obj,OriGPairBCB,Xover)
            InputCell=OriGPairBCB;
            if isempty(Xover)
                NOriGPairBCB=OriGPairBCB;
                return
            end
            FourPoints=[Xover(:,1:3);Xover(:,4:6)];
            FourPoints=sortrows(FourPoints);
            whichcell=zeros(1,4);
            n=1;
            for iP=1:4   %search for which two cell contains these four points
                PP=FourPoints(iP,:);
                for jCell=1:length(OriGPairBCB)
                    Cycle=OriGPairBCB{jCell};
                    if ismember(PP,Cycle,'rows')
                        whichcell(n)= jCell;
                        n=n+1;
                    else
                        
                        Mix=sortrows([Cycle;PP]);
                        [~,indexofPP]= ismember(PP,Mix,'rows');
                        if mod(indexofPP,2)==0  && indexofPP~=1 && indexofPP~=size(Mix,1)
                            if size(union(  Mix(indexofPP-1:indexofPP+1, 1:2),[-1 -1 ],'rows') ,1)==2 %check on same edge
                                if (Mix(indexofPP,3)-Mix(indexofPP+1,3))*(Mix(indexofPP,3)-Mix(indexofPP-1,3))<0  %z value between
                                    whichcell(n)= jCell;
                                    n=n+1;
                                end
                            end
                        end
                    end
                end
            end
            whichcell=setdiff(whichcell,0);
            takeU=unique(whichcell);
            NumberOfCylceInvolved= length(  takeU  )    ;
            if NumberOfCylceInvolved==1
                CellA=OriGPairBCB{takeU};
                %----------Optimize performance 09/11/2019, base by base
                %sort
                BaseByBase_BCB = interpolateBase_ThreeCol( CellA ) ;
                [AA,BB] = ismember(FourPoints ,BaseByBase_BCB  ,'rows' )       ;
                
                [a,b] = sortrows(BB) ;
                FourPoints=FourPoints(b,:)  ;
                BB=BB(b,:)  ;
                
                
                [AA2,BB2] = ismember(CellA ,  BaseByBase_BCB ,'rows' ) ;
                
                CC1_new= [ BaseByBase_BCB( BB2(BB2<min(BB)) ,:) ;  FourPoints(1,: ) ; FourPoints(4,: ) ;  BaseByBase_BCB( BB2(BB2>max(BB)) ,:) ] ;
                CC2_new= [  FourPoints(2,: ) ; BaseByBase_BCB( BB2( and(BB2<max(BB),BB2>min(BB) )   ) ,:)  ;  FourPoints(3,: ) ] ;
                
                CC1_new= unique(CC1_new,'rows','stable');
                CC2_new= unique(CC2_new,'rows','stable');
                
                InputCell(takeU)=[];
                InputCell{end+1}=CC1_new;
                InputCell{end+1}=CC2_new;
                NOriGPairBCB=InputCell;
                %                 UseNew =1
                return
                %-----------------------
                
                %                 CheckPoint=Xover(:,4:6);    %--------
                InsertPosition=zeros(2,2);
                for iraw=1:size(CellA,1)-1
                    MM1=[CellA(iraw:iraw+1,:);Xover(1,4:6)];
                    MM2=[CellA(iraw:iraw+1,:);Xover(2,4:6)];
                    NN1=[CellA(iraw:iraw+1,:);Xover(1,1:3)];
                    NN2=[CellA(iraw:iraw+1,:);Xover(2,1:3)];
                    
                    if size(unique(MM1(:,1:2),'rows') ,1)==1  && (MM1(1,3)-MM1(3,3))*(MM1(2,3)-MM1(3,3))<0
                        InsertPosition(1,2)= iraw+0.5;
                    end
                    if size(unique(MM2(:,1:2),'rows') ,1)==1  && (MM2(1,3)-MM2(3,3))*(MM2(2,3)-MM2(3,3))<0
                        InsertPosition(2,2)= iraw+0.5;
                    end
                    if size(unique(NN1(:,1:2),'rows') ,1)==1  && (NN1(1,3)-NN1(3,3))*(NN1(2,3)-NN1(3,3))<0
                        InsertPosition(1,1)= iraw+0.5;
                    end
                    if size(unique(NN2(:,1:2),'rows') ,1)==1  && (NN2(1,3)-NN2(3,3))*(NN2(2,3)-NN2(3,3))<0
                        InsertPosition(2,1)= iraw+0.5 ;
                    end
                    
                end
                for checkside=1:2
                    for prelate=1:2
                        [A,ind]= ismember(Xover(checkside,3*prelate-2:3*prelate),CellA,'rows');
                        if A==1
                            InsertPosition(checkside,prelate)=ind;
                        end
                    end
                end
                Xover;
                XoverOnside=zeros(2,2);
                XoverV2=-1*ones(size(Xover));
                for qq=1:2
                    for ww=1:3:4
                        if  ismember(Xover(qq,ww:ww+2),CellA,'rows')
                            XoverOnside(qq,(ww+2)/3)=1;
                        else
                            XoverV2(qq,ww:ww+2) =  XoverV2(qq,ww:ww+2);
                        end
                    end
                end
                
                
                
                sdsf=23424;
                
                sortAllIP=sort(InsertPosition(:));
                halfIndex= mod(sortAllIP,1)==0.5;
                Add=union( ceil(sortAllIP(halfIndex)),floor(sortAllIP(halfIndex)));
                sortAllIP=sortAllIP(~halfIndex);
                sortAllIP=union(sortAllIP,Add);
                
                XoverCoverRange= abs(InsertPosition(:,1)-InsertPosition(:,2));
                %                 XoverCoverRange(1) < XoverCoverRange(2)  %
                if                   InsertPosition(1,2) < InsertPosition(2,2)
                    %                     C1=[ CellA(1:InsertPosition(1,2),:) ; Xover(1,4:6) ; Xover(1,1:3) ; CellA(InsertPosition(2,2)+1:end,:)] ;
                    %                     C2=[  Xover(2,1:3); CellA(InsertPosition(1,2)+1:InsertPosition(2,2),:) ;   Xover(2,4:6)];
                    
                    C1=[ CellA(1:sortAllIP(1),:) ; Xover(1,4:6) ; Xover(1,1:3) ; CellA(sortAllIP(4):end,:)] ;
                    
                    
                    C2=[  Xover(2,1:3); CellA( sortAllIP(2):sortAllIP(3),:) ;   Xover(2,4:6)];
                    UseFcnXTo1C1=1;
                else
                    %                      setdiff(  Xover(2,1:3),CellA(InsertPosition(1),:) ,'rows')
                    
                    
                    %                     C1=[ CellA(1:floor(InsertPosition(1,2)),:) ; Xover(1,4:6) ; Xover(1,1:3) ; CellA(ceil(InsertPosition(1,1)):end,:)] ;
                    %                     C2=[  Xover(2,1:3); CellA( ceil(InsertPosition(2,1)):floor(InsertPosition(2,2)),:) ;   Xover(2,4:6)];
                    %
                    %                     %----test
                    C1=[ CellA(1:sortAllIP(1),:) ; Xover(2,4:6) ; Xover(2,1:3) ; CellA(sortAllIP(4):end,:)] ;
                    C2=[  Xover(1,1:3); CellA( sortAllIP(2):sortAllIP(3),:) ;   Xover(1,4:6)];
                    %
                    %
                    
                    %                     seFcnXTo1C2=2;
                end
                CC1= unique(C1,'rows','stable');
                CC2= unique(C2,'rows','stable');
                if mod(size(CC1,1),2)==1 || mod(size(CC2,1),2)==1
                    sdfsf=234 ;
                end
                %                 XXX{1}=C1
                %                 obj.drawBCBdata2(InputCell{2},5)   obj.drawBCBdata2({CellA},5)
                
                InputCell(takeU)=[];
                InputCell{end+1}=CC1;   % weird
                InputCell{end+1}=CC2;
                
                CC1_new
                CC2_new
                ssdfs=3
            end
            NOriGPairBCB=InputCell;
        end   % End of func AddXover2OneCycle
        
        
        
        %         function NOriGPairBCB=AddXover2OneCycle(obj,OriGPairBCB,Xover)
        %             InputCell=OriGPairBCB;
        %             if isempty(Xover)
        %                 NOriGPairBCB=OriGPairBCB;
        %                 return
        %             end
        %             FourPoints=[Xover(:,1:3);Xover(:,4:6)];
        %             FourPoints=sortrows(FourPoints);
        %             whichcell=zeros(1,4);
        %             n=1;
        %             for iP=1:4   %search for which two cell contains these four points
        %                 PP=FourPoints(iP,:);
        %                 for jCell=1:length(OriGPairBCB)
        %                 Cycle=OriGPairBCB{jCell};
        %                    if ismember(PP,Cycle,'rows')
        %                    whichcell(n)= jCell;
        %                    n=n+1;
        %                    else
        %
        %                     Mix=sortrows([Cycle;PP]);
        %                     [~,indexofPP]= ismember(PP,Mix,'rows');
        %                         if mod(indexofPP,2)==0  && indexofPP~=1 && indexofPP~=size(Mix,1)
        %                              if size(union(  Mix(indexofPP-1:indexofPP+1, 1:2),[-1 -1 ],'rows') ,1)==2 %check on same edge
        %                                  if (Mix(indexofPP,3)-Mix(indexofPP+1,3))*(Mix(indexofPP,3)-Mix(indexofPP-1,3))<0  %z value between
        %                                     whichcell(n)= jCell;
        %                                     n=n+1;
        %                                  end
        %                              end
        %                         end
        %                     end
        %                 end
        %             end
        %             whichcell=setdiff(whichcell,0);
        %             takeU=unique(whichcell);
        %             NumberOfCylceInvolved= length(  takeU  )    ;
        %               if NumberOfCylceInvolved==1
        %                Xover;
        %            % special care of compound joint,  subs into over connection
        %             AllFCColumn=[obj.ForcedConnectList(1:2:end,:),obj.ForcedConnectList(2:2:end,:)];
        %             XoverColumn=[Xover(1,:),Xover(2,:)];
        %             [~,XoverThisTime]=ismember(XoverColumn,AllFCColumn,'rows');
        % %             VVCB=obj.CompJointInfo.BinRepofCJ(:,1);
        % %
        %             XoverIsSide=zeros(2,2);
        %             for xi=1:2
        %                 for yj=1:2
        %                  PointBCB= Xover(xi,3*yj-2:3*yj);
        %                     if PointBCB(3)==obj.containBundle{PointBCB(1)}.Zbase1(PointBCB(2))-obj.OverRoute
        %                      XoverIsSide( xi,  yj)=1;
        %                     elseif PointBCB(3)==obj.containBundle{PointBCB(1)}.Zbase2(PointBCB(2))+obj.OverRoute
        %                      XoverIsSide( xi,  yj)=1;
        %                     end
        %                 end
        %             end
        %
        %
        % %            if  XoverThisTime==0
        % %                sdfsf=34
        % %            end
        %             if  XoverThisTime~=0
        %             if min(find(VVCB==VVCB(XoverThisTime)))~=XoverThisTime && sum(sum(XoverIsSide))==4
        %                 WWW= AllFCColumn(find(VVCB==VVCB(XoverThisTime)),:);
        %                 CompoudJointNodes=unique([WWW(:,1:3);WWW(:,4:6);WWW(:,7:9);WWW(:,10:12)],'rows');
        %                 WWW2=[WWW(:,1:6);WWW(:,7:12)];
        %                 CellA= OriGPairBCB{takeU(1)};
        %                   CellA2= OriGPairBCB{takeU(1)};
        %                 PreCompJ= setdiff(find(VVCB==VVCB(XoverThisTime)), XoverThisTime);
        %                 PreCompJBCB=   AllFCColumn(PreCompJ,:);
        %                 PreCompJBCB=unique([PreCompJBCB(:,1:3);PreCompJBCB(:,4:6);PreCompJBCB(:,7:9);PreCompJBCB(:,10:12)],'rows');
        %                CirShiftStopHere= intersect(PreCompJBCB,FourPoints,'rows');
        %
        %
        %                [~, GA0(1)]=ismember(Xover(1,1:3) ,CellA,'rows');  [~, GA0(2)]=ismember(Xover(2,4:6) ,CellA,'rows');
        %                [~, GB0(1)]=ismember(Xover(2,1:3) ,CellA,'rows');  [~, GB0(2)]=ismember(Xover(1,4:6) ,CellA,'rows');
        %                nCC=1;
        %                 while nCC<200
        %                    CellA2=circshift(CellA2,1) ;
        %                     CellA2=circshift(CellA2,1) ;
        %                    [~, GA(1)]=ismember(Xover(1,1:3) ,CellA2,'rows');  [~, GA(2)]=ismember(Xover(2,4:6) ,CellA2,'rows');
        %                    [~, GB(1)]=ismember(Xover(2,1:3) ,CellA2,'rows');  [~, GB(2)]=ismember(Xover(1,4:6) ,CellA2,'rows');
        %                    [A,B]= ismember(CirShiftStopHere,CellA2,'rows') ;
        %                    B;
        %                    if ismember(CellA2(1,:),CirShiftStopHere,'rows')
        %                        sdsf=4324;
        %                        CCW=diff(sort([GA,GB]));
        %                        if (abs( diff(GA))==1 || abs( diff(GB))==1)   %%&& CellA2(1,1)==CellA2(2,1) && CellA2(1,2)==CellA2(2,2)
        %                            CCase=1;
        %                         break
        %                        elseif max([GA,GB])==size(CellA2,1)  && CCW(2)==1
        %                            CCase=2;
        %                         break
        %                        end
        %                    end
        %                     if nCC==195
        %                         CCase=2;
        %                     end
        %
        %
        %                    nCC=nCC+1;
        %                 end
        %                 if nCC>195
        %                     sdff=234;
        %                 end
        %
        % %                [~, GA(1)]=ismember(Xover(1,1:3) ,CellA,'rows');  [~, GA(2)]=ismember(Xover(2,4:6) ,CellA,'rows');
        % %                [~, GB(1)]=ismember(Xover(2,1:3) ,CellA,'rows');  [~, GB(2)]=ismember(Xover(1,4:6) ,CellA,'rows');
        %
        %                switch CCase
        %                    case 1
        %                        if  abs( diff(GA))==1
        %                            NNCyl1=CellA2(1:min(GA),:);
        %                            NNCyl2=CellA2(max(GA):end,:);
        %                        elseif  abs( diff(GB))==1
        %                            NNCyl1=CellA2(1:min(GB),:);
        %                            NNCyl2=CellA2(max(GB):end,:);
        %                        else
        %                            sdf=234;
        %                        end
        %                    case 2
        %                   NOriGPairBCB=InputCell;
        %                   return
        %                        sdfsf=234;
        %
        %                end
        %                InputCell(takeU)=[];
        %                InputCell{end+1}=NNCyl1;
        %                InputCell{end+1}=NNCyl2;
        %                NOriGPairBCB=InputCell;
        %                return
        %             end
        %            end
        %             %-------------
        %                  CellA=OriGPairBCB{takeU};
        %                  CheckPoint=Xover(:,4:6);    %--------
        %
        %                  InsertPosition=zeros(2,2);
        %                  for iraw=1:size(CellA,1)-1
        %                      MM1=[CellA(iraw:iraw+1,:);CheckPoint(1,:)];
        %                      MM2=[CellA(iraw:iraw+1,:);CheckPoint(2,:)];
        %                      NN1=[CellA(iraw:iraw+1,:);Xover(1,1:3)];
        %                      NN2=[CellA(iraw:iraw+1,:);Xover(2,1:3)];
        %
        %                      if size(unique(MM1(:,1:2),'rows') ,1)==1  && (MM1(1,3)-MM1(3,3))*(MM1(2,3)-MM1(3,3))<=0
        %                         InsertPosition(1,2)= iraw+0.5;
        %                      end
        %                      if size(unique(MM2(:,1:2),'rows') ,1)==1  && (MM2(1,3)-MM2(3,3))*(MM2(2,3)-MM2(3,3))<=0
        %                         InsertPosition(2,2)= iraw+0.5;
        %                      end
        %                      if size(unique(NN1(:,1:2),'rows') ,1)==1  && (NN1(1,3)-NN1(3,3))*(NN1(2,3)-NN1(3,3))<=0
        %                         InsertPosition(1,1)= iraw+0.5;
        %                      end
        %                      if size(unique(NN2(:,1:2),'rows') ,1)==1  && (NN2(1,3)-NN2(3,3))*(NN2(2,3)-NN2(3,3))<=0
        %                         InsertPosition(2,1)= iraw+0.5 ;
        %                      end
        %
        %                  end
        %                  for checkside=1:2
        %                      for prelate=1:2
        %                      [A,ind]= ismember(Xover(checkside,3*prelate-2:3*prelate),CellA,'rows');
        %                      if A==1
        %                          InsertPosition(checkside,prelate)=ind;
        %                      end
        %                      end
        %                  end
        %                  Xover;
        %                  XoverOnside=zeros(2,2);
        %                  XoverV2=-1*ones(size(Xover));
        %                  for qq=1:2
        %                      for ww=1:3:4
        %                          if  ismember(Xover(qq,ww:ww+2),CellA,'rows')
        %                              XoverOnside(qq,(ww+2)/3)=1;
        %                          else
        %                              XoverV2(qq,ww:ww+2) =  XoverV2(qq,ww:ww+2);
        %                          end
        %                      end
        %                  end
        %                  sortAllIP=sort(InsertPosition(:));
        %                 halfIndex= mod(sortAllIP,1)==0.5;
        %                 Add=union( ceil(sortAllIP(halfIndex)),floor(sortAllIP(halfIndex)));
        %                 sortAllIP=sortAllIP(~halfIndex);
        %                 sortAllIP=union(sortAllIP,Add);
        %
        %                 XoverCoverRange= abs(InsertPosition(:,1)-InsertPosition(:,2));
        % %                 XoverCoverRange(1) < XoverCoverRange(2)  %
        %                  if                   InsertPosition(1,2) < InsertPosition(2,2)
        % %                     C1=[ CellA(1:InsertPosition(1,2),:) ; Xover(1,4:6) ; Xover(1,1:3) ; CellA(InsertPosition(2,2)+1:end,:)] ;
        % %                     C2=[  Xover(2,1:3); CellA(InsertPosition(1,2)+1:InsertPosition(2,2),:) ;   Xover(2,4:6)];
        %
        %                     C1=[ CellA(1:sortAllIP(1),:) ; Xover(1,4:6) ; Xover(1,1:3) ; CellA(sortAllIP(4):end,:)] ;
        %
        %
        %                     C2=[  Xover(2,1:3); CellA( sortAllIP(2):sortAllIP(3),:) ;   Xover(2,4:6)];
        %                     UseFcnXTo1C1=1;
        %                  else
        % %                      setdiff(  Xover(2,1:3),CellA(InsertPosition(1),:) ,'rows')
        %
        %
        % %                     C1=[ CellA(1:floor(InsertPosition(1,2)),:) ; Xover(1,4:6) ; Xover(1,1:3) ; CellA(ceil(InsertPosition(1,1)):end,:)] ;
        % %                     C2=[  Xover(2,1:3); CellA( ceil(InsertPosition(2,1)):floor(InsertPosition(2,2)),:) ;   Xover(2,4:6)];
        % %
        % %                     %----test
        %                     C1=[ CellA(1:sortAllIP(1),:) ; Xover(2,4:6) ; Xover(2,1:3) ; CellA(sortAllIP(4):end,:)] ;
        %                     C2=[  Xover(1,1:3); CellA( sortAllIP(2):sortAllIP(3),:) ;   Xover(1,4:6)];
        % %
        % %
        %                     seFcnXTo1C2=2;
        %                  end
        %                CC1= unique(C1,'rows','stable');
        %                CC2= unique(C2,'rows','stable');
        %                 if mod(size(CC1,1),2)==1 || mod(size(CC2,1),2)==1
        %                     sdfsf=234;
        %                 end
        % %                 XXX{1}=C1
        % %                 obj.drawBCBdata2(InputCell{2},5)
        %
        %                InputCell(takeU)=[];
        %                InputCell{end+1}=CC1;
        %                InputCell{end+1}=CC2;
        %               end
        %             NOriGPairBCB=InputCell;
        %         end   % End of func AddXover2OneCycle
        
        
        
        function [NOriGPairBCB, UnUsedXover,takeU]=AddXover22cycles(obj,OriGPairBCB,Xover)
            UnUsedXover=[];
            obj.OverRoute=2;
            InputCell=OriGPairBCB;
            FourPoints=[Xover(1,1:3);Xover(2,4:6);Xover(2,1:3);Xover(1,4:6)];
            %             FourPoints=[Xover(:,1:3);Xover(:,4:6)];   %Dec30
            %             FourPoints=sortrows(FourPoints);
            whichcell=zeros(1,4);
            n=1;
            for iP=1:4   %search for which two cell contains these four points
                PP=FourPoints(iP,:);
                for jCell=1:length(OriGPairBCB)
                    Cycle=OriGPairBCB{jCell};
                    Mix=sortrows([Cycle;PP]);
                    [~,indexofPP]= ismember(PP,Mix,'rows');
                    if sum( ismember(PP, Cycle,'rows') )==1 %Xover locates on side
                        whichcell(n)= jCell;
                        n=n+1;
                    elseif mod(indexofPP,2)==0  && indexofPP~=1 && indexofPP~=size(Mix,1)
                        if size(unique(  Mix(indexofPP-1:indexofPP+1, 1:2),'rows') ,1)==1 %check on same edge
                            if (Mix(indexofPP,3)-Mix(indexofPP+1,3))*(Mix(indexofPP,3)-Mix(indexofPP-1,3))<=0  %z value between
                                whichcell(n)= jCell;
                                n=n+1;
                            end
                        end
                    end
                end
            end
            %----finish search
            takeU=setdiff(unique(whichcell),0) ;
            NumberOfCylceInvolved= length(  takeU  )    ;
            
            if NumberOfCylceInvolved==2
                
                % special care of compound joint,  subs into over connection
                AllFCColumn=[obj.ForcedConnectList(1:2:end,:),obj.ForcedConnectList(2:2:end,:)];
                XoverColumn=[Xover(1,:),Xover(2,:)];
                [~,XoverThisTime]=ismember(XoverColumn,AllFCColumn,'rows');
                %             if XoverThisTime==0
                %                 sdfd=234
                %             end
                %             VVCB=obj.CompJointInfo.BinRepofCJ(:,1);
                
                XoverIsSide=zeros(2,2);
                for xi=1:2
                    for yj=1:2
                        PointBCB= Xover(xi,3*yj-2:3*yj);
                        if PointBCB(3)==obj.containBundle{PointBCB(1)}.Zbase1(PointBCB(2))-obj.OverRoute
                            XoverIsSide( xi,  yj)=1;
                        elseif PointBCB(3)==obj.containBundle{PointBCB(1)}.Zbase2(PointBCB(2))+obj.OverRoute
                            XoverIsSide( xi,  yj)=1;
                        end
                    end
                end
                
                
                
                
                %
                %             if XoverThisTime~=0 &&  sum(sum(XoverIsSide))==4
                %                 if min(find(VVCB==VVCB(XoverThisTime)))~=XoverThisTime
                %                     WWW= AllFCColumn(find(VVCB==VVCB(XoverThisTime)),:);
                %                     CompoudJointNodes=unique([WWW(:,1:3);WWW(:,4:6);WWW(:,7:9);WWW(:,10:12)],'rows');
                %                     WWW2=[WWW(:,1:6);WWW(:,7:12)];
                %                       CellA2= OriGPairBCB{takeU(1)};
                %                       CellB2= OriGPairBCB{takeU(2)};
                %                       nA=1;
                %                     while nA<200
                %                        CellA2=circshift(CellA2,1) ;
                %                        nA=nA+1;
                %                        if ismember(CellA2(1,:),CompoudJointNodes,'rows') && CellA2(1,1)==CellA2(2,1) && CellA2(1,2)==CellA2(2,2)
                %                           break
                %                        end
                %                     end
                %                     nB=1;
                %                     while nB<200
                %                        CellB2=circshift(CellB2,1) ;
                %                        nB=nB+1;
                %                        if ismember(CellB2(1,:),CompoudJointNodes,'rows') && CellB2(1,1)==CellB2(2,1) && CellB2(1,2)==CellB2(2,2)
                %                           break
                %                        end
                %                     end
                %                     if nA>195 || nB>195
                %                         sdfsf=234;
                %                     end
                %
                %                          InputCell(takeU)=[];
                %                         InputCell{end+1}=[CellA2;CellB2];
                %                          NOriGPairBCB=InputCell;
                %                          return
                %                 end
                %             end
                %-------------
                
                
                
                
                Xover;
                CellA= OriGPairBCB{takeU(1)};
                CellB= OriGPairBCB{takeU(2)};
                %---insert 4points
                for iP2=1:2   %2point add in one time
                    if sum(ismember(FourPoints(2*iP2-1:2*iP2,:),CellA,'rows'))==2  %xover locates on side
                        [~,B]=ismember(FourPoints(2*iP2-1:2*iP2,:),CellA,'rows');
                        while length( unique([B;1  ;size(CellA,1)]  ) )==2
                            CellA= circshift(CellA,2);
                            [~,B]=ismember(FourPoints(2*iP2-1:2*iP2,:),CellA,'rows');
                        end
                        NcellA=CellA;
                    else
                        for icellA=1:size(CellA,1)-1
                            QQA= [  CellA(icellA,:) ; FourPoints(2*iP2-1:2*iP2,:) ; CellA(icellA+1,:)];
                            if  and(  sum(QQA(1,1)*ones(1,4) == QQA(:,1))==4 ,  sum(QQA(1,2)*ones(1,4) == QQA(:,2))==4 )
                                %                           if size(unique(QQA(:,1:2),'rows'),1)==1
                                Insettest=[ CellA(icellA,:) ; mean(FourPoints(2*iP2-1:2*iP2,:)) ; CellA(icellA+1,:)];
                                if size(union(  Insettest(:,1:2) ,[-1 -1 ],'rows') ,1)==2 %check on same edge
                                    if  (Insettest(2,3)- Insettest(1,3))*(Insettest(2,3)- Insettest(3,3))<=0
                                        FirstChoice=[ CellA(icellA,:) ; FourPoints(2*iP2-1:2*iP2,:) ; CellA(icellA+1,:)];
                                        SecondChoice=[ CellA(icellA,:) ; flip(FourPoints(2*iP2-1:2*iP2,:)) ; CellA(icellA+1,:)];
                                        if  sum( abs(diff(FirstChoice(:,3)))) <  sum( abs(diff(SecondChoice(:,3))))
                                            NcellA =[ CellA(1:icellA,:) ; FourPoints(2*iP2-1:2*iP2,:) ;  CellA(icellA+1:end,:)];
                                        else
                                            NcellA =[ CellA(1:icellA,:) ;flip( FourPoints(2*iP2-1:2*iP2,:)) ;  CellA(icellA+1:end,:)];
                                        end
                                    end
                                end
                            end
                            
                        end
                    end
                    
                    if sum(ismember(FourPoints(2*iP2-1:2*iP2,:),CellB,'rows'))==2  %xover locates on side
                        [~,B]=ismember(FourPoints(2*iP2-1:2*iP2,:),CellB,'rows');
                        while length( unique([B;1  ;size(CellB,1)]  ) )==2
                            CellB= circshift(CellB,2);
                            [~,B]=ismember(FourPoints(2*iP2-1:2*iP2,:),CellB,'rows');
                        end
                        NcellB=CellB;
                    else
                        for icellB=1:size(CellB,1)-1
                            QQB= [ CellB(icellB,:) ; FourPoints(2*iP2-1:2*iP2,:) ; CellB(icellB+1,:)];
                            %                             if   sum(CellB(icellB,:) ==FourPoints(2*iP2-1,:) )==2
                            if and(  sum(QQB(1,1)*ones(1,4) == QQB(:,1))==4 ,  sum(QQB(1,2)*ones(1,4) == QQB(:,2))==4 )
                                %                           if size(unique(QQB(:,1:2),'rows'),1)==1
                                Insettest=[ CellB(icellB,:) ; mean(FourPoints(2*iP2-1:2*iP2,:)) ; CellB(icellB+1,:)];
                                if size(union(  Insettest(:,1:2) ,[-1 -1 ],'rows') ,1)==2 %check on same edge
                                    if  (Insettest(2,3)- Insettest(1,3))*(Insettest(2,3)- Insettest(3,3))<=0
                                        FirstChoice=[ CellB(icellB,:) ; FourPoints(2*iP2-1:2*iP2,:) ; CellB(icellB+1,:)];
                                        SecondChoice=[ CellB(icellB,:) ; flip(FourPoints(2*iP2-1:2*iP2,:)) ; CellB(icellB+1,:)];
                                        if  sum( abs(diff(FirstChoice(:,3)))) <  sum( abs(diff(SecondChoice(:,3))))
                                            NcellB =[ CellB(1:icellB,:) ; FourPoints(2*iP2-1:2*iP2,:) ;  CellB(icellB+1:end,:)];
                                        else
                                            NcellB =[ CellB(1:icellB,:) ;flip( FourPoints(2*iP2-1:2*iP2,:)) ;  CellB(icellB+1:end,:)]  ;
                                        end
                                    end
                                end
                            end
                            
                        end
                    end
                end
                
                
                if  ~exist('NcellA') ||  ~exist('NcellB')
                    NOriGPairBCB=InputCell;
                    UnUsedXover=Xover;
                    return
                end
                Xover;  NcellA;      NcellB;
                XoverIndexInNcell=zeros(4,2);   % 1st Column:1->in A, 2-> in B ; 2nd Column: which index
                for iX=1:2
                    Test=Xover(iX,1:3);
                    [TFA,IIndexA]=ismember(Test,NcellA,'rows');
                    [TFB,IIndexB]=ismember(Test,NcellB,'rows');
                    if TFA==1;  XoverIndexInNcell(2*iX-1,1:2)=[1,IIndexA];
                    elseif TFB==1  ;  XoverIndexInNcell(2*iX-1,1:2)=[2,IIndexB]; end
                    
                    Test=Xover(iX,4:6);
                    [TFA,IIndexA]=ismember(Test,NcellA,'rows');
                    [TFB,IIndexB]=ismember(Test,NcellB,'rows');
                    if TFA==1;  XoverIndexInNcell(2*iX,1:2)=[1,IIndexA];
                    elseif TFB==1  ;  XoverIndexInNcell(2*iX,1:2)=[2,IIndexB]; end
                end
                %                 XoverIndexInNcell=[1 3; 2 6; 2 7; 1 2];
                %FFF
                %=---------
                OO=XoverIndexInNcell;
                ssXIIN=sortrows(XoverIndexInNcell);
                if ssXIIN(2,2)~=ssXIIN(1,2)+1
                    ssXIIN2=ssXIIN;
                    ssXIIN2(2,2)=ssXIIN2(1,2)+1;
                    [~,ind]=ismember(XoverIndexInNcell,ssXIIN,'rows');
                    XoverIndexInNcell=ssXIIN2(ind,:);
                elseif ssXIIN(4,2)~=ssXIIN(3,2)+1
                    ssXIIN2=ssXIIN;
                    if mod(ssXIIN2(3,2),2)==0
                        ssXIIN2(4,2)=ssXIIN2(3,2)+1;
                    else
                        ssXIIN2(4,2)=ssXIIN2(3,2)-1;
                    end
                    [~,ind]=ismember(XoverIndexInNcell,ssXIIN,'rows');
                    XoverIndexInNcell=ssXIIN2(ind,:);
                else
                    
                end
                %----------^^^^^^new  Feb 2
                
                TodoWithA= XoverIndexInNcell(XoverIndexInNcell(:,1)==1,2);
                NcellAV2=[NcellA(max(TodoWithA):end,:) ; NcellA(1:min(TodoWithA),:)];
                
                TodoWithB= XoverIndexInNcell(XoverIndexInNcell(:,1)==2,2);
                NcellBV2=[NcellB(max(TodoWithB):end,:) ; NcellB(1:min(TodoWithB),:)];
                InputCell(takeU)=[];
                InputCell{end+1}=[NcellAV2;NcellBV2];
            else
                UnUsedXover=Xover;
            end
            
            NOriGPairBCB=InputCell;
        end   % end of function AddXover22cycles
        
        
        function NCList=CListsortCListNotFollow53(obj,ConnectList)
            KK=zeros(size(ConnectList));
            for si=1:size(ConnectList)
                TwoPoints= ConnectList(si,:);
                Tunr90=[TwoPoints(1:3);TwoPoints(4:6)];
                N90=sortrows( Tunr90);
                
                KK(si,1:6)=[N90(1,:),N90(2,:)];
            end
            NCList=sortrows(KK);
            
        end
        
        function drawBCBdata(obj,varargin)
            %             pause(1)
            %                         figure(123);clf;hold on ;
            %                         aH=gca;
            %             NNN=nargin
            
            cycleList=varargin{1} ;
            if nargin==2   % default target axes
                aH= findobj(fH,'Tag','MechScaffold2D');
                
            elseif nargin==3  % if assign axe handle
                aH= varargin{2} ;
            end
            
            %             aH= findobj(fH,'Tag','MechScaffold2D');
            axes(aH);  cla; hold on;
            %             aH=gca;
            %             set(aH,'ytick',[]);
            aH.YColor=[1,1,1];
            %             set(gcf,'Color',[1 1 1 ]);
            n=0;
            minmin=9000;
            for k=1:length(cycleList)
                n=n+1;
                Cycle= cycleList{k};
                xydata=zeros( size(Cycle,1),2);
                for k2=1: size(Cycle,1)
                    xydata( k2,1)=Cycle(k2,3)  ;    %uint base ,draw in Xaxis
                    [~,GCylIndex]=ismember(Cycle(k2,1:2),obj.GlobalCylinderIndex,'rows');
                    xydata( k2,2)=  GCylIndex + 5*Cycle(k2,1)  ;
                    %                                    xydata( k2,2)=  GCylIndex ;
                    
                end
                xydata(end+1,:) =xydata(1,:) ;
                
                
                ditt=0*rand();
                %                plot([xydata(:,1);xydata(1,1)]+ditt*ones(length(xydata(:,1))+1,1),-[xydata(:,2);xydata(1,2)]);
                plot([xydata(:,1)]+ditt*ones(length(xydata(:,1)),1),-[xydata(:,2)]);
                %                 text(xydata(1,1)+ditt,-xydata(1,2), num2str(k));
                %                text(xydata(2,1)+ditt,-xydata(2,2), num2str(k+0.1));
                minmin=min([minmin; xydata(:,1)]);
            end
            title(strcat(' Number of Cycles = ', num2str(n)));
            
            %             for show=1:size(obj.GlobalCylinderIndex,1)
            %                 text(minmin-10, -show -5*obj.GlobalCylinderIndex(show,1) , strcat( num2str(obj.GlobalCylinderIndex(show,1)), '-' ,  num2str(obj.GlobalCylinderIndex(show,2))));
            %                 text(minmin-20,-show-5*obj.GlobalCylinderIndex(show,1) ,strcat( num2str(mod(obj.RelateTable(show,4),2)==0) ) ,'Color','red');
            %             end
            axis off
        end
        
        
        
        function drawBCBdata2(obj,cycleList,num)
            %             pause(1)
            figure(num);clf;hold on ; aH=gca;
            %             set(aH,'ytick',[]);
            aH.YColor=[1,1,1];
            set(gcf,'Color',[1 1 1 ]);
            n=0;
            minmin=9000;
            for k=1:length(cycleList)
                n=n+1;
                Cycle= cycleList{k};
                xydata=zeros( size(Cycle,1),2);
                for k2=1: size(Cycle,1)
                    xydata( k2,1)=Cycle(k2,3);    %uint base ,draw in Xaxis
                    [~,GCylIndex]=ismember(Cycle(k2,1:2),obj.GlobalCylinderIndex,'rows');
                    xydata( k2,2)=  GCylIndex;
                end
                ditt=0*rand();
                plot([xydata(:,1);xydata(1,1)]+ditt*ones(length(xydata(:,1))+1,1),-[xydata(:,2);xydata(1,2)]);
                text(xydata(1,1)+ditt,-xydata(1,2), num2str(k));
                text(xydata(2,1)+ditt,-xydata(2,2), num2str(k+0.1));
                minmin=min([minmin; xydata(:,1)]);
            end
            title(strcat(' Number of Cycle = ', num2str(n)));
            
            for show=1:size(obj.GlobalCylinderIndex,1)
                text(minmin-10, -show, strcat( num2str(obj.GlobalCylinderIndex(show,1)), '-' ,  num2str(obj.GlobalCylinderIndex(show,2))));
                text(minmin-20,-show ,strcat( num2str(mod(obj.RelateTable(show,4),2)==0) ) ,'Color','red');
            end
        end
        function obj=UpdataConnect(obj,targetAxeHandle)
            axes(targetAxeHandle)
            if ~isempty(obj.UConPlotHandle)
                delete(obj.UConPlotHandle{:});
            end
            PlotList= obj.ForcedConnectList;
            plotXYZ=zeros(size(PlotList,1),6);
            for linei=1:    size(PlotList,1)
                K=PlotList(linei,:);
                alpha=K(3)-obj.containBundle{K(1)}.Zbase1(K(2));
                beta=obj.containBundle{K(1)}.Zbase2(K(2))-K(3);
                twoEnds=obj.containBundle{K(1)}.CylinderXYZGlobal(K(2),:);
                P1=( beta*twoEnds(1:3)+  alpha*twoEnds(4:6))/(alpha+beta);
                
                alpha2=K(6)-obj.containBundle{K(4)}.Zbase1(K(5));
                beta2=obj.containBundle{K(4)}.Zbase2(K(5))-K(6);
                twoEnds2=obj.containBundle{K(4)}.CylinderXYZGlobal(K(5),:);
                P2=( beta2*twoEnds2(1:3)+  alpha2*twoEnds2(4:6))/(alpha2+beta2);
                
                plotXYZ(linei,1:6)=[P1 , P2];
            end
            
            obj.UConPlotHandle{1}=plot3([plotXYZ(:,1), plotXYZ(:,4)]' , [plotXYZ(:,2),plotXYZ(:,5)]' ,[plotXYZ(:,3),plotXYZ(:,6)]' );
        end
        
        
        
        function  [TDoable,CellPairList,CellUnpair]=checkpairableH(obj,preference)
            
            % prefercase=1 ; % XY ramdon
            % prefercase=2 ; % Y
            % prefercase=3 ; % X
            preferenceSS=preference ;
            NumberOfBundle=length(obj.containBundle);
            CellPairList=cell(1,NumberOfBundle);
            CellUnpair=cell(1,NumberOfBundle);
            TDoable=true;
            
            for Nbundle=1:NumberOfBundle
                adjM=obj.containBundle{Nbundle}.NarrowAdj;
                [U,V]=find(adjM==1);
                %                 if ismember(Nbundle, [5,4] )  % hard code
                %                     preference=2  ;
                %                    fprintf('hard code for pairing \n')
                %                 else
                %                     preference=preferenceSS;
                %                 end
                %
                % only square lattice has preferred pairing options.
                %                 Sp_hard=[9 10 11 12] ;  %hard code
                if preference==3  && strcmp(obj.containBundle{Nbundle}.type,'SQ')     %|| ismember(Nbundle,[3,4])
                    for deletList=1:length(U)
                        Uxy=obj.containBundle{Nbundle}.CylInplanePosition(U(deletList),:);
                        Vxy=obj.containBundle{Nbundle}.CylInplanePosition(V(deletList),:);
                        if Uxy(1)==Vxy(1)
                            adjM(  U(deletList),V(deletList))=0;
                        end
                    end
                elseif preference==2  && strcmp(obj.containBundle{Nbundle}.type,'SQ')
                    for deletList=1:length(U)
                        Uxy=obj.containBundle{Nbundle}.CylInplanePosition(U(deletList),:);
                        Vxy=obj.containBundle{Nbundle}.CylInplanePosition(V(deletList),:);
                        if Uxy(2)==Vxy(2)  % && ~ismember( U(deletList) ,Sp_hard)
                            adjM(  U(deletList),V(deletList))=0;
                        end
                    end
                end
                ResetAdjM = adjM ;
                
                %-------hard code for pairing
                %                 if Nbundle==3
                %                     adjM(2,3)=0;  adjM(3,2)=0;
                %                     adjM(4,5)=0;  adjM(5,4)=0;
                %                     adjM(18,19)=0;  adjM(19,18)=0;
                %                     adjM(16,17)=0;  adjM(17,16)=0;
                %                end
                
                
                for mulifooLoop = 1 :10
                    adjM= ResetAdjM ;
                    PairList=zeros(ceil((size(obj.containBundle{Nbundle}.CylAdjMat,1))/2) ,2);
                    nP=0; N=0; Doable=false;      Count=0;
                    while N<100
                        Anything=0;
                        N=N+1;
                        QQ=sum(adjM); rr=nnz(QQ);
                        if rr==0; Doable=false;
                            break
                        end
                        WW=find(QQ==1);
                        cc=rand;
                        if ~isempty(WW) && cc<0.9
                            tt=WW(randi([1 length(WW)]));
                            potent=find(adjM(tt(1),:)==1); %  length(potent)
                            wa=potent(randi([1 length(potent)]));
                            if abs(obj.containBundle{Nbundle}.Zbase1(wa)-obj.containBundle{Nbundle}.Zbase1(tt))<=obj.containBundle{Nbundle}.Tol && abs(obj.containBundle{Nbundle}.Zbase2(wa)-obj.containBundle{Nbundle}.Zbase2(tt))<=obj.containBundle{Nbundle}.Tol
                                nP=nP+1;
                                PairList(nP,1)=tt;
                                PairList(nP,2)=wa;
                                adjM(tt(1),:)=0;adjM(:,tt)=0;
                                adjM(wa,:)=0;adjM(:,wa)=0;
                                Anything=1;
                            end
                        else
                            yy=find(sum(adjM)==2) ;
                            if ~isempty(yy)
                                yy=yy(randi([1 length(yy)]));
                                xx=find(adjM(yy,:)==1);
                                xx=xx(randi([1 length(xx)])) ;
                                if abs(obj.containBundle{Nbundle}.Zbase1(yy)-obj.containBundle{Nbundle}.Zbase1(xx))<=obj.containBundle{Nbundle}.Tol && abs(obj.containBundle{Nbundle}.Zbase2(yy)-obj.containBundle{Nbundle}.Zbase2(xx))<=obj.containBundle{Nbundle}.Tol
                                    nP=nP+1;
                                    PairList(nP,1)=yy;
                                    PairList(nP,2)=xx;
                                    adjM(yy,:)=0;adjM(:,yy)=0;
                                    adjM(xx,:)=0;adjM(:,xx)=0;
                                    Anything=1;
                                end
                            end
                        end
                        
                        if Anything==0
                            Count=Count+1;
                        end
                        if nnz(PairList)==size(PairList,1)*size(PairList,2)
                            Doable=true;    break; end
                    end   %end of while
                    
                    if   Doable==true;  break ;end  % break for loop
                    
                end  %end of for
                
                
                used=PairList(:);
                Unpair=setdiff(setdiff(union(obj.containBundle{Nbundle}.AGroup,obj.containBundle{Nbundle}.BGroup),used),path);
                PairList(sum(PairList,2)==0,: )=[];
                CellPairList{Nbundle}=PairList;
                CellUnpair{Nbundle}=Unpair;
                TDoable=and(TDoable,Doable);
            end
            
        end  % end of checkpairableH
        
        function obj=FindStap(obj,type)   % Get properties: StapList , stapBP ,stapBPinCell
            InitialStappEnds=zeros( 2*max(obj.RelateTable(:,3)),2); %hard
            
            if strcmp(obj.ScafOption.prevStack,'scafloop')   % remember
                for i=1:size(InitialStappEnds,1)/2
                    Bundle=obj.containBundle{obj.RelateTable(i,  1)};
                    InitialStappEnds(2*i-1:2*i,1)= obj.RelateTable(i,  3);   % C3 notation
                    Cyl=obj.RelateTable(i,  2);
                    InitialStappEnds(2*i-1,2)= Bundle.Zbase1(Cyl);
                    
                    InitialStappEnds(2*i,2)=Bundle.Zbase2(Cyl);
                    %                     InitialStappEnds(2*i-1,2)= Bundle.Zbase1(Cyl);
                    %                     InitialStappEnds(2*i,2)=Bundle.Zbase2(Cyl);
                    
                end
            elseif strcmp(obj.ScafOption.prevStack,'polyT')
                ForceConect =zeros(length(obj.ssOption.ForceConnUpdate) , 6) ;
                for k=1: size(ForceConect,1)
                    ForceConect(k,:) =obj.ssOption.ForceConnUpdate{k}.BCB0 ;
                end
                for i=1:size(InitialStappEnds,1)/2
                    InitialStappEnds(2*i-1:2*i,1)= obj.RelateTable(i,  3);   % C3 notation
                    Bundle=obj.containBundle{obj.RelateTable(i,  1)};
                    Cyl=obj.RelateTable(i,  2);
                    BCBZ1 =  [obj.RelateTable(i,  1) , Cyl , Bundle.Zbase1(Cyl) ] ;
                    BCBZ2 =  [obj.RelateTable(i,  1) , Cyl , Bundle.Zbase2(Cyl) ] ;
                    
                    if ismember( BCBZ1 ,ForceConect(:,1:3) ,'rows') || ismember( BCBZ1 ,ForceConect(:,4:6) ,'rows')
                        InitialStappEnds(2*i-1,2)= Bundle.Zbase1(Cyl) ;
                    else
                        InitialStappEnds(2*i-1,2)= Bundle.Zbase1(Cyl)-obj.ssOption.OverRoute;
                    end
                    
                    if ismember( BCBZ2 ,ForceConect(:,1:3) ,'rows') || ismember( BCBZ2 ,ForceConect(:,4:6) ,'rows')
                        InitialStappEnds(2*i,2)=Bundle.Zbase2(Cyl);
                    else
                        InitialStappEnds(2*i,2)=Bundle.Zbase2(Cyl)+obj.ssOption.OverRoute;
                    end
                end
            end
            
            
            for newIniStaEndi =1:size(InitialStappEnds,1)
                BCBrep=[obj.RelateTable( find( obj.RelateTable(:,3)==InitialStappEnds(newIniStaEndi,1)) ,1:2),InitialStappEnds(newIniStaEndi,2)];
                if BCBrep(3)==obj.containBundle{BCBrep(1)}.Zbase1(BCBrep(2))
                    InitialStappEnds(newIniStaEndi,2)=InitialStappEnds(newIniStaEndi,2)+obj.ssOption.BundleShiftTwoSide(BCBrep(1),1) ;
                elseif BCBrep(3)==obj.containBundle{BCBrep(1)}.Zbase2(BCBrep(2))
                    InitialStappEnds(newIniStaEndi,2)=InitialStappEnds(newIniStaEndi,2)-obj.ssOption.BundleShiftTwoSide(BCBrep(1),2) ;
                end
            end
            
            %------ask if want to use overhangs
            table_OH=findobj(gcf,'Tag','OHTable') ;
            if ~isempty(table_OH.Data)
                answer = questdlg('Would you like to use overhangs ?','Overhang option', ...
                    'Yes','No' ,'Yes');
                if strcmp(answer,'Yes')
                    
                    Info=  table_OH.Data; QQ= cell2mat(Info(:,6:7) ) ; 
                    %--------correct impratical inputs
                    for s_inf =   1:size(Info,1)
                        if  contains( Info{s_inf,8} , 'Double crossover')
                            fprintf('For double crossover, the input of directionality in Row %i will be ignored. \n' , s_inf);
                            table_OH.Data{s_inf,2} = '3''' ;
                            table_OH.Data{s_inf,5} = '3''' ;
                        elseif contains( Info{s_inf,8} , 'Double overhang1')
                            fprintf('For Double overhang1, the input of directionality in Row %i will be ignored. \n' , s_inf);
                            table_OH.Data{s_inf,2} = '3''' ;
                            table_OH.Data{s_inf,5} = '5''' ;
                        elseif contains( Info{s_inf,8} , 'Double overhang2')
                            fprintf('For Double overhang2, the input of directionality in Row %i will be ignored. \n' , s_inf);
                            table_OH.Data{s_inf,2} = '5''' ;
                            table_OH.Data{s_inf,5} = '3''' ;         
                        end
                    end
                    %---------
                    
                    QQ=QQ(:) ; QQ2 = repelem(QQ, 2) ;% for 0-nt overhang
                    IndEff= cellfun(@isequal,Info(:,3) , num2cell(ones(size(Info(:,3) ) ) ) ) ;
                    InfoHide= table_OH.UserData ;
                    Nick= InfoHide(IndEff,[1 4]);   %Nick= InfoHide(:,[1 4]);
                    Nick=Nick(:);
                    ConvertC3= zeros(length(Nick)*2 , 2) ;
                    for k=1:length(Nick)
                        Arr= str2num(Nick{k} ) ; %#ok<ST2NM>
                        [~,Ind]= ismember(  Arr(1:2), obj.RelateTable(:,1:2),'rows') ;
                        C3_ind=  obj.RelateTable(Ind,3) ;
                        ConvertC3(2*k-1,:) =[C3_ind, Arr(3)] ;
                        ConvertC3(2*k,:) =[C3_ind, Arr(4)] ;
                    end
%                     ConvertC3=ConvertC3(QQ2~=0 ,:) ;
                    
                    InitialStappEndsOri=InitialStappEnds;
                    InitialStappEnds=[InitialStappEnds; ConvertC3];  InitialStappEnds=sortrows(InitialStappEnds);
                    obj.UserWantOH='Yes' ;
                    %------hard coding scripts 
%                     stepAssign=2 ;  Ex_DecoratePlatesWithOverhangs ;
                    
                else
                    obj.UserWantOH='No' ;
                end
            else
                obj.UserWantOH='No' ;
                fprintf(' No overhang is used \n')
            end
            
            n=size(InitialStappEnds,1)/2;
            StapleList=cell(n,1);
            for u=1:size(InitialStappEnds,1)/2
                StapleList{u}=InitialStappEnds(2*u-1:2*u,:) ;
            end
            
            StapleCell=StapleList;
            % convert to C5
            for k=1:length(StapleCell)
                TT=StapleCell{k};
                for k2=1:size(   TT,1)
                    TT(k2,1)= obj.RelateTable( TT(k2,1),5);
                end
                StapleCell{k}=TT;
            end
            [ StapleCell ] = CalibStapDir( StapleCell,obj.RelateVec);
            
            %            if strcmp(obj.StapOption.halfXover,'yes')  && strcmp(obj.StapOption.type,'straightuncut')
            %             scriptcase=1;%----------------halfCrossover
            %             AddHalfXoverV3LargeCrossSec;   %input: StapleCell
            %            end
            
            
            obj.StapList=StapleCell;  %cylinder index use C3;
            
            %          GetHyperB.StapList =StapleCell ;
            
            %             obj.StapList=StapleList;  %cylinder index use C3;
            %                 return
            
            if type==1  %no xover between cyl (initial stap)
                return;
            else
                AdjM=obj.CylinderC5AdjM;    %use CylAdjM in hyperbundle
                if strcmp(obj.ScafOption.prevStack ,'polyT')
%                     clearance=15;
%                     clearanceNew = 15*ones(length(obj.containBundle) ,2) ;
                       FFh=gcf;
                      clearanceNew = QueryStapleClr(obj,FFh)   
                else
%                     clearance=3;   % 8 change 6/19 : default 8

%                     opts.Interpreter = 'tex';
% 
%                     prompt = {' \fontsize{12} Enter the clearance that staple Xovers are ignored from ends:  (Min=2)'};
%                     dlgtitle = 'Input';
%                     dims = [1 60];
%                     definput = {num2str(8)};
%                     answer = inputdlg(prompt,dlgtitle,dims,definput,opts) ;
%                     QQ=  round(str2double(answer{1} ) ) ;
%                    
%                    if QQ>1
%                      clearance=QQ ;  
%                    else
%                     clearance=8 ;   
%                    end
%                    clearance
                   
                   FFh=gcf;
                   clearanceNew = QueryStapleClr(obj,FFh) 
                   clearanceNew(:,2)=clearanceNew(:,2)-1 ;
                   clearanceNew(:,1)=clearanceNew(:,1)-1 ;
                end
                
                [U,V]=find(AdjM~=0);   %means cylinder in C5
                U2=U;  V2=V; %means cylinder in C5
                Edgeof5Index=  unique([U2 V2],'rows');  %C5
                k=1;
                EdgeList=cell(size(Edgeof5Index,1)/2,2);
                EdgeType=zeros(size(Edgeof5Index,1)/2,1);   % 1->SQ ,  2->HC
                for i=1:size(Edgeof5Index,1)
                    Ind=find(Edgeof5Index(i,1)==obj.RelateTable(:,5));Ind=Ind(1);
                    Bundle=obj.containBundle{obj.RelateTable(Ind,1)};
                    if ~xor(ismember( obj.RelateTable(Ind,2)   ,Bundle.AGroup) ,Bundle.AGroupGoUp==1)
                        EdgeList{k,1}=[ Edgeof5Index(i,1) , Edgeof5Index(i,2)];   %in 5rd index;
                        Bundlei=obj.RelateTable(Ind,1);
                        %                     obj.containBundle{Bundlei}.
                        if  strcmp( obj.containBundle{Bundlei}.Lattice, 'Square')
                            EdgeType(k)=1;
                        elseif strcmp( obj.containBundle{Bundlei}.Lattice, 'HoneyComb')
                            EdgeType(k)=2;
                        end
                        k=k+1;
                    end
                    
                end
                InitialStappEndsC5=InitialStappEnds;
                for iInt=1:size(InitialStappEnds,1)
                    InitialStappEndsC5(iInt,1)= obj.RelateTable(InitialStappEndsC5(iInt,1),5);
                end
                
                DCellList=obj.ScafdigitSQ;
                G2=[15 16]  ; %West
                G4=[31 32]   ;%East
                G3=[7 8]   ;%South
                G1=[23 24]   ; %North
                if strcmp(obj.UserWantOH,'Yes') % 05072019 when using OH,
                    for iInt=1:size(InitialStappEndsOri,1)
                        InitialStappEndsOri(iInt,1)= obj.RelateTable(InitialStappEndsOri(iInt,1),5);
                    end
                end
                TwoDExpre=zeros(20000,2); nTwo=1;
                for edge=1:size(EdgeList,1)                    
                    PosCyl=EdgeList{edge,1}(1);    %refers global index, C5
                    NegCyl=EdgeList{edge,1}(2);
                    BundleInd = unique( union(obj.RelateTable(  obj.RelateTable(:,5)==PosCyl ,1),obj.RelateTable(  obj.RelateTable(:,5)==NegCyl ,1) )) ;
                    if length(BundleInd)~=1
                        error('Something wrong!!')
                    end
                    
                    if strcmp(obj.UserWantOH,'Yes') % 05072019 when using OH,
                        EndsOfPosCyl=InitialStappEndsOri(find(InitialStappEndsOri(:,1)==PosCyl),2);
                        EndsOfNegCyl=InitialStappEndsOri(find(InitialStappEndsOri(:,1)==NegCyl),2);
                    else
                        EndsOfPosCyl=InitialStappEndsC5(find(InitialStappEndsC5(:,1)==PosCyl),2);
                        EndsOfNegCyl=InitialStappEndsC5(find(InitialStappEndsC5(:,1)==NegCyl),2);
                    end
                    minP=min(EndsOfPosCyl);  maxP=max(EndsOfPosCyl);
                    minN=min(EndsOfNegCyl);  maxN=max(EndsOfNegCyl);
                    
%                     clearanceNew; Allow locally assign staple clearance
%                     for controlling sharp and round corners. 11142020,
%                     Curved-linear

                    minAll=max(minP,minN)+ clearanceNew(BundleInd,1);
                    maxAll=min(maxP,maxN)- clearanceNew(BundleInd,2) ;
%                     minAll=max(minP,minN)+ clearance;
%                     maxAll=min(maxP,maxN)- clearance ;

                    EndsOfPosCyl2=EndsOfPosCyl;
                    for kk2=1:length(EndsOfPosCyl2)
                        if mod(kk2,2)==1
                            EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)+ clearanceNew(BundleInd,1);
%                             EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)+ clearance;
                        else
                            EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)- clearanceNew(BundleInd,2);
%                             EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)- clearance;
                        end
                    end
                    
                    EndsOfNegCyl2=EndsOfNegCyl;
                    for kk2=1:length(EndsOfNegCyl2)
                        if mod(kk2,2)==1
                            EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)+clearanceNew(BundleInd,1);
%                             EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)+clearance;
                        else
                            EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)-clearanceNew(BundleInd,2);
%                             EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)-clearance;
                        end
                    end
                    
                    Pind=find(obj.RelateTable(:,5)==PosCyl);Pind=Pind(1);
                    Nind=find(obj.RelateTable(:,5)==NegCyl);Nind=Nind(1);
                    
                    if EdgeType(edge)==1    % this edge is square lattice
                        NWSE=FindNWSE(Pind,Nind,obj.RelateTable(:,6:7));    % for SQ
                        period=32;
                        switch NWSE
                            case 1
                                StartOfFor=G1(1);
                            case 2
                                StartOfFor=G2(1);
                            case 3
                                StartOfFor=G3(1);
                            case 4
                                StartOfFor=G4(1);
                        end
                        for zz=StartOfFor:period:maxAll
                            zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                            zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                            if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                                EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                                TwoDExpre(nTwo:nTwo+3 ,:) =[PosCyl zz ;NegCyl zz;PosCyl zz+1;NegCyl zz+1] ;nTwo=nTwo+4 ;
                                
                                %                                 TwoDExpre(end+1,:)=[PosCyl zz];    TwoDExpre(end+1,:)=[NegCyl zz];
                                %                                 TwoDExpre(end+1,:)=[PosCyl zz+1];  TwoDExpre(end+1,:)=[NegCyl zz+1];
                            end
                        end
                    else
                        Bundle= obj.containBundle{ obj.RelateTable(Pind,1)};
                        
                        Pind2=obj.RelateTable(Pind,2);
                        Nind2=obj.RelateTable(Nind,2);
                        %                         obj.RelateTable(Pind,1)
                        %                         [Pind,Nind]
                        [ Dir ] = FindHCDir(Pind2,Nind2,Bundle.CylInplanePosition) ;
                        %                         [ Dir ] = FindHCDir(Pind,Nind,Bundle.CylInplanePosition) ;
                        %                         [ Dir ] = FindHCDir(Pind,Nind,obj.RelateTable(:,6:7)) ;
                        period=21;
                        switch Dir
                            case 1   %+90
                                StartOfFor2=13;
                            case 2  %-150
                                StartOfFor2=20;
                            case 3  %-30
                                StartOfFor2=6;
                        end
                        for zz=StartOfFor2:period:maxAll
                            zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                            zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                            if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                                EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                                TwoDExpre(nTwo:nTwo+3 ,:) =[PosCyl zz ;NegCyl zz;PosCyl zz+1;NegCyl zz+1] ;nTwo=nTwo+4 ;
                                %                                 TwoDExpre(end+1,:)=[PosCyl zz];    TwoDExpre(end+1,:)=[NegCyl zz];
                                %                                 TwoDExpre(end+1,:)=[PosCyl zz+1];  TwoDExpre(end+1,:)=[NegCyl zz+1];
                            end
                        end
                    end
                end
                TwoDExpre=TwoDExpre(1:nTwo-1,:) ;
                
                obj.stapBP=TwoDExpre;
                obj.OstapBP=TwoDExpre;
                obj.stapBPinCell=EdgeList;
            end
        end %end of FindStap
        
        function getScafXoverAsStap(obj)
            ConnectList=obj.ForcedConnectList;
            ConnectList2 =[ConnectList(:,1:3) ;ConnectList(:,4:6)] ;
            InitialStappEnds=zeros(2*max(obj.RelateTable(:,3)),2);
            for i=1:size(InitialStappEnds,1)/2
                InitialStappEnds(2*i-1:2*i,1)= obj.RelateTable(i,  3);   % C3 notation
                Bundle=obj.containBundle{obj.RelateTable(i,  1)};
                Cyl=obj.RelateTable(i,  2);
                InitialStappEnds(2*i-1,2)= Bundle.Zbase1(Cyl);
                InitialStappEnds(2*i,2)=Bundle.Zbase2(Cyl);
            end
            
            AdjM=obj.CylinderC5AdjM;    %use CylAdjM in hyperbundle
            if strcmp(obj.ScafOption.prevStack ,'polyT')
                clearance=15;  clearance=8;
            else
                clearance=8;
            end
            [U,V]=find(AdjM~=0);   %means cylinder in C5
            U2=U;  %means cylinder in C5
            V2=V;
            Edgeof5Index=  unique([U2 V2],'rows');  %C5
            k=1;
            EdgeList=cell(size(Edgeof5Index,1)/2,3);
            EdgeType=zeros(size(Edgeof5Index,1)/2,1);   % 1->SQ ,  2->HC
            for i=1:size(Edgeof5Index,1)
                Ind=find(Edgeof5Index(i,1)==obj.RelateTable(:,5));Ind=Ind(1);  % Ind
                Ind2=find(Edgeof5Index(i,2)==obj.RelateTable(:,5));Ind2=Ind2(1);
                Bundle=obj.containBundle{obj.RelateTable(Ind,1)};
                if ~xor(ismember( obj.RelateTable(Ind,2)   ,Bundle.AGroup) ,Bundle.AGroupGoUp==1)
                    EdgeList{k,1}=[ Edgeof5Index(i,1) , Edgeof5Index(i,2)];   %in 5rd index;
                    EdgeList{k,3}=[ obj.RelateTable(Ind,1:2)  ; obj.RelateTable(Ind2,1:2)    ];   %in 5rd index;
                    
                    Bundlei=obj.RelateTable(Ind,1);
                    %                     obj.containBundle{Bundlei}.
                    if  strcmp( obj.containBundle{Bundlei}.Lattice, 'Square')
                        EdgeType(k)=1;
                    elseif strcmp( obj.containBundle{Bundlei}.Lattice, 'HoneyComb')
                        EdgeType(k)=2;
                    end
                    k=k+1;
                end
                
            end
            %--------------
            InitialStappEndsC5=InitialStappEnds;
            for iInt=1:size(InitialStappEnds,1)
                InitialStappEndsC5(iInt,1)= obj.RelateTable(InitialStappEndsC5(iInt,1),5);
            end
            %-----------------
            DCellList=obj.ScafdigitSQ;
            G2=[15 16]+16  ; %West
            G4=[31 32]+16   ;%East
            G3=[7 8]+16   ;%South
            G1=[23 24]+16   ; %North
            TwoDExpre=[]; Tol = 8 ;
            for edge=1:size(EdgeList,1)
                PosCyl=EdgeList{edge,1}(1);    %refers global index, C5
                NegCyl=EdgeList{edge,1}(2);
                
                EndsOfPosCyl=InitialStappEndsC5(find(InitialStappEndsC5(:,1)==PosCyl),2);
                EndsOfNegCyl=InitialStappEndsC5(find(InitialStappEndsC5(:,1)==NegCyl),2);
                minP=min(EndsOfPosCyl);
                maxP=max(EndsOfPosCyl);
                minN=min(EndsOfNegCyl);
                maxN=max(EndsOfNegCyl);
                
                minAll=max(minP,minN)+clearance;
                maxAll=min(maxP,maxN)-clearance ;
                EndsOfPosCyl2=EndsOfPosCyl;
                for kk2=1:length(EndsOfPosCyl2)
                    if mod(kk2,2)==1
                        EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)+clearance;
                    else
                        EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)-clearance;
                    end
                end
                
                EndsOfNegCyl2=EndsOfNegCyl;
                for kk2=1:length(EndsOfNegCyl2)
                    if mod(kk2,2)==1
                        EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)+clearance;
                    else
                        EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)-clearance;
                    end
                end
                
                Pind=find(obj.RelateTable(:,5)==PosCyl);Pind=Pind(1);
                Nind=find(obj.RelateTable(:,5)==NegCyl);Nind=Nind(1);
                
                if EdgeType(edge)==1    % this edge is square lattice
                    NWSE=FindNWSE(Pind,Nind,obj.RelateTable(:,6:7));    % for SQ
                    period=32;
                    switch NWSE
                        case 1
                            StartOfFor=G1(1);
                        case 2
                            StartOfFor=G2(1);
                        case 3
                            StartOfFor=G3(1);
                        case 4
                            StartOfFor=G4(1);
                    end
                    for zz=StartOfFor:period:maxAll
                        zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                        zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                        if   mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            %                         if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            
                            EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz];    TwoDExpre(end+1,:)=[NegCyl zz];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz+1];  TwoDExpre(end+1,:)=[NegCyl zz+1];
                            [ia,ib]=ismember( EdgeList{1,3}(1,:) ,obj.RelateTable(:,1:2) , 'rows') ;
                            if mod(obj.RelateTable(ib,4) ,2)==0
                                if sum(and(ismember(ConnectList2(:,1:2),EdgeList{edge,3}(2,:),'rows') , abs(zz-ConnectList2(:,3))<Tol ) )==0 % min distance from foce connection
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz , EdgeList{edge,3}(1,:) zz  ];
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz+1 , EdgeList{edge,3}(2,:) zz+1 ];
                                    %                                else
                                    %                                    sdfsf=3
                                end
                            else
                                if sum(and(ismember(ConnectList2(:,1:2),EdgeList{edge,3}(1,:),'rows') , abs(zz-ConnectList2(:,3))<Tol) )==0  % min distance from foce connection
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz , EdgeList{edge,3}(2,:) zz  ];
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz+1 , EdgeList{edge,3}(1,:) zz+1 ];
                                end
                            end
                            
                        end
                    end
                else
                    Bundle= obj.containBundle{ obj.RelateTable(Pind,1)};
                    Pind2=obj.RelateTable(Pind,2);
                    Nind2=obj.RelateTable(Nind,2);
                    [ Dir ] = FindHCDir(Pind2,Nind2,Bundle.CylInplanePosition) ;
                    period=21;
                    switch Dir
                        case 1   %+90
                            StartOfFor2=round(13+5);
                        case 2  %-150
                            StartOfFor2=round(20+5);
                        case 3  %-30
                            StartOfFor2=round(6+5);
                    end
                    for zz=StartOfFor2:period:maxAll
                        zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                        zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                        
                        if mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0
                            %                             (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz];    TwoDExpre(end+1,:)=[NegCyl zz];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz+1];  TwoDExpre(end+1,:)=[NegCyl zz+1];
                            TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz , EdgeList{edge,3}(1,:) zz  ];
                            TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz+1 , EdgeList{edge,3}(2,:) zz+1 ];
                            
                        end
                    end
                end
            end
            obj.AllScafXover= TwoDExpre ;
        end  % end of function getScafXoverAsStap
        
        function getScafXoverAsStap_every13(obj)
            ConnectList=obj.ForcedConnectList;
            ConnectList2 =[ConnectList(:,1:3) ;ConnectList(:,4:6)] ;
            InitialStappEnds=zeros(2*max(obj.RelateTable(:,3)),2);
            for i=1:size(InitialStappEnds,1)/2
                InitialStappEnds(2*i-1:2*i,1)= obj.RelateTable(i,  3);   % C3 notation
                Bundle=obj.containBundle{obj.RelateTable(i,  1)};
                Cyl=obj.RelateTable(i,  2);
                InitialStappEnds(2*i-1,2)= Bundle.Zbase1(Cyl);
                InitialStappEnds(2*i,2)=Bundle.Zbase2(Cyl);
            end
            
            AdjM=obj.CylinderC5AdjM;    %use CylAdjM in hyperbundle
            if strcmp(obj.ScafOption.prevStack ,'polyT')
                clearance=15;
            else
                clearance=8;
            end
            [U,V]=find(AdjM~=0);   %means cylinder in C5
            U2=U;  %means cylinder in C5
            V2=V;
            Edgeof5Index=  unique([U2 V2],'rows');  %C5
            k=1;
            EdgeList=cell(size(Edgeof5Index,1)/2,3);
            EdgeType=zeros(size(Edgeof5Index,1)/2,1);   % 1->SQ ,  2->HC
            for i=1:size(Edgeof5Index,1)
                Ind=find(Edgeof5Index(i,1)==obj.RelateTable(:,5));Ind=Ind(1);  % Ind
                Ind2=find(Edgeof5Index(i,2)==obj.RelateTable(:,5));Ind2=Ind2(1);
                Bundle=obj.containBundle{obj.RelateTable(Ind,1)};
                if ~xor(ismember( obj.RelateTable(Ind,2)   ,Bundle.AGroup) ,Bundle.AGroupGoUp==1)
                    EdgeList{k,1}=[ Edgeof5Index(i,1) , Edgeof5Index(i,2)];   %in 5rd index;
                    EdgeList{k,3}=[ obj.RelateTable(Ind,1:2)  ; obj.RelateTable(Ind2,1:2)    ];   %in 5rd index;
                    
                    Bundlei=obj.RelateTable(Ind,1);
                    %                     obj.containBundle{Bundlei}.
                    if  strcmp( obj.containBundle{Bundlei}.Lattice, 'Square')
                        EdgeType(k)=1;
                    elseif strcmp( obj.containBundle{Bundlei}.Lattice, 'HoneyComb')
                        EdgeType(k)=2;
                    end
                    k=k+1;
                end
                
            end
            %--------------
            InitialStappEndsC5=InitialStappEnds;
            for iInt=1:size(InitialStappEnds,1)
                InitialStappEndsC5(iInt,1)= obj.RelateTable(InitialStappEndsC5(iInt,1),5);
            end
            %-----------------
            DCellList=obj.ScafdigitSQ;
            G2=[26 27]+16  ; %West
            %             G2=[15 16]+16  ; %West Ori
            
            G4=[31 32]+16-32   ;%East
            
            G3=[40 41]+16 -1  ;%South
            %              G3=[7 8]+16   ;%South Ori
            
            G1=[13 14]+16   ; %North
            %             G1=[23 24]+16   ; %North Ori
            
            TwoDExpre=[]; Tol = 8 ;
            for edge=1:size(EdgeList,1)
                PosCyl=EdgeList{edge,1}(1);    %refers global index, C5
                NegCyl=EdgeList{edge,1}(2);
                
                EndsOfPosCyl=InitialStappEndsC5(find(InitialStappEndsC5(:,1)==PosCyl),2);
                EndsOfNegCyl=InitialStappEndsC5(find(InitialStappEndsC5(:,1)==NegCyl),2);
                minP=min(EndsOfPosCyl);
                maxP=max(EndsOfPosCyl);
                minN=min(EndsOfNegCyl);
                maxN=max(EndsOfNegCyl);
                
                minAll=max(minP,minN)+clearance;
                maxAll=min(maxP,maxN)-clearance ;
                EndsOfPosCyl2=EndsOfPosCyl;
                for kk2=1:length(EndsOfPosCyl2)
                    if mod(kk2,2)==1
                        EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)+clearance;
                    else
                        EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)-clearance;
                    end
                end
                
                EndsOfNegCyl2=EndsOfNegCyl;
                for kk2=1:length(EndsOfNegCyl2)
                    if mod(kk2,2)==1
                        EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)+clearance;
                    else
                        EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)-clearance;
                    end
                end
                
                Pind=find(obj.RelateTable(:,5)==PosCyl);Pind=Pind(1);
                Nind=find(obj.RelateTable(:,5)==NegCyl);Nind=Nind(1);
                
                if EdgeType(edge)==1    % this edge is square lattice
                    NWSE=FindNWSE(Pind,Nind,obj.RelateTable(:,6:7));    % for SQ
                    %                     period=32;
                    period=32/3*5 ;
                    switch NWSE
                        case 1
                            StartOfFor=G1(1);
                        case 2
                            StartOfFor=G2(1);
                        case 3
                            StartOfFor=G3(1);
                        case 4
                            StartOfFor=G4(1);
                    end
                    for zz=StartOfFor:period:maxAll
                        zz=round(zz) ;
                        zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                        zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                        if   mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            %                         if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            
                            EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz];    TwoDExpre(end+1,:)=[NegCyl zz];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz+1];  TwoDExpre(end+1,:)=[NegCyl zz+1];
                            [ia,ib]=ismember( EdgeList{1,3}(1,:) ,obj.RelateTable(:,1:2) , 'rows') ;
                            if mod(obj.RelateTable(ib,4) ,2)==0
                                if sum(and(ismember(ConnectList2(:,1:2),EdgeList{edge,3}(2,:),'rows') , abs(zz-ConnectList2(:,3))<Tol ) )==0 % min distance from foce connection
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz , EdgeList{edge,3}(1,:) zz  ];
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz+1 , EdgeList{edge,3}(2,:) zz+1 ];
                                    %                                else
                                    %                                    sdfsf=3
                                end
                            else
                                if sum(and(ismember(ConnectList2(:,1:2),EdgeList{edge,3}(1,:),'rows') , abs(zz-ConnectList2(:,3))<Tol) )==0  % min distance from foce connection
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz , EdgeList{edge,3}(2,:) zz  ];
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz+1 , EdgeList{edge,3}(1,:) zz+1 ];
                                end
                            end
                            
                        end
                    end
                else
                    Bundle= obj.containBundle{ obj.RelateTable(Pind,1)};
                    Pind2=obj.RelateTable(Pind,2);
                    Nind2=obj.RelateTable(Nind,2);
                    [ Dir ] = FindHCDir(Pind2,Nind2,Bundle.CylInplanePosition) ;
                    period=21;
                    switch Dir
                        case 1   %+90
                            StartOfFor2=round(13+5);
                        case 2  %-150
                            StartOfFor2=round(20+5);
                        case 3  %-30
                            StartOfFor2=round(6+5);
                    end
                    for zz=StartOfFor2:period:maxAll
                        zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                        zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                        
                        if mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0
                            %                             (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz];    TwoDExpre(end+1,:)=[NegCyl zz];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz+1];  TwoDExpre(end+1,:)=[NegCyl zz+1];
                            TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz , EdgeList{edge,3}(1,:) zz  ];
                            TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz+1 , EdgeList{edge,3}(2,:) zz+1 ];
                            
                        end
                    end
                end
            end
            obj.AllScafXover= TwoDExpre ;
        end  % end of function getScafXoverAsStap
        
        function getScafXoverAsStap_every18(obj)
            ConnectList=obj.ForcedConnectList;
            ConnectList2 =[ConnectList(:,1:3) ;ConnectList(:,4:6)] ;
            InitialStappEnds=zeros(2*max(obj.RelateTable(:,3)),2);
            for i=1:size(InitialStappEnds,1)/2
                InitialStappEnds(2*i-1:2*i,1)= obj.RelateTable(i,  3);   % C3 notation
                Bundle=obj.containBundle{obj.RelateTable(i,  1)};
                Cyl=obj.RelateTable(i,  2);
                InitialStappEnds(2*i-1,2)= Bundle.Zbase1(Cyl);
                InitialStappEnds(2*i,2)=Bundle.Zbase2(Cyl);
            end
            
            AdjM=obj.CylinderC5AdjM;    %use CylAdjM in hyperbundle
            if strcmp(obj.ScafOption.prevStack ,'polyT')
                clearance=15;
            else
                clearance=8;
            end
            [U,V]=find(AdjM~=0);   %means cylinder in C5
            U2=U;  %means cylinder in C5
            V2=V;
            Edgeof5Index=  unique([U2 V2],'rows');  %C5
            k=1;
            EdgeList=cell(size(Edgeof5Index,1)/2,3);
            EdgeType=zeros(size(Edgeof5Index,1)/2,1);   % 1->SQ ,  2->HC
            for i=1:size(Edgeof5Index,1)
                Ind=find(Edgeof5Index(i,1)==obj.RelateTable(:,5));Ind=Ind(1);  % Ind
                Ind2=find(Edgeof5Index(i,2)==obj.RelateTable(:,5));Ind2=Ind2(1);
                Bundle=obj.containBundle{obj.RelateTable(Ind,1)};
                if ~xor(ismember( obj.RelateTable(Ind,2)   ,Bundle.AGroup) ,Bundle.AGroupGoUp==1)
                    EdgeList{k,1}=[ Edgeof5Index(i,1) , Edgeof5Index(i,2)];   %in 5rd index;
                    EdgeList{k,3}=[ obj.RelateTable(Ind,1:2)  ; obj.RelateTable(Ind2,1:2)    ];   %in 5rd index;
                    
                    Bundlei=obj.RelateTable(Ind,1);
                    %                     obj.containBundle{Bundlei}.
                    if  strcmp( obj.containBundle{Bundlei}.Lattice, 'Square')
                        EdgeType(k)=1;
                    elseif strcmp( obj.containBundle{Bundlei}.Lattice, 'HoneyComb')
                        EdgeType(k)=2;
                    end
                    k=k+1;
                end
                
            end
            %--------------
            InitialStappEndsC5=InitialStappEnds;
            for iInt=1:size(InitialStappEnds,1)
                InitialStappEndsC5(iInt,1)= obj.RelateTable(InitialStappEndsC5(iInt,1),5);
            end
            %-----------------
            DCellList=obj.ScafdigitSQ;
            G2=[51 52]+1 ; %West
            %             G2=[15 16]+16  ; %West Ori
            
            G4=[15 16]   ;%East
            
            G3=[33 34]   ;%South
            %              G3=[7 8]+16   ;%South Ori
            
            G1=[70 71]+1  ; %North
            %             G1=[23 24]+16   ; %North Ori
            
            TwoDExpre=[]; Tol = 8 ;
            for edge=1:size(EdgeList,1)
                PosCyl=EdgeList{edge,1}(1);    %refers global index, C5
                NegCyl=EdgeList{edge,1}(2);
                
                EndsOfPosCyl=InitialStappEndsC5(find(InitialStappEndsC5(:,1)==PosCyl),2);
                EndsOfNegCyl=InitialStappEndsC5(find(InitialStappEndsC5(:,1)==NegCyl),2);
                minP=min(EndsOfPosCyl);
                maxP=max(EndsOfPosCyl);
                minN=min(EndsOfNegCyl);
                maxN=max(EndsOfNegCyl);
                
                minAll=max(minP,minN)+clearance;
                maxAll=min(maxP,maxN)-clearance ;
                EndsOfPosCyl2=EndsOfPosCyl;
                for kk2=1:length(EndsOfPosCyl2)
                    if mod(kk2,2)==1
                        EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)+clearance;
                    else
                        EndsOfPosCyl2(kk2)=EndsOfPosCyl2(kk2)-clearance;
                    end
                end
                
                EndsOfNegCyl2=EndsOfNegCyl;
                for kk2=1:length(EndsOfNegCyl2)
                    if mod(kk2,2)==1
                        EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)+clearance;
                    else
                        EndsOfNegCyl2(kk2)=EndsOfNegCyl2(kk2)-clearance;
                    end
                end
                
                Pind=find(obj.RelateTable(:,5)==PosCyl);Pind=Pind(1);
                Nind=find(obj.RelateTable(:,5)==NegCyl);Nind=Nind(1);
                
                if EdgeType(edge)==1    % this edge is square lattice
                    NWSE=FindNWSE(Pind,Nind,obj.RelateTable(:,6:7));    % for SQ
                    %                     period=32;
                    period=32/3*7 ;
                    switch NWSE
                        case 1
                            StartOfFor=G1(1);
                        case 2
                            StartOfFor=G2(1);
                        case 3
                            StartOfFor=G3(1);
                        case 4
                            StartOfFor=G4(1);
                    end
                    for zz=StartOfFor:period:maxAll
                        zz=round(zz) ;
                        zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                        zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                        if   mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            %                         if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            
                            EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz];    TwoDExpre(end+1,:)=[NegCyl zz];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz+1];  TwoDExpre(end+1,:)=[NegCyl zz+1];
                            [ia,ib]=ismember( EdgeList{1,3}(1,:) ,obj.RelateTable(:,1:2) , 'rows') ;
                            if mod(obj.RelateTable(ib,4) ,2)==0
                                if sum(and(ismember(ConnectList2(:,1:2),EdgeList{edge,3}(2,:),'rows') , abs(zz-ConnectList2(:,3))<Tol ) )==0 % min distance from foce connection
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz , EdgeList{edge,3}(1,:) zz  ];
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz+1 , EdgeList{edge,3}(2,:) zz+1 ];
                                    %                                else
                                    %                                    sdfsf=3
                                end
                            else
                                if sum(and(ismember(ConnectList2(:,1:2),EdgeList{edge,3}(1,:),'rows') , abs(zz-ConnectList2(:,3))<Tol) )==0  % min distance from foce connection
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz , EdgeList{edge,3}(2,:) zz  ];
                                    TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz+1 , EdgeList{edge,3}(1,:) zz+1 ];
                                end
                            end
                            
                        end
                    end
                else
                    Bundle= obj.containBundle{ obj.RelateTable(Pind,1)};
                    Pind2=obj.RelateTable(Pind,2);
                    Nind2=obj.RelateTable(Nind,2);
                    [ Dir ] = FindHCDir(Pind2,Nind2,Bundle.CylInplanePosition) ;
                    period=21;
                    switch Dir
                        case 1   %+90
                            StartOfFor2=round(13+5);
                        case 2  %-150
                            StartOfFor2=round(20+5);
                        case 3  %-30
                            StartOfFor2=round(6+5);
                    end
                    for zz=StartOfFor2:period:maxAll
                        zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                        zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                        
                        if mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0
                            %                             (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied
                            EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz];    TwoDExpre(end+1,:)=[NegCyl zz];
                            %                             TwoDExpre(end+1,:)=[PosCyl zz+1];  TwoDExpre(end+1,:)=[NegCyl zz+1];
                            TwoDExpre(end+1,:)=[EdgeList{edge,3}(2,:) zz , EdgeList{edge,3}(1,:) zz  ];
                            TwoDExpre(end+1,:)=[EdgeList{edge,3}(1,:) zz+1 , EdgeList{edge,3}(2,:) zz+1 ];
                            
                        end
                    end
                end
            end
            obj.AllScafXover= TwoDExpre ;
        end  % end of function getScafXoverAsStap
        
        
        
        
        function obj=ConvertStap(obj,LaticeType)
            %              RelateTable=obj.RTable;
            if strcmp(LaticeType,'Square')
                MM=size(obj.ScafdigitSQ{1},1);
            elseif strcmp(LaticeType,'Honeycomb')
                MM=size(obj.ScafdigitHC{1},1);
            end
            
            NoOfCyls=max(obj.RelateTable(:,5));
            StappCellList=cell(NoOfCyls,1);
            for i=1:NoOfCyls   %create shell to store information
                StappCellList{i}=-1*ones(MM,4);
            end
            QQ=obj.RelateTable(:,6:7);
            heads=[];
            for j=1:size(obj.StapList3,1)
                C0=-1  ;Z0=-1;
                StappStrand2=  obj.StapList3{j};
                FirstCyl=StappStrand2(1,1);
                NANOCylnum=obj.RelateVec(FirstCyl);
                if mod(NANOCylnum,2)==1 && StappStrand2(2,2)>StappStrand2(1,2)  %% classify into 5' & 3'
                    StappStrand=StappStrand2;
                elseif mod(NANOCylnum,2)==0 && StappStrand2(2,2)<StappStrand2(1,2)
                    StappStrand=StappStrand2;
                elseif mod(NANOCylnum,2)==0 && StappStrand2(2,2)>StappStrand2(1,2)
                    StappStrand=flip(StappStrand2);
                elseif mod(NANOCylnum,2)==1 && StappStrand2(2,2)<StappStrand2(1,2)
                    StappStrand=flip(StappStrand2);
                end
                heads=[heads;StappStrand(1,1) StappStrand(1,2)] ;
                for k=1:2:size(StappStrand,1)-1
                    CylinderIndex=StappStrand(k,1);   %edit 12/2
                    if StappStrand(k,2)>StappStrand(k+1,2)     %move down
                        for z=StappStrand(k,2):-1:StappStrand(k+1,2)
                            if z== StappStrand(k+1,2)
                                if k+2>size(StappStrand,1)
                                    StappCellList{CylinderIndex}(z,:)=[C0,Z0,-1,-1];
                                    continue
                                end
                                C1 =StappStrand(k+2,1)  ;
                                ZNext=StappStrand(k+2,2) ;
                                StappCellList{CylinderIndex}(z,:)=[C0,Z0,C1,ZNext];
                                C0=CylinderIndex;
                                Z0=z;
                            else
                                C1=CylinderIndex;
                                ZNext=z-1;
                                StappCellList{CylinderIndex}(z,:)=[C0,Z0,C1,ZNext];
                                C0=CylinderIndex;
                                Z0=z;
                            end
                        end
                    else    %move up
                        for z=StappStrand(k,2):1:StappStrand(k+1,2)
                            if z== StappStrand(k+1,2)
                                if k+2>size(StappStrand,1)
                                    StappCellList{CylinderIndex}(z,:)=[C0,Z0,-1,-1];
                                    continue
                                end
                                C1 =StappStrand(k+2,1)  ;
                                ZNext=StappStrand(k+2,2) ;
                                StappCellList{CylinderIndex}(z,:)=[C0,Z0,C1,ZNext];
                                C0=CylinderIndex;
                                Z0=z;
                            else
                                C1=CylinderIndex;
                                ZNext=z+1;
                                StappCellList{CylinderIndex}(z,:)=[C0,Z0,C1,ZNext];
                                C0=CylinderIndex;
                                Z0=z    ;
                            end
                        end
                    end
                end
            end
            if strcmp(LaticeType,'Square')
                obj.DigitStapSQ=StappCellList;
                obj.HeadOfStap=heads;
            elseif strcmp(LaticeType,'Honeycomb')
                obj.DigitStapHC=StappCellList;
            end
        end %end of ConvertStap
        
        
        function [h_bingraph,Graph,ST_list ]=FindStapGraph(obj)
            %             tic
            ST_list = zeros(10000,2) ; cc= 1 ;
            for i= 1 : length(obj.StapList3)
                I_stap = obj.StapList3{i} ;
                for j= i+1: length(obj.StapList3)
                    J_stap = obj.StapList3{j} ;
                    if J_stap(1,1) == I_stap(end,1) && abs( J_stap(1,2)-I_stap(end,2))==1
                        ST_list(cc,:) = [i ,j] ; cc=cc+1;
                    end
                    if I_stap(1,1) == J_stap(end,1) && abs( I_stap(1,2)-J_stap(end,2))==1
                        ST_list(cc,:) = [j ,i] ;cc=cc+1;
                    end
                end
            end
            ST_list=ST_list(1:cc-1, :) ;
            Graph = digraph(ST_list(:,1),ST_list(:,2)) ;
            NodesTobeAdded =max(max(ST_list))+1 : i ;
            for k=1: length(NodesTobeAdded)
                Graph = addnode(Graph,1) ;
            end
            
            f552=figure(552);clf ; h_bingraph = plot(Graph); axis off;
            f552.Name='Staple Graph' ; set(f552,'NumberTitle','off');
            bins = conncomp(Graph) ;
            obj.StapGraphBin = bins ;
            %             toc
        end
        
        function obj=FindStapStep2(obj)  %---------
            %             StapleCell=obj.StapList;   % already convert to C5 in FindStapStep  Nov 14 2018
            %             % convert to C5
            %             for k=1:length(StapleCell)
            %                 TT=StapleCell{k};
            %                 for k2=1:size(   TT,1)
            %                     TT(k2,1)= obj.RelateTable( TT(k2,1),5);
            %                 end
            %                 StapleCell{k}=TT;
            %             end
            %---
            
            ScafAndStapClearance=8; % hard, may change for overhang design
            GetScafXOver(obj);
            obj=UpdateStapBP(obj,ScafAndStapClearance);
            
            %             AppOrder=randperm(size(obj.stapBP,1)/4);
            NofOne=0;
            %             HaveAppliedXovers=zeros(size(AppOrder));
            nwhile=1; TakeEasy=0;
            while  1
                StapleCell=obj.StapList;   % already convert to C5 in FindStapStep  Nov 14 2018
                AppOrder=randperm(size(obj.stapBP,1)/4);  HaveAppliedXovers=zeros(size(AppOrder));
                FailedBP = [];
                for i=1:size(obj.stapBP,1)/4
                    if HaveAppliedXovers(AppOrder(i))==0
                        InputM= obj.stapBP(4*AppOrder(i)-3:4*AppOrder(i),:);
                        [ whichstrands,NStapleCell,NoOfBreakStrand] = AddXoverInStapple(StapleCell,InputM );
                        NoOfBreakStrand;
                        
                        if NoOfBreakStrand==2  || TakeEasy==1  % accept the applying Xover
                            StapleCell=  NStapleCell;
                            HaveAppliedXovers(AppOrder(i))=1;
                            NofOne=NofOne+1  ;
                        else
                            FailedBP= union(FailedBP,i) ;
                        end
                        
                    end
                end
                
                
                for reTry =1 :10
                    HaveNotApplied = find(HaveAppliedXovers==0) ; 
                    AppOrder2=randperm(length(HaveNotApplied));
                    for i=1:length(HaveNotApplied)
                            BPind = 4*HaveNotApplied(AppOrder2(i))-3:4*HaveNotApplied(AppOrder2(i)) ;
                            InputM= obj.stapBP(BPind,:);
                            [ whichstrands,NStapleCell,NoOfBreakStrand] = AddXoverInStapple(StapleCell,InputM );
                            NoOfBreakStrand;                            
                            if NoOfBreakStrand==2  || TakeEasy==1  % accept the applying Xover
                                StapleCell=  NStapleCell;
                                HaveAppliedXovers(AppOrder(i))=1;
                                NofOne=NofOne+1  ;
                            else
                                FailedBP= union(FailedBP,i) ;
                            end
                    end                   
                end
                                
                if sum(HaveAppliedXovers)==length(HaveAppliedXovers)  ||nwhile> 10
                    break;
                end
                
                if nwhile==6
                    TakeEasy=1;
                end
                nwhile=nwhile+1 ;
            end
            tt=StapleCell;
            for ii=1:5
                [ StapleCell ] = CalibStapDir( StapleCell,obj.RelateVec);
                [ StapleCell ] = OrganizeStapList( StapleCell,obj ,tt);
            end
            %----check loop staple
            Savek =[];
            for k=1:length(StapleCell)
                stp= StapleCell{k};  stpO=stp ;
                [~,indA]=ismember(stp(1,:),obj.stapBP,'rows');
                [~,indB]=ismember(stp(end,:),obj.stapBP,'rows');
                if ceil(indA/4)==ceil(indB/4) && indA~=0  && indB~=0
%                     k
                    nw=1;
                    while 1
                        stp=circshift(stp,2) ;
                        nw=nw+1;
                        if abs(stp(1,2)-stp(2,2))<8  || nw>20
                            break
                        end
                    end
                    
                    if stp(2,2)>stp(1,2)
                        KK=[ stp(1,:); round(mean(stp(1:2,:))+[0 -1] ) ;round(mean(stp(1:2,:))) ;stp(2:end,:)]  ;
                    else
                        KK=[ stp(1,:); round(mean(stp(1:2,:))+[0 0] ) ;round(mean(stp(1:2,:)))+[0 -1] ;stp(2:end,:)];
                    end
                    StapleCell{k}=[KK(3:end,:)  ;KK(1:2,:) ];
                    
                    if  StapleCell{k}(1,1)==StapleCell{k}(2,1) && StapleCell{k}(1,2)==StapleCell{k}(2,2)
% %                         StapleCell{k}= circshift(StapleCell{k}(2:end-1,:) ,1) ;  % 
                           k
%                         stpO
%                         Savek= union(Savek, k) ;
                        StapleCell{k} =stpO ;
%                                          fprintf(' k = %i ------------\n' , k) ;
                            obj.printStp_C4( StapleCell{k}) ;
                    end
                        
                end
            end
            %-------finish loop
            
            %             ddfgdg=1234
            %             [ NewStapList] = CheckCycleStapAndBreak( NewStapList,obj.stapBP );
            %             GetHyperB.ConnectStapleIf_0ntOnScaf ;  % optional, if ssDNA scaffold assigned as 0, connect them on staple
            scriptcase=1;%----------------halfCrossover
            if strcmp(obj.ScafOption.prevStack,'scafloop')  && strcmp(obj.StapOption.halfXover,'yes')
                AddHalfXoverV3LargeCrossSec;   %input: StapleCell
            end
            %               AddHalfXoverV3LargeCrossSec;   %input: StapleCell
            %             output:StapleCell%   Hard Code
            

             
            %             stapBPExtra=[];
            %             scriptcase=2;%--------------
            %             AddHalfXoverV3LargeCrossSec ;
            
            
            obj.StapList2=StapleCell;   %use C5
            %             obj.StapList2=StapleCell;
        end %end of FindStapStep2
        function printStp_C4(obj, C5)
%             C5= obj.StapList3{stpInd} ;
            
            for k = 1: size(C5 ,1)
                CC = obj.RelateTable(  obj.RelateTable(:,5)==C5(k,1) ,4) ;
               fprintf(' %i %i \n' ,CC ,C5(k,2) ) ;
            end
            
        end
        
        
        
        function ConnectStapleIf_0ntOnScaf(obj)
            fprintf('Connecting staples if ssDNA on Scaf=0 nt. \n ')
            ScafRC5 =  obj.scafC5{1} ;
            for k = 2 : length(obj.scafC5)
                ScafRC5= [ScafRC5 ;  obj.scafC5{k} ] ;
            end
            [aa,bb] = cellfun(@size,obj.scafC5) ;
            
            ConnectedStapleList = [-1,-1];
            for i=1: length(obj.StapList)
                for j=1:length(obj.StapList)
                    if i==j ; continue;end
                    Headi_Tailj = [obj.StapList{i}(1,:) ; obj.StapList{j}(end,:) ] ;
                    [a,b] = ismember(Headi_Tailj, ScafRC5 ,'rows' ) ;   % check head of i and tail of j staple are consecutive in scaf R
                    if sum(a)==2 && diff(b)==1 &&  sum( ismember(b, cumsum(aa)) )==0  % make  sure is not connected due to multi-scaffold. 08/12/2019
                        ConnectedStapleList=union(ConnectedStapleList , [i j],'rows') ;
                    end
                end
            end
            ConnectedStapleList=setdiff(ConnectedStapleList,[-1 -1],'rows') ;
            G = digraph(ConnectedStapleList(:,1)', ConnectedStapleList(:,2)') ;
            %             figure; p=plot(G) ;
            bins= conncomp(G,'Type','weak') ;
            bin= conncomp(G,'OutputForm','cell','Type','weak') ;  ApplyOrder=[];
            for k=1:length(bin)
                if length(bin{k})>1
                    Ntemp=zeros(length(bin{k}),1 ) ;
                    for j=1:length(bin{k})
                        Ntemp(j) = length(dfsearch(G,bin{k}(j)) )  ;  % length of dfsearch to find head and correct order
                    end
                    Root= bin{k}(find(Ntemp==max(Ntemp) ) ) ;  Root=Root(1) ; % prevent loops
                    OneBranch = dfsearch(G,Root) ;
                    ApplyOrder=[ApplyOrder ;[OneBranch(1:end-1),OneBranch(2:end)] ];
                end
            end
            %             MissingDueToLoop =setdiff(ConnectedStapleList,ApplyOrder,'rows') ; % loops with decreace by one
            %             ApplyOrder=[ApplyOrder;MissingDueToLoop] ;
            OriSaveStap=obj.StapList ;
            SaveStap=obj.StapList ;
            for k=1:size(ApplyOrder,1)  % two way merging
                iStap=SaveStap{ApplyOrder(k,1)} ;   jStap=SaveStap{ApplyOrder(k,2)} ;
                NewMergedStap=[jStap ;iStap];
                SaveStap{ApplyOrder(k,1)}=NewMergedStap ;      SaveStap{ApplyOrder(k,2)}=NewMergedStap ;
            end
            IndNodes = zeros(size(bin));  % two-way's merging staples. Use this to extract non-repeated cells.
            for k2=1:length(IndNodes)
                Inds=find(bins==k2) ;
                [aa,~]= cellfun(@size,SaveStap(Inds)) ;
                TakeLargest=find(aa==max(aa));
                IndNodes(k2)=Inds(TakeLargest(1) );
            end
            MissOriginalStap = setdiff(1:length(OriSaveStap), 1:max(max(ConnectedStapleList)));
            IndNodes=union(IndNodes,MissOriginalStap);
            NewStapSet= SaveStap(IndNodes)  ;
            obj.StapList=NewStapSet;
            %----checking
            %             for k=1:length(obj.StapList)
            %             BaseRout = interpolateBase( obj.StapList{k} ) ;
            %            NoRepeat=  size(BaseRout,1)==size(unique(BaseRout,'rows') ,1);
            %            fprintf('Staple %i has  NoRepeat= %i \n',k , NoRepeat)
            %             end
        end % end of ConnectStapleIf_0ntOnScaf
        
        function ConnectStaple3If_0ntOnScaf(obj)
            fprintf('Connecting staples if ssDNA on Scaf=0 nt. \n ')
            ScafRC5 =  obj.scafC5 ;
            ConnectedStapleList = [-1,-1];
            for i=1: length(obj.StapList3)
                for j=1:length(obj.StapList3)
                    if i==j ; continue;end
                    Headi_Tailj = [obj.StapList3{i}(1,:) ; obj.StapList3{j}(end,:) ] ;
                    [a,b] = ismember(Headi_Tailj, ScafRC5 ,'rows' ) ;   % check head of i and tail of j staple are consecutive in scaf R
                    if sum(a)==2 && diff(b)==1
                        ConnectedStapleList=union(ConnectedStapleList , [i j],'rows') ;
                    end
                end
            end
            ConnectedStapleList=setdiff(ConnectedStapleList,[-1 -1],'rows') ;
            G = digraph(ConnectedStapleList(:,1)', ConnectedStapleList(:,2)') ;
            %             figure; p=plot(G) ;
            bins= conncomp(G,'Type','weak') ;
            bin= conncomp(G,'OutputForm','cell','Type','weak') ;  ApplyOrder=[];
            for k=1:length(bin)
                if length(bin{k})>1
                    Ntemp=zeros(length(bin{k}),1 ) ;
                    for j=1:length(bin{k})
                        Ntemp(j) = length(dfsearch(G,bin{k}(j)) )  ;  % length of dfsearch to find head and correct order
                    end
                    Root= bin{k}(find(Ntemp==max(Ntemp) ) ) ;  Root=Root(1) ; % prevent loops
                    OneBranch = dfsearch(G,Root) ;
                    ApplyOrder=[ApplyOrder ;[OneBranch(1:end-1),OneBranch(2:end)] ];
                end
            end
            %             MissingDueToLoop =setdiff(ConnectedStapleList,ApplyOrder,'rows') ; % loops with decreace by one
            %             ApplyOrder=[ApplyOrder;MissingDueToLoop] ;
            OriSaveStap=obj.StapList3 ;
            SaveStap=obj.StapList3 ;
            for k=1:size(ApplyOrder,1)  % two way merging
                iStap=SaveStap{ApplyOrder(k,1)} ;   jStap=SaveStap{ApplyOrder(k,2)} ;
                NewMergedStap=[jStap ;iStap];
                SaveStap{ApplyOrder(k,1)}=NewMergedStap ;      SaveStap{ApplyOrder(k,2)}=NewMergedStap ;
            end
            IndNodes = zeros(size(bin));  % two-way's merging staples. Use this to extract non-repeated cells.
            for k2=1:length(IndNodes)
                Inds=find(bins==k2) ;
                [aa,~]= cellfun(@size,SaveStap(Inds)) ;
                TakeLargest=find(aa==max(aa));
                IndNodes(k2)=Inds(TakeLargest(1) );
            end
            MissOriginalStap = setdiff(1:length(OriSaveStap), 1:max(max(ConnectedStapleList)));
            IndNodes=union(IndNodes,MissOriginalStap);
            NewStapSet= SaveStap(IndNodes)  ;
            obj.StapList3=NewStapSet;
            %----checking
            %             for k=1:length(obj.StapList)
            %             BaseRout = interpolateBase( obj.StapList{k} ) ;
            %            NoRepeat=  size(BaseRout,1)==size(unique(BaseRout,'rows') ,1);
            %            fprintf('Staple %i has  NoRepeat= %i \n',k , NoRepeat)
            %             end
        end % end of ConnectStaple3If_0ntOnScaf  %input: obj.StapList3 output: obj.StapList3
        
        function obj=GetScafXOver(obj)    %Get property: ScafXover  %edit from PartBundle 4/8
            
            PseudoOneScaf =  obj.ScafRouting{1} ;
            for k = 2: length(obj.ScafRouting)
                PseudoOneScaf= [PseudoOneScaf ;  obj.ScafRouting{k} ] ;
            end
            
            Seq=[PseudoOneScaf(1:end-1,:) PseudoOneScaf(2:end,:)];
            SeqGind=zeros(size(Seq,1),4);
            for i2=1:size(Seq,1)
                [~,A1]=ismember(Seq(i2,1:2), obj.RelateTable(:,1:2),'rows');
                SeqGind(i2,1) =  A1 ;
                SeqGind(i2,2) =  Seq(i2,3)  ;
                [~,A2]=ismember(Seq(i2,4:5), obj.RelateTable(:,1:2),'rows');
                SeqGind(i2,3) = A2 ;
                SeqGind(i2,4) = Seq(i2,6)  ;
            end
            
            BPL=obj.stapBPinCell;
            MCyl=max(max(SeqGind(:,[1 3])));
            Corner=cell(MCyl,1);
            MidCorner=cell(MCyl,1);
            SeqP=[SeqGind(:,[1 2]) ;SeqGind(end,3) SeqGind(end,4)];
            for i=1:size(SeqGind,1)
                cly= SeqGind(i,1);
                Position=SeqGind(i,2);
                Corner{cly}=union(Corner{cly}, Position);
                if i==size(SeqGind,1)
                    cly= SeqGind(i,3);
                    Position=SeqGind(i,4);
                    Corner{cly}=union(Corner{cly}, Position);
                end
            end
            
            for j=1:size(MidCorner,1)
                if  mod(length(Corner{j}),2)==0 &&  length(Corner{j})>2
                    MidCorner{j}=Corner{j}(2:end-1);
                end
            end
            
            Result=cell(size(BPL));
            for jj=1:size(BPL,1)
                Result{jj,1}=BPL{jj,1};
            end
            for w=2:2:size(SeqP,1)-2
                cyl1= SeqP(w,1);
                cyl2= SeqP(w+1,1);
                Posi=SeqP(w,2);
                if ismember(Posi,MidCorner{cyl1})
                    for q=1:size(Result,1)
                        if length(union([cyl1 cyl2],Result{q,1}) )==2
                            Result{q,2}=union( Result{q,2},Posi);
                        end
                    end
                end
            end
            obj.ScafXover=Result;
            
        end %end of getscafXover
        
        
        function obj=UpdateStapBP(obj,ScafAndStapClearance)
            OSBP=obj.stapBP; size(OSBP);
            HaveScafXover=[];
            for ScafEdge=1:size(obj.ScafXover,1)
                if ~isempty(obj.ScafXover{ScafEdge,2})~=0
                    HaveScafXover=union(HaveScafXover,ScafEdge);
                end
            end
            k=0;
            for i=1:length(HaveScafXover)
                BreakPonEdge=obj.ScafXover{HaveScafXover(i),2};
                Edge=obj.ScafXover{HaveScafXover(i),1};
                for j=1:4:size(obj.stapBP,1)
                    XMatrix=obj.stapBP(j:j+3,:);
                    XEdge=union(XMatrix(:,1) ,[]);
                    ZPosition=mean(XMatrix(:,2));
                    if nnz(ismember(XEdge,Edge))==2
                        delta= sort(abs(ZPosition-BreakPonEdge));
                        if delta(1)<=ScafAndStapClearance
                            OSBP(j:j+3,:)=0;
                            k=k+1  ;
                        end
                    end
                end
            end
            OSBP((OSBP(:,1)==0),:)=[];
            obj.stapBP=OSBP;
        end %end of UpdatStapBP
        
        function BundleGroupWithSameCylPosition= InplanePairingMatch(obj)
            AdjB = zeros( length(obj.containBundle), length(obj.containBundle)) ;
            for Bi = 1: length(obj.containBundle)
                for Bj = Bi+1: length(obj.containBundle)
                    if obj.BundleHasCylinder(Bi) == obj.BundleHasCylinder(Bj)
                        ClyBi =  obj.containBundle{Bi}.CylInplanePosition ;
                        ClyBj =  obj.containBundle{Bj}.CylInplanePosition ;
                        Compare =  ClyBi==ClyBj ;
                        if sum(sum(Compare)) == numel(Compare)
                            AdjB(Bi, Bj ) =1 ;           AdjB(Bj, Bi ) =1 ;
                        end
                    end
                end
            end
            G = graph(AdjB);
            BundleGroupWithSameCylPosition = conncomp(G) ;
        end
        
        function obj=CreateTFWindow(obj,PremPair,varargin)
            %             tic
            if ~isempty(varargin)
                AssignAdjM=varargin{1};nargin;
                if nargin>3
                    if  ~isempty(varargin{3})
                        ftab= varargin{3} ; figure(ftab.Parent.Parent.Parent.Parent) ;
                        fH=gcf;
                    end
                else
                    fH=gcf;
                end
            else
                AssignAdjM=[];
            end
            
            rotate3d off;
            ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
            %             % initialize uicontrol in ss_Assembly
            ax2=findobj(fH,'Tag','AssemblyGraph') ;
            %             axes(ax2) ;
            ax2.Parent= ss_Assembly ;
            
            btn_IntLoad=findobj(gcf,'String','Initialize') ; btn_IntLoad.Visible='off';
            GO_inthistab= findobj( findobj(gcf,'Tag','ss_Assembly'),'Type','UIcontrol' ) ;
            for k=1:length(GO_inthistab)
                if ~isfield(GO_inthistab(k).UserData,'keepme');   delete(GO_inthistab(k)) ;   end
            end
            BG = findobj( findobj(gcf,'Tag','ss_Assembly'),'Type','uibuttongroup' ) ;
            for k=1:length(BG)
                delete(BG(k)) ;
            end
            TG=  findobj( findobj(gcf,'Tag','ss_Assembly'),'Type','uitable' ) ;
            for k=1:length(TG)
                delete(TG(k)) ;
            end
            
            if ~isempty(varargin)
                AssignAdjM=varargin{1};
            else
                AssignAdjM=[];
            end
            obj.SavePremPair=PremPair ;
            
            %         fH=figure('name','Assembly window','numbertitle','off');
            fH=ss_Assembly ;  axes(findobj(gcf,'Tag','AssemblyMain'));
            cltab;        axes(findobj(gcf,'Tag','AssemblyMain'));
            hold on; xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
            Str={'--'};
            for k=1:length(obj.containBundle) ;  Str{k}=strcat('Bundle  ',num2str(k));    end
            %---------
            MovingGroup = uibuttongroup(fH,'Position',[0.7 0.85 0.28 0.13],'Title','Bundle manipulation','FontSize',12,'BorderWidth',2) ;
            %             MovingGroup = uibuttongroup(fH,'Position',[0.7 0.85 0.28 0.13],'Title','Moving panel','FontSize',12,'BorderWidth',2) ;
            AddEditSLGroup = uibuttongroup(fH,'Position',[0.7 0.02 0.28 0.12],'Title','Modify components and Save/Load mechanism','FontSize',12,'BorderWidth',2) ;
            BiBjPairGroup = uibuttongroup(fH,'Position',[0.7 0.15 0.28 0.12],'Title','Query and assign individual connection','FontSize',12,'BorderWidth',2) ;
            AdjMGroup = uibuttongroup(fH,'Position',[0.7 0.28 0.28 0.56],'Title','Connectivity','FontSize',12,'BorderWidth',2) ;
            uistack(AdjMGroup,'bottom') ;
            %             MovingGroup.BackgroundColor=0.9*ones(1,3) ;   AddEditSLGroup.BackgroundColor=0.98*ones(1,3) ;
            align([MovingGroup AdjMGroup BiBjPairGroup AddEditSLGroup ],'Left','distribute');
            popupH = uicontrol(MovingGroup,'Style', 'listbox',...
                'String', Str,'Unit','normalized','Position', [0.05 0.6 0.23 0.35],'Max',length(obj.containBundle));    % [0.7 0.9 0.1 0.05]   %change to listbox
            s = sprintf('Button tooltip line 1\nButton tooltip line 2');
            txtBundle = uicontrol(fH,'Style','text','Unit','normalized','Position', [0.82 0.9 0.18 0.08],....
                'String','Select Bundle'  ,'FontUnits','normalized','FontSize',0.5,'Tooltip',s,'Visible','off');
            %         txtBundle = uicontrol(fH,'Style','text','Unit','normalized','FontSize',14,'Position', [0.82 0.9 0.1 0.08],....
            %         'String','Select Bundle'   );
            
            checkH = uicontrol(MovingGroup,'Style', 'checkbox','String', 'Trans. / Rotate','Unit','normalized','Position', [0.05 0.1 0.25 0.4] ,'FontSize',12); % [0.7 0.85 0.07 0.05]
            checkH2_GL = uicontrol(MovingGroup,'Style', 'checkbox','String', 'Global / Local','Unit','normalized','Position', [0.68 0.1 0.25 0.4],'FontSize',12); %[0.78 0.85 0.07 0.05]
            editH = uicontrol(MovingGroup,'Style', 'edit','String', '1','Unit','normalized','Position', [0.32 0.12 0.1 0.35]);  % [0.86 0.85 0.05 0.05]
            txtH2 = uicontrol(MovingGroup,'Style','text','Unit','normalized','FontSize',14,'Position', [0.42 0.1 0.2 0.4],....
                'String','Unit: nm','HorizontalAlignment','left');  % [0.92 0.85 0.05 0.05]
            
            checkH.Callback=@(src,evn)checkFcn(obj,src,evn,editH,txtH2)  ;
            tH = uitable(fH,'Data',zeros(6,length(obj.containBundle)),'Unit','normalized','Position',[0.7 0.01 0.25 0.2]);
            tH.Visible='off';
            
            aH=findobj(gcf,'Tag','AssemblyMain')  ;
            axes(aH);
            %                  aH.Position(3)=0.5; aH.Position(1)=0.05;
            %             aH.UserData='ax1';
            
            
            
            %---------
            patchH=cell(3, length(obj.containBundle)) ;
            V0CylinderXYZ=cell(1, length(obj.containBundle)) ;
            N_slice = 100 ;
            for bundlei=1:length(obj.containBundle)
                Bundle=obj.containBundle{bundlei};
                CylInBundles=1:length(Bundle.Z1);
                TwoSide=Bundle.CylinderXYZGlobal ;
                %                 QQ=[];
                BaseIndex=Bundle.findNeiborBaseIndexofCylinder(CylInBundles);
                Hex=Bundle.findHelixWithR(1,CylInBundles,BaseIndex);
                uX= Bundle.TransformMatrix2(1:3,1)';
                uY= Bundle.TransformMatrix2(1:3,2)';
                %              Hex=Bundle.findHelixWithR(3,CylInBundles,BaseIndex);
                
                if strcmp(Bundle.Lattice, 'Square')
                    QQ=zeros(length(Bundle.Zbase1)*N_slice*4,3);
                    dist=0.8;    %distance to define alphashape from cylinder center  2.2
                    
                    for cyli=1:length(Bundle.Zbase1)
                        BottomEnd=TwoSide(cyli,1:3)  ;
                        TopEnd=TwoSide(cyli,4:6)  ;
                        x=[0 ;1];
                        v=[BottomEnd ;TopEnd];
                        xq = linspace(0,1,N_slice);   %each cylinder genearates  4*100 external points
                        vq = interp1(x,v,xq,'linear','extrap');
                        nV=size(vq,1);
                        %                  QQ=[QQ;   [vq + dist*ones(nV,1)*[uX+uY]  ]; [vq+ dist*ones(nV,1)*[uX-uY]];[vq+ dist*ones(nV,1)*[-uX+uY]]  ;[vq+ dist*ones(nV,1)*[-uX-uY]]];
                        QQ(1+4*N_slice*(cyli-1):4*N_slice*cyli,:)=[ (vq + dist*ones(nV,1)*(uX+uY) ) ; (vq + dist*ones(nV,1)*(-uX+uY) ); (vq + dist*ones(nV,1)*(-uX-uY) ) ;(vq + dist*ones(nV,1)*(uX-uY) )];
                    end
                elseif strcmp(Bundle.Lattice, 'HoneyComb')
                    coeff=sqrt(3)/2 ;
                    dist=0.75;    %distance to define alphashape from cylinder center  2.2
                    QQ=zeros(length(Bundle.Zbase1)*N_slice*6,3);
                    for cyli=1:length(Bundle.Zbase1)
                        BottomEnd=TwoSide(cyli,1:3)  ;
                        TopEnd=TwoSide(cyli,4:6)  ;
                        x=[0 ;1];
                        v=[BottomEnd ;TopEnd];
                        xq = linspace(0,1,N_slice);   %each cylinder genearates  4*100 external points
                        vq = interp1(x,v,xq,'linear','extrap');
                        nV=size(vq,1);
                        %                  QQ=[QQ;   [vq + dist*ones(nV,1)*[uX+uY]  ]; [vq+ dist*ones(nV,1)*[uX-uY]];[vq+ dist*ones(nV,1)*[-uX+uY]]  ;[vq+ dist*ones(nV,1)*[-coeff*uX-uY]];  ];
                        %                  QQ(1+600*(cyli-1):400*cyli,:)=[ [vq + dist*ones(nV,1)*[uX+uY]  ]; [vq+ dist*ones(nV,1)*[uX-uY]];[vq+ dist*ones(nV,1)*[-uX+uY]]  ;[vq+ dist*ones(nV,1)*[-uX-uY]]];
                        QQ(1+6*N_slice*(cyli-1):6*N_slice*cyli,:)=  [(vq+ dist*ones(nV,1)*(-uY));(vq+ dist*ones(nV,1)*(+uY))  ;(vq+ dist*ones(nV,1)*(-coeff*uX-0.5*uY));....
                            (vq+ dist*ones(nV,1)*(coeff*uX-0.5*uY));(vq+ dist*ones(nV,1)*(-coeff*uX+0.5*uY));(vq+ dist*ones(nV,1)*(coeff*uX+0.5*uY))];
                    end
                end
                [patchH{2,bundlei},V0CylinderXYZ{bundlei}] =DrawPartH( Bundle,[1 0 1],0,aH,bundlei ) ;
                %             shp=alphaShape(QQ,2,'RegionThreshold',200);
                %           bundlei
                %           try
                %                 t2=toc ;
                %                 size(QQ)
                %                 k = boundary(QQ,0.75);
                %                 QQ=  QQ(unique(k) ,:) ; size(QQ) ;
                %                 shp=alphaShape(QQ,Inf);
                
                
                shp=alphaShape(QQ,2);
                
                %                 shp=alphaShape(QQ,1.3);
                %                 crit = criticalAlpha(shp, 'one-region') ;
                
                %                 t3=toc ;
                patchH{1,bundlei}=plot(shp);   %plot Volume
                %                 t4=toc;
                %---reduce number of faces to increase graphic performance
                TargetFaces= 400;
                R=TargetFaces/size(patchH{1,bundlei}.Faces ,1);
                
                reducepatch(patchH{1,bundlei},R) ;
                
                %------------  May 10
                
                %                  d32 =t3-t2
                %                  d43 =t4-t3
                patchH{1,bundlei}.EdgeColor='none' ;
                patchH{1,bundlei}.FaceAlpha=0.5;
                
                patchH{1,bundlei}.FaceColor=[0.3,0.3,0.3]; % bundle color
                %                 patchH{1,bundlei}.FaceColor(1)=0;
                patchH{1,bundlei}.UserData= patchH{1,bundlei}.FaceColor;
                patchH{1,bundlei}.Tag=num2str(bundlei);
                %-----------find bundle reasonable Connect points
                %                 if strcmp(Bundle.Lattice, 'HoneyComb')
                template=Bundle.template;
                %                 else
                %                 template=[2 ,3, 4, 5,7 8,10 ,11, 12,13,15 16,18 ,19,20 21,23, 24,26,27,28,29,31 ,32];
                %                 end
                ExaHelixTable=[];
                ConnectPoint=[];
                for cylinH=1: length(Hex)
                    HexofCyl=Hex{cylinH};   % [ x, y,z]
                    BaseIndexCyl=BaseIndex{cylinH};
                    %                     scatter3( Hex{cylinH}(:,1),Hex{cylinH}(:,2),Hex{cylinH}(:,3) ,'.b');  %--------------------
                    InorOut= inShape(shp,HexofCyl(:,1),HexofCyl(:,2),HexofCyl(:,3));   %1-> in, 0->out
                    HexofCyl=HexofCyl(InorOut==0,:);
                    BaseIndexCyl=BaseIndexCyl(InorOut==0);
                    
                    YY=mod(BaseIndexCyl,Bundle.period(1)); YY(YY==0)=Bundle.period(1);
                    [~,indY]=ismember(YY,template);
                    TakeYY=zeros(size(indY));
                    for iY=1:length(TakeYY)-1    %check paired
                        if mod(indY(iY),2)==1 && mod(indY(iY+1),2)==0
                            TakeYY(iY:iY+1)=1;
                        end
                    end
                    BaseIndexCyl=BaseIndexCyl(TakeYY==1);
                    HexofCyl=HexofCyl(TakeYY==1,:);
                    %                 plot3( HexofCyl(:,1),HexofCyl(:,2),HexofCyl(:,3) ,'.r');%
                    %                 real base location , in case of explanation
                    fract=[cylinH*ones(length(BaseIndexCyl),1), BaseIndexCyl,HexofCyl];
                    GroupAdd=zeros(size(fract,1)/2,6);   %   [Cyl, Base1,Base12, Gx,Gy,Gz];
                    for nG=1:size(GroupAdd,1)
                        GroupAdd( nG,:)=[ fract(2*nG-1,1),fract(2*nG-1,2),fract(2*nG,2),mean( fract(2*nG-1:2*nG,3:5))];
                    end
                    GroupAdd2=[GroupAdd(:,1:2),GroupAdd(:,1),GroupAdd(:,3),GroupAdd(:,4:6)];
                    %                 ExaHelixTable=[ExaHelixTable;GroupAdd]; % old
                    ExaHelixTable=[ExaHelixTable;GroupAdd2];
                end
                
                PrevPairRes= PremPair.CellPairList{bundlei};
                AddToExaHT=zeros(size(PrevPairRes,1)*2,7);  %side connect points [base, -Cyl1,-Cyl2, Gx,Gy,Gz];
                for iPre=1:size(PrevPairRes)
                    Cly1=PrevPairRes(iPre,1);
                    Cly2=PrevPairRes(iPre,2); tol2=Bundle.Tol;
                    if abs(Bundle.Z1(Cly1)-Bundle.Z1(Cly2))<tol2 && abs(Bundle.Z2(Cly1)-Bundle.Z2(Cly2))<tol2
                        OverBase=2;
                        Z1A=Bundle.Zbase1(Cly1)-OverBase;  Z1B=Bundle.Zbase1(Cly2)-OverBase;
                        Z2A=Bundle.Zbase2(Cly1)+OverBase;  Z2B=Bundle.Zbase2(Cly2)+OverBase;
                        pseudoRadius=0.5;
                        BottomXYZ1=Bundle.findHelixWithR(pseudoRadius,Cly1,{Z1A});BottomXYZ1=BottomXYZ1{1};
                        BottomXYZ2=Bundle.findHelixWithR(pseudoRadius,Cly2,{Z1B});BottomXYZ2=BottomXYZ2{1};
                        TopXYZ1=Bundle.findHelixWithR(pseudoRadius,Cly1,{Z2A});TopXYZ1=TopXYZ1{1};
                        TopXYZ2=Bundle.findHelixWithR(pseudoRadius,Cly2,{Z2B});TopXYZ2=TopXYZ2{1};
                        
                        AddToExaHT(2*iPre-1,:) =[-Cly1,Z1A,-Cly2,Z1B,mean([BottomXYZ1;BottomXYZ2])] ;
                        AddToExaHT(2*iPre,:) =[-Cly1,Z2A,-Cly2,Z2B,mean([TopXYZ1;TopXYZ2])]  ;
                        %                  plot3( BottomXYZ1(:,1),BottomXYZ1(:,2),BottomXYZ1(:,3) ,'.r');%
                        %                   plot3( BottomXYZ2(:,1),BottomXYZ2(:,2),BottomXYZ2(:,3) ,'.r');%
                        %                    plot3( TopXYZ1(:,1),TopXYZ1(:,2),TopXYZ1(:,3) ,'.r');%
                        %                     plot3( TopXYZ2(:,1),TopXYZ2(:,2),TopXYZ2(:,3) ,'.r');%
                    end
                end
                ExaHelixTable=[ExaHelixTable ;AddToExaHT];
                
                ExaHelixTable=ExaHelixTable( ~ismember(ExaHelixTable,zeros(1,7),'rows') ,:) ;
                
                
                patchH{3,bundlei}=scatter3( ExaHelixTable(:,5),ExaHelixTable(:,6),ExaHelixTable(:,7),'k' ,'Marker','o','MarkerFaceColor',[0.3 0.3 0.3],'HitTest','off');
                patchH{3,bundlei}.MarkerFaceAlpha=1;
                %             patchH{3,bundlei}.SizeData=20;
                obj.containBundle{bundlei}.ExternalXoverAsFB=ExaHelixTable ;
            end
            
            
            
            
            axis auto;grid on;
            nBumdle=length(obj.containBundle);
            if nBumdle~=1
                AllCombination= nchoosek(1:length(obj.containBundle),2) ; %
            else
                AllCombination=[];
            end
            fH.UserData.pH=cell(size(AllCombination,1),1);
            %--------- -----
            btn1 = uicontrol(AdjMGroup,'Style', 'pushbutton', 'String', 'FindXover','Unit','normalized', 'Position', [0.4 0.83 0.2 0.15] ,...
                'Callback', {@(src,evn)SortDist(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ)});  % [0.7 0.75 0.05 0.05]
            btn2 = uicontrol(AdjMGroup,'Style', 'pushbutton', 'String', 'ssDNA','Unit','normalized', 'Position', [0.78 0.83 0.2 0.15] ,...
                'Callback', {@(src,evn)ConnectFC(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ,PremPair)}); % [0.77 0.75 0.05 0.05]
            Tg1 = uicontrol(AdjMGroup,'Style', 'togglebutton', 'String', 'Graph','Unit','normalized', 'Position', [0.02 0.83 0.1 0.15] ) ;
            Tg2 = uicontrol(AdjMGroup,'Style', 'togglebutton', 'String', 'Table1','Unit','normalized', 'Position', [0.2 0.83 0.1 0.15] ) ;
            Tg3 = uicontrol(AdjMGroup,'Style', 'togglebutton', 'String', 'Table2','Unit','normalized', 'Position', [0.35 0.83 0.1 0.15] ) ;
            Tg1.Callback=@(src,evn)SwitchVisualGraph(obj,src,evn,Tg1,Tg2,Tg3 ) ;
            Tg2.Callback=@(src,evn)SwitchVisualGraph(obj,src,evn,Tg1,Tg2,Tg3 ) ;
            Tg3.Callback=@(src,evn)SwitchVisualGraph(obj,src,evn,Tg1,Tg2,Tg3 ) ;
            align([Tg1 Tg2 Tg3 btn1 btn2],'distribute','bottom');
            
            
            btn3 = uicontrol(MovingGroup,'Style', 'pushbutton', 'String', 'Orthogonal R','Unit','normalized', 'Position', [0.35 0.55 0.2 0.43] ,...
                'Callback', {@(src,evn)GetOrthogonal(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ,PremPair)},'FontSize',12); % [0.84 0.75 0.05 0.05]
            %             btn4 = uicontrol(fH,'Style', 'pushbutton', 'String', 'Simulate Assembly','Unit','normalized', 'Position', [0.91 0.75 0.05 0.05] ,...
            %                 'Callback', {@(src,evn)SimulateAssembly(obj,src,evn)});
            
            btn5 = uicontrol(AddEditSLGroup,'Style', 'pushbutton', 'String', 'Save','Unit','normalized', 'Position', [0.5 0.1 0.2 0.8] ,...
                'Callback', {@(src,evn)SaveT(obj,src,evn)});    %[0.8 0.55 0.05 0.05]
            btn6 = uicontrol(AddEditSLGroup,'Style', 'pushbutton', 'String', 'Load','Unit','normalized', 'Position', [0.75 0.1 0.2 0.8] ,...
                'Callback', {@(src,evn)LoadT(obj,src,evn)});   % [0.87 0.55 0.05 0.05]
            
            popupBuni = uicontrol(BiBjPairGroup,'Style', 'popup','String', Str,'Unit','normalized','Position', [0.02 0.55 0.15 0.3]);   % [0.7 0.2 0.05 0.05]
            popupBunj = uicontrol(BiBjPairGroup,'Style', 'popup','String', Str,'Unit','normalized','Position', [0.2 0.55 0.15 0.3]); %[0.77 0.2 0.05 0.05]
            popupBuni.ButtonDownFcn=@(src,evn)IgnoreStrBundleLeaveIndex(obj,src,evn) ; popupBuni.UserData.Mode= 1 ;
            popupBunj.ButtonDownFcn=@(src,evn)IgnoreStrBundleLeaveIndex(obj,src,evn) ; popupBunj.UserData.Mode= 1 ;
            
            %             strNum = cellfun(@num2str,num2cell(0:10),'uniformoutput',0);
            textConnQS = uicontrol(BiBjPairGroup,'Style', 'text','String', '# of Connections =','Unit','normalized','Position', [0.15 0.1 0.3 0.3],'Tag','textConnQS','HorizontalAlignment','left'); %[0.73 0.15 0.05 0.05]
            %             if want to directly click on bundle to select
%             for k=1:length(obj.containBundle)
%                 patchH{1,k}.ButtonDownFcn=@(src,evn)ClickPatch(obj,src,evn,popupH,patchH,popupBuni , popupBunj,textConnQS ) ;
%             end
            
            btn7 = uicontrol(BiBjPairGroup,'Style', 'pushbutton', 'String', 'Set(left_C) ','Unit','normalized', 'Position', [0.5 0.1 0.2 0.8]); % [0.84 0.2 0.05 0.05]
            btn8 = uicontrol(BiBjPairGroup,'Style', 'pushbutton', 'String', 'Delete(right_C) ','Unit','normalized', 'Position', [0.75 0.1 0.2 0.8] ); %  [0.91 0.2 0.05 0.05]
            %             btn7.Callback=@(src,evn)QueryConn(obj,src,evn, popupBuni,popupBunj,textConnQS) ;
            popupBuni.Callback=@(src,evn)QueryConn(obj,src,evn, popupBunj,textConnQS,popupH,patchH) ;  %update text about # of Conn
            popupBunj.Callback=@(src,evn)QueryConn(obj,src,evn, popupBuni,textConnQS,popupH,patchH) ;
            
            
            btn7.Callback=@(src,evn)SetConn(obj,src,evn, popupBuni,popupBunj,textConnQS,patchH) ;
            btn8.Callback=@(src,evn)DeleteConn(obj,src,evn, popupBuni,popupBunj,textConnQS,patchH) ;
            
            btn9 = uicontrol(AddEditSLGroup,'Style', 'pushbutton', 'String', 'Delete/Insert bundle ','Unit','normalized', 'Position',  [0.05 0.1 0.2 0.8] );
            btn9.Callback=@(src,evn)DelAddBundle(src,evn,obj,fH,patchH,popupH) ;  %  [0.7 0.03 0.08 0.08]
            
            btn10 = uicontrol(AddEditSLGroup,'Style', 'pushbutton', 'String', 'Edit bundle ','Unit','normalized', 'Position', [0.25 0.1 0.2 0.8] );
            btn10.Callback=@(src,evn)EditBundle(src,evn,obj,fH,patchH ,popupH ) ; % [0.8 0.03 0.08 0.08]
            align([btn9 btn10 btn5 btn6],'distribute','bottom');
            %----------------
            box on
            h1 = light('Position',[-20 20 20],'Style','infinite'); h2 = light('Position',[20 -20 20],'Style','infinite');
            h3 = light('Position',[-20 -20 20],'Style','infinite');h4 = light('Position',[20 20 20],'Style','infinite');
            
            
            
            ax2=findobj(gcf,'Tag','AssemblyGraph') ;
            axes(ax2) ;
            ax2.Parent= AdjMGroup ;  ax2.Position=[ 0.1, 0.1, 0.8 ,0.58];
            %           ax2.Parent=fH; ax2.Position=[0.7,0.3,0.23,0.23];
            %            ax2=axes(fH,'Position',[0.7,0.3,0.25,0.25]);
            hold on;ax2.UserData='ax2';
            xlim([-1,1]); ylim([-1,1])  ;   ax2.XTick=[] ;ax2.YTick=[]; ax2.Box='on';
            theta=linspace(0,2*pi,nBumdle+1); theta=theta(1:end-1);
            [x,y] = pol2cart(theta,0.8);     [x2,y2] = pol2cart(theta,0.9);
            scatter(ax2,x,y,'O');   str=num2cell(1:length(obj.containBundle));
            text(x2,y2,str,'FontSize',16 );
            
            %           AllCombination= nchoosek(1:length(obj.containBundle),2) ; %
            if isempty(AssignAdjM)
                obj.BundleAdjM=zeros(nBumdle,nBumdle);   %default root is 1st bundle
                obj.BundleAdjM(1,2:end)=2;
                obj.BundleAdjM(2:end,1)=2;
            else
                obj.BundleAdjM=AssignAdjM;
            end
            
            obj.LineCombination=AllCombination;
            fH.UserData.BundleConnectivity=[];
            fH.UserData.textNumberHandle=[];
            fH.UserData.HyperBundle=obj;
            set(fH,'defaultLegendAutoUpdate','off');
            for lineAllC=size(AllCombination,1):-1:1
                xx=[x(AllCombination(lineAllC,1)), x(AllCombination(lineAllC,2))];
                yy=[y(AllCombination(lineAllC,1)), y(AllCombination(lineAllC,2))];
                if obj.BundleAdjM(AllCombination(lineAllC,1) ,AllCombination(lineAllC,2)) >=1
                    %                     if strcmp(obj.TakeSeveralV3Type(lineAllC),'auto')
                    fH.UserData.BundleConnectivity{lineAllC}= plot( xx , yy,'-r'  );
                    %                     else
                    %                     fH.UserData.BundleConnectivity{lineAllC}= plot( xx , yy,'Color',[1 0 1 ] );
                    %                     end
                else
                    fH.UserData.BundleConnectivity{lineAllC}= plot(xx ,yy ,':k' );
                end
                fH.UserData.BundleConnectivity{lineAllC}.ButtonDownFcn=@(src,evn)lineselect(obj,src,evn,lineAllC,patchH);
                fH.UserData.textNumberHandle{lineAllC}=text(0.3*xx(1)+ 0.7*xx(2) ,0.3*yy(1)+ 0.7*yy(2) , num2str(obj.BundleAdjM(AllCombination(lineAllC,1),AllCombination(lineAllC,2))),'Color','blue' );
            end
            
            obj.LineOption=ones(size(obj.LineCombination,1),2);
            cellstr=cell(size(obj.LineOption,1),5);
            for k3=1:size(cellstr,1)  % scanning connection setting
                cellstr{k3,1}=obj.LineCombination(k3,1); cellstr{k3,2}='Both';
                cellstr{k3,3}=obj.LineCombination(k3,2); cellstr{k3,4}='Both';
                cellstr{k3,5} = obj.BundleAdjM(obj.LineCombination(k3,1),obj.LineCombination(k3,2) ) ;
            end
            SearchOps={'Both','Side only','End only','Manual'} ;
            tH2 = uitable(AdjMGroup,'Data',cellstr,'Unit','normalized','Position',[0.02 0.05 0.96 0.65],'ColumnEditable',[false true false true true],...
                'ColumnName',{'Bi'; 'opt' ;'Bj';'opt';'# conn' } ,'ColumnFormat',({[] SearchOps [] SearchOps []}) ,'FontSize',10 );  %[0.7 0.6 0.25 0.1]
            tH3 = uitable(AdjMGroup,'Data',cellstr,'Unit','normalized','Position',[0.02 0.05 0.96 0.65],'ColumnEditable',[false true false true true],...
                'ColumnName',{'Bi'; 'opt' ;'Bj';'opt';'# conn' } ,'ColumnFormat',({[] SearchOps [] SearchOps []}),...
                'CellEditCallback',@(src,evn)tH2Change(obj,src,evn,patchH , tH2 ) ,'FontSize',10  );  %[0.7 0.6 0.25 0.1]
            %--------update tH3 Data
            IndsRemain= zeros(size(tH2.Data ,1),1) ;
            for k=1:size(tH2.Data ,1)
                if tH2.Data{k,5}~=0
                    IndsRemain(k)=1;
                end
            end
            tH3.Data=tH2.Data(IndsRemain==1,:) ;
            %------------
            
            tH2.CellEditCallback =@(src,evn)tH2Change(obj,src,evn,patchH,tH3 ) ;
            
            
            fH.UserData.tb2=tH2;
            fH.UserData.tb3=tH3;
            
            axes(aH);
            QQ=[ aH.XLim ; aH.YLim ; aH.ZLim];
            h5 = light('Position', mean(QQ,2),'Style','local');
            checkH_View = uicontrol(MovingGroup, 'Style', 'checkbox','String', 'Axis auto/equal','Unit','normalized','Position', [0.7 0.65 0.28 0.35],'FontSize',12); %[0.7 0.55 0.1 0.05]
            checkH_View.Callback=@(src,evn)checkFcn2(src,evn,aH)  ;
            checkH_Rcenter = uicontrol(MovingGroup, 'Style', 'checkbox','String', 'Rcenter','Unit','normalized','Position', [0.7 0.35 0.28 0.2],'FontSize',12); %[0.7 0.55 0.1 0.05]
            popupH.Callback=@(src,evn)SelectBundlePop(obj,src,evn,patchH)  ;
            textBundle=cell(1,size(patchH,2));
            for k=1:length(textBundle)
                XYZALL= [  patchH{3,k}.XData; patchH{3,k}.YData; patchH{3,k}.ZData];
                xyz=mean(XYZALL,2);
                textBundle{k}=text(xyz(1),xyz(2),xyz(3)+2, strcat(num2str(k)),'FontSize',22,'clipping','on','Color',[0 0.9 1] );
                set( textBundle{k},'Parent',findobj(gcf,'Tag', 'HidenAssemblyMain'));
            end
             align([checkH_View checkH_Rcenter checkH2_GL ],'left','distribute');
            for k=1:length(obj.containBundle)
                patchH{1,k}.ButtonDownFcn=@(src,evn)ClickPatch(obj,src,evn,popupH,patchH,popupBuni , popupBunj,textConnQS ,checkH_Rcenter) ;
            end
            %---unable all toolbars
            %             aT = findall(gcf) ;
            %             bT = findall(aT,'ToolTipString','Save Figure') ;
            %             cT = findall(bT(1).Parent  ,'-not','Tag','FigureToolBar')  ;
            %             for ii = 1:length(cT)
            %                if  strcmp(cT(ii).Enable ,'on')
            %                    if isfield(cT(ii) ,'State')
            %                     cT(ii).State='off' ;
            %                    end
            %                end
            %             end
            rotate3d on ;
            rotate3d off ;
            set(gcf,'KeyPressFcn',{@(src,evn)keyMove(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ,editH,textBundle,checkH2_GL ,txtH2,checkH_Rcenter )}) ;
            fHH=gcf ;
            fHH.UserData.saveKeyMove = get(gcf,'KeyPressFcn') ;
            
            %           setGOfontsize( gctab , 10 , {'UIControl'} )  % set fontsize for uicontrol in this tab
            %--------------create legend for indication & callback for
            %boxing
            axes(findobj(gcf,'Tag','AssemblyMain')) ; ax=gca;
            if isfield(fHH.UserData,'LgIndication')
                if isvalid(fHH.UserData.LgIndication)
                    delete(fHH.UserData.LgIndication) ;
                end
            end
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            ForLegend2=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            ForLegend3=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');    %ForLegend3=scatter3(mean(ax.XLim),mean(ax.YLim),mean(ax.ZLim));
            fHH.UserData.LgIndication= legend([ForLegend,ForLegend2,ForLegend3],'x','y','z' ,'Location','northwest' ) ;
            fHH.UserData.LgIndication.String={'[Q][A] =x';'[W][S] =y';'[E][D] =z'};
            fHH.UserData.LgIndication.Interpreter='tex';  %latex
            fHH.UserData.LgIndication.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            ForLegend2.Marker='.' ; ForLegend2.Marker='none';
            ForLegend3.Marker='.' ; ForLegend3.Marker='none';
            fHH.UserData.LgIndication.ButtonDownFcn=@(src,evn)LegendBoxing( src,evn,ax );
            %             fHH.UserData.LgIndication.Title.Interpreter='latex';
            
            %             drawnow;
            fHH.UserData.LgIndication.Units='normalized'; fHH.UserData.LgIndication.AutoUpdate ='off';
            fHH.UserData.LgIndication.Title.String='Click me for instructions' ;
            fHH.UserData.LgIndication.Position=[0.0063 0.9528 0.1569 0.0387];
            
            %-----------------
            
            ax= findobj(gcf,'Tag','MechOH3D');
            axes(ax); cltab;
            t=findobj(gcf,'Tag','OHTable');  t.Visible='off';
            t.Data='';
            
            TypeStr="auto" ;
            obj.TakeSeveralV3Type= repmat(TypeStr, size(obj.LineCombination,1), 1) ;
            
            AssignIcon( btn1,'FindXOver.jpg' ); btn1.TooltipString = 'Find scaffold connections to assemble bundles';
            AssignIcon( btn2,'ssDNA.jpg' ); btn2.TooltipString = 'Use helical model to specify ssDNA lengths on scaffold connections';
            AssignIcon( btn3,'OrthogonalR.jpg' ); btn3.TooltipString = 'Snap to the closest X-Y-Z axis in Global coordinate';
            AssignIcon( btn9,'AddOrRemove.jpg' );  btn9.TooltipString = 'Add or remove bundles or mechanisms.';
            AssignIcon( btn10,'EditBundle.jpg' );  btn10.TooltipString = 'Edit the cylinder model of the selected bundle.If changed, remember to export again before ssDNA. ';
            AssignIcon( btn7,'LeftClick.jpg' ); btn7.TooltipString = 'Equivalent to left-click the bundle connectivity graph with mouse, assigning connectivity between two bundles. Useful for more-bundle cases ';
            AssignIcon( btn8,'RightClick.jpg' ); btn8.TooltipString = 'Equivalent to right-click the bundle connectivity graph with mouse, canceling the connectivity. Useful for more-bundle cases ';
            AssignIcon( btn5,'SaveT.jpg' ); btn5.TooltipString = 'Save the entire mechanism for future use.';
            AssignIcon( btn6,'LoadT.jpg' ); btn6.TooltipString = 'Load a mechanism.';
            AssignIcon( Tg1,'GraphRep.jpg' ); Tg1.TooltipString = 'Switch to graph reprensentation';
            AssignIcon( Tg2,'Table1.jpg' ); Tg2.TooltipString = 'Switch to table reprensentation with all combinations';
            AssignIcon( Tg3,'Table2.jpg' ); Tg3.TooltipString = 'Switch to table reprensentation only with connected edges';
            
            Tg1.Callback=@(src,evn)SwitchVisualGraph(obj,src,evn,Tg1,Tg2,Tg3,ax2,tH2,tH3 ) ;
            Tg2.Callback=@(src,evn)SwitchVisualGraph(obj,src,evn,Tg1,Tg2,Tg3,ax2,tH2,tH3 ) ;
            Tg3.Callback=@(src,evn)SwitchVisualGraph(obj,src,evn,Tg1,Tg2,Tg3,ax2,tH2,tH3 ) ;
            checkH2_GL.TooltipString = 'Transform the selected bundle by the global XYZ(fixed) or the local xyz coordinate(moving).';
            checkH.TooltipString = 'Switch between translation and rotation mode. Use keyboard "QWEASD"="+xyz -xyz" to control.';
            Tg1.Value=1;
            uistack(ax2,'top'); tH2.Visible='off'; tH3.Visible='off';
            axH= findobj(gcf,'Tag', 'HidenAssemblyMain') ;
            uistack(axH.Children,'top') ;
            evn.IntersectionPoint=[0.3,0.3] ;evn.Button=3;
            LegendBoxing( [],evn,ax ) ;
            
            
%             sdsf=3
            
            %             drawnow;
            %           LineCombination
            %           gctab
            %                t6=toc
        end %end of CreateTF
        
        function IgnoreStrBundleLeaveIndex(obj,src,evn)
            %             evn
            %             src
            Str={'--'};
            if src.UserData.Mode == 0
                for k=1:length(obj.containBundle) ;  Str{k}=strcat('Bundle  ',num2str(k));    end
                src.UserData.Mode=1 ;
            else
                for k=1:length(obj.containBundle) ;  Str{k}=num2str(k);    end
                src.UserData.Mode=0 ;
            end
            src.String= Str;
            
        end
        
        function  SwitchVisualGraph(obj,src,evn,Tg1,Tg2,Tg3,ax2,tH2,tH3 )
            %             [Tg1.Value Tg2.Value Tg3.Value ]
            if Tg1.Value==1
                ax2.Visible='on' ;ax2.Position(3)=0.8; ax2.Position(4)=0.58; tH2.Visible='off' ; tH3.Visible='off' ;
            elseif Tg2.Value==1
                ax2.Visible='off' ; ax2.Position(3:4)=0;       tH2.Visible='on' ; tH3.Visible='off' ;
            elseif Tg3.Value==1
                ax2.Visible='off' ; ax2.Position(3:4)=0;       tH2.Visible='off' ; tH3.Visible='on' ;
            end
            
            %             sdfsf3=34
        end
        
        function QueryConn(obj,src,evn, popupBuni,textConnQS,popupH,patchH)
            Bi=popupBuni.Value;  Bj=src.Value;
            textConnQS.String= strcat('# of Connections = ', num2str( obj.BundleAdjM(Bi,Bj)) );
            
            popupH.Value =Bj;
            SelectBundlePop(obj,popupH,'',patchH) ;
            
        end
        function SetConn(obj,src,evn, popupBuni,popupBunj,textConnQS,patchH)
            Bi=popupBuni.Value;     Bj=popupBunj.Value;
            Edge=union(Bi,Bj);
            LineList= obj.LineCombination;
            if Bi~=Bj
                if ismember(Edge,LineList,'rows')
                    [~,ind]= ismember(Edge,LineList,'rows');
                    evnz.Button=1 ;
                    lineselect(obj,src,evnz,ind,patchH)  ;
                end
            end
        end
        function DeleteConn(obj,src,evn, popupBuni,popupBunj,textConnQS,patchH)
            Bi=popupBuni.Value;     Bj=popupBunj.Value;
            Edge=union(Bi,Bj);
            LineList= obj.LineCombination;
            if Bi~=Bj
                if ismember(Edge,LineList,'rows')
                    [~,ind]= ismember(Edge,LineList,'rows');
                    evnz.Button=2 ;
                    lineselect(obj,src,evnz,ind,patchH)  ;
                end
            end
        end
        
        function SelectBundlePop(obj,src,evn,patchH)
            sBundle=src.Value;
            for k=1:size(patchH,2)
                patchH{1,k}.FaceColor= patchH{1,k}.UserData;
            end
            for k2= 1:length(sBundle)
                patchH{1,sBundle(k2)}.FaceColor=[1,0,0];
            end
            if length(sBundle)==1
                src.UserData.FirstOne = sBundle;
            end
            
        end
        function SaveT(obj,src,evn)
            %             SaveHB=obj;
            S=obj.saveobj;
            ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
            names = fieldnames(ss_Assembly.UserData );
            SavepH2 = cell(size(ss_Assembly.UserData.pH) );  % after 05/09/2019, cadDOM save more
            for k =1:length(SavepH2)
                
                if ~isempty(ss_Assembly.UserData.pH{k})
                    %                     subcell = cell(size(ss_Assembly.UserData.pH{k}) );
                    subXYZ=[-1,-1,-1,-1,-1,-1] ;
                    for k2=1:length(ss_Assembly.UserData.pH{k})
                        if isvalid(ss_Assembly.UserData.pH{k}(k2))
                            %                             subcell{k2}=[ ss_Assembly.UserData.pH{k}(k2).XData' , ss_Assembly.UserData.pH{k}(k2).YData' ,ss_Assembly.UserData.pH{k}(k2).ZData'];
                            subXYZ=[subXYZ; [ss_Assembly.UserData.pH{k}(k2).XData,  ss_Assembly.UserData.pH{k}(k2).YData ,  ss_Assembly.UserData.pH{k}(k2).ZData] ];
                        end
                    end
                    SavepH2{k}=subXYZ;
                end
            end
            
            OHTable= findobj(gcf,'Tag','OHTable') ;
            axOH= findobj(gcf,'Tag','MechOH3D') ;
            SaveOHTable=[];
            if ~isempty(OHTable.Data)
                SaveOHTable.UserData=OHTable.UserData;
                SaveOHTable.Data=OHTable.Data;
                %                 ~isempty(axOH.UserData)
                
                IndConn=  ~cellfun('isempty',axOH.UserData.Conn);   %IndConn=zeros(size(IndConn))==1 ;
                LeaveConn=axOH.UserData.Conn(IndConn) ;
                
                IsvalidCheck = zeros(size(LeaveConn)) ;
                for ci=1:length(IsvalidCheck)
                    if   isvalid(LeaveConn{ci}) ;  IsvalidCheck(ci)=1 ;end
                end
                
                Indtext=  ~cellfun('isempty',axOH.UserData.text) ;  %Indtext=zeros(size(Indtext))==1 ;
                Leavetext=axOH.UserData.text(Indtext) ;
                ConnSaveXYZ=zeros(sum(IsvalidCheck),3,2) ;
                %                 TextSave=[];
                %                 for k = length(Leavetext):-1:1    % Backwards!
                %                     TextSave(k).String = 'xx';
                %                     TextSave(k).Color = [1 0 0];
                %                     TextSave(k).Position = [0 0 0];
                %                 end
                ccLC=1;
                for k=1: length(LeaveConn)
                    %                     [k,LeaveConn{k}.XData]
                    if isvalid(LeaveConn{k})
                        ConnSaveXYZ(ccLC,:,:)= [LeaveConn{k}.XData ; LeaveConn{k}.YData;LeaveConn{k}.ZData] ; ccLC=ccLC+1 ;
                    end
                end
                
                %                 for k=1: length(Leavetext)
                %                     if  isvalid(Leavetext{k})
                %                         TextSave(k).String=Leavetext{k}.String ;
                %                         TextSave(k).Color=Leavetext{k}.Color ;
                %                         TextSave(k).Position=Leavetext{k}.Position ;
                %                     end
                %                 end
                scatterColorCode =axOH.UserData.NicksScatH.CData ;
                SaveOHTable.ConnSaveXYZ2=ConnSaveXYZ;             % SaveOHTable.TextSave=TextSave;
                SaveOHTable.scatterColorCode=scatterColorCode;
            end
            AssemFolder=[ pwd '\Previous Mechanism\'];
            
            if isfield(S ,'ObjFigure')
               if ~isempty(S.ObjFigure)
                S.ObjFigure=[];
               end
            end
                
            uisave({'S','SavepH2','SaveOHTable'} ,[AssemFolder 'Mechanism']);
            %             uisave({'SaveHB','fH'},'SaveHB');            
        end
        
        function SaveT_wTarget(obj,src,evn,Target)
            %             SaveHB=obj;
            S=obj.saveobj;
            ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
            names = fieldnames(ss_Assembly.UserData );
            SavepH2 = cell(size(ss_Assembly.UserData.pH) );  % after 05/09/2019, cadDOM save more
            for k =1:length(SavepH2)
                if ~isempty(ss_Assembly.UserData.pH{k})
                    subXYZ=[-1,-1,-1,-1,-1,-1] ;
                    for k2=1:length(ss_Assembly.UserData.pH{k})
                        if isvalid(ss_Assembly.UserData.pH{k}(k2))
                            subXYZ=[subXYZ; [ss_Assembly.UserData.pH{k}(k2).XData,  ss_Assembly.UserData.pH{k}(k2).YData ,  ss_Assembly.UserData.pH{k}(k2).ZData] ];
                        end
                    end
                    SavepH2{k}=subXYZ;
                end
            end
            S.ObjFigure =[];
            
            OHTable= findobj(gcf,'Tag','OHTable') ;
            axOH= findobj(gcf,'Tag','MechOH3D') ;
            SaveOHTable=[];
            if ~isempty(OHTable.Data)
                SaveOHTable.UserData=OHTable.UserData;
                SaveOHTable.Data=OHTable.Data;
                IndConn=  ~cellfun('isempty',axOH.UserData.Conn);   %IndConn=zeros(size(IndConn))==1 ;
                LeaveConn=axOH.UserData.Conn(IndConn) ;
                IsvalidCheck = zeros(size(LeaveConn)) ;
                for ci=1:length(IsvalidCheck)
                    if   isvalid(LeaveConn{ci}) ;  IsvalidCheck(ci)=1 ;end
                end
                Indtext=  ~cellfun('isempty',axOH.UserData.text) ;  %Indtext=zeros(size(Indtext))==1 ;
                Leavetext=axOH.UserData.text(Indtext) ;
                ConnSaveXYZ=zeros(sum(IsvalidCheck),3,2) ;

                ccLC=1;
                for k=1: length(LeaveConn)
                    %                     [k,LeaveConn{k}.XData]
                    if isvalid(LeaveConn{k})
                        ConnSaveXYZ(ccLC,:,:)= [LeaveConn{k}.XData ; LeaveConn{k}.YData;LeaveConn{k}.ZData] ; ccLC=ccLC+1 ;
                    end
                end
                
                scatterColorCode =axOH.UserData.NicksScatH.CData ;
                SaveOHTable.ConnSaveXYZ2=ConnSaveXYZ;             % SaveOHTable.TextSave=TextSave;
                SaveOHTable.scatterColorCode=scatterColorCode;
            end
            
            
%             AssemFolder=[ pwd '\Previous Mechanism\'];
% %             uisave({'S','SavepH2','SaveOHTable'} ,[AssemFolder 'Mechanism']);
%             save(Target , {'S','SavepH2','SaveOHTable'}) ;
            save(strcat(Target , '.mat') , 'S','SavepH2','SaveOHTable') ;
        end
        
        
        function obj=LoadT(obj,src,evn)
            AssemFolder=[ pwd '\Previous Mechanism\'];
            uiopen([AssemFolder '*.mat']);
            ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
            if exist('S')   %prevent cancel loading
                obj=loadobj(obj,S) ;
                obj.CreateTFWindow(obj.SavePremPair,obj.BundleAdjM);
            end
            
            if exist('SavepH2')   % after 05/09/2019, cadDOM save more
                axes(findobj(gcf,'Tag','AssemblyMain')) ;
                for k=1:length(SavepH2)
                    if ~isempty(SavepH2{k})
                        %                         LL= find((cellfun('length',SavepH{k})~=0) ); XYZ_Together=zeros()
                        %                         for k2=1:length(LL)
                        %                             %                             [k,k2]
                        %                           ss_Assembly.UserData.pH{k}=  plot3( SavepH{k}{LL(k2)}(:,1) ,SavepH{k}{LL(k2)}(:,2) ,SavepH{k}{LL(k2)}(:,3),'LineWidth',2) ;  % recover line
                        %                         end
                        ss_Assembly.UserData.pH{k}=plot3( SavepH2{k}(:,1:2)' ,SavepH2{k}(:,3:4)'  ,SavepH2{k}(:,5:6)' ,'LineWidth',2) ;  % recover line
                        
                        
                        %                         sdfsf=3
                    end
                end
            end
            
            axCurrent=gca;  axOH= findobj(gcf,'Tag','MechOH3D') ; axes(axOH);
            if exist('SaveOHTable')  % after 05/08/2019, cadDOM save more, all data saved before that can't be loaded for OH.
                if ~isempty(SaveOHTable)
                    if isfield(SaveOHTable,'ConnSaveXYZ2' ) % after 09/28/2018, cadDOM save overhang graphics
                        axOH =findobj(gcf,'Tag','MechOH3D') ;
                        OHTable= findobj(gcf,'Tag','OHTable') ; OHTable.Visible='on';
                        OHbtn=findobj(gcf,'Tag','OHfunction') ;
                        
                        OverhangInitial(OHbtn,[]);
                        OHTable.Data=SaveOHTable.Data ;
                        OHTable.UserData=SaveOHTable.UserData ;
                        
                        %                     sdfsd=3;
                        
                        axOH.UserData.NicksScatH.CData=  SaveOHTable.scatterColorCode ;
                        nn=size(SaveOHTable.ConnSaveXYZ2 ,1);
                        SaveLineXYZ=zeros(3*size(SaveOHTable.ConnSaveXYZ2 ,1) ,3); cc=1;
                        
                        for k=1: size(SaveOHTable.ConnSaveXYZ2 ,1)
                            MM=reshape(SaveOHTable.ConnSaveXYZ2(k,:,:),3,2 ) ;
                            %                             axOH.UserData.Conn{k} = line(MM(1,:) ,MM(2,:),MM(3,:));
                            SaveLineXYZ(cc:cc+2 ,:) =[MM' ; nan nan nan ]; cc=cc+3;
                            axOH.UserData.Conn{k} = plot3(MM(1,:) ,MM(2,:),MM(3,:),'LineWidth',2);
                            axOH.UserData.Conn{k}.ButtonDownFcn=@(src,evn)CancelOH(obj,src,evn,axOH) ; % function in hyperbudnle, synchronize with overhangInitial.m
                            %                             PP= SaveOHTable.TextSave(k).Position;
                            %                             axOH.UserData.text{k,1} = text(PP(1) ,PP(2),PP(3),SaveOHTable.TextSave(k).String,'Color',SaveOHTable.TextSave(k).Color );
                            %                             PP2 = SaveOHTable.TextSave(k+nn).Position;
                            %                             axOH.UserData.text{k,2} = text(PP2(1) ,PP2(2),PP2(3),SaveOHTable.TextSave(k+nn).String,'Color',SaveOHTable.TextSave(k+nn).Color );
                        end
                        %                         axOH.UserData.Conn{1} =line(SaveLineXYZ(:,1)' ,SaveLineXYZ(:,2)' ,SaveLineXYZ(:,3)');  % 05072019
                    end
                end
            end
            
            
            
            axes(axCurrent);
            uistack(findobj(gcf,'Tag', 'HidenAssemblyMain'),'top') ;
            
            %             if exist('Sgraphics')   % after 08/31/2018, cadDOM save more
            %                 ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
            %                 SaveFieldsExceptGO={ 'LineCombination','LineOption'} ;
            %                 for k=1 :length(SaveFieldsExceptGO)
            %                     ss_Assembly.UserData.(SaveFieldsExceptGO{k})=   Sgraphics.(SaveFieldsExceptGO{k}) ;
            %                 end
            %             end
            
            
        end
        
        function  CancelOH(obj,src,evn,ax)
            % make sure synchronize with overhangInitial.m
            t=findobj(gcf,'Tag','OHTable');
            BAori=     [ones(1,3); 0.94*ones(1,3)] ;
            IndsValid= cellfun(@isempty,ax.UserData.Conn); IndsValid=find(~IndsValid );
            for k=1:length(IndsValid)
                if isequal(src,ax.UserData.Conn{k} )
                    NthConn=k ;
                end
            end
            if  src.LineWidth ==2 % hightlight to cancel
                if evn.Button==1
                    opts.Interpreter = 'tex';      opts.Default = 'Yes';
                    quest = 'Are you sure to cancel this ??';
                    answer = questdlg(quest,'Cancel this pair of overhangs',...
                        'Yes','No',opts) ;
                    if strcmp(answer,'Yes')
                        ax.UserData.NicksScatH.CData(t.UserData{NthConn,6} , : )   =zeros(2,3) ;
                        delete( ax.UserData.Conn{NthConn}) ;
                        ax.UserData.Conn(NthConn)=[];
                        t.UserData(NthConn,:) =[];
                        t.Data(NthConn,:) =[];
                        t.BackgroundColor(NthConn,:)=[];
                    end
                else
                    if evn.Button==2
                        src.LineWidth =0.3 ;% change to thin
                        t.Data{NthConn,9}=false;
                    end
                end
            else
                if evn.Button==1
                    BaColor2 =  repmat(BAori,k,1); BaColor2=BaColor2(1:k+1,:);
                    
                    t.BackgroundColor=BaColor2;
                    t.BackgroundColor(NthConn,:)=[0.9,0.9,0.6] ;
                elseif evn.Button==3
                    src.LineWidth =2 ;% change to bold, ready to be able to delete.
                    t.Data{NthConn,9}=true;
                end
            end
        end
        
        
        
        
        function checkFcn(obj,src,evn,editH,txtH2)
            if  src.Value==0
                txtH2.String = 'Unit : nm';
                editH.String = '5' ;
            else
                txtH2.String = 'Unit : deg';
                editH.String = '30' ;
            end
        end
        
        function ClickPatch(obj,src,evn, PopBundle,patchH ,popupBuni , popupBunj ,textConnQS ,checkH_Rcenter)
            %             src
%             evn.Button
                        
             %----------           
            BundleInd= 0; BundleInd2=[];
            for k= 1:size(patchH,2)
                if isequal(patchH{1,k},src)
                    BundleInd=k;
                    BundleInd2=union(PopBundle.Value,BundleInd) ;
                end
            end
            if evn.Button==2
                PopBundle.Value=BundleInd2 ;  % multiple selection
            elseif evn.Button==1
                PopBundle.Value=BundleInd ;
            elseif checkH_Rcenter.Value==1 
                if isempty(obj.Rcenter_scatterH)
                    obj.Rcenter_scatterH = scatter3(evn.IntersectionPoint(1),evn.IntersectionPoint(2),evn.IntersectionPoint(3),36,'ob','filled') ;
                else
                   obj.Rcenter_scatterH.XData =  evn.IntersectionPoint(1) ;
                   obj.Rcenter_scatterH.YData =  evn.IntersectionPoint(2) ;
                   obj.Rcenter_scatterH.ZData =  evn.IntersectionPoint(3) ;
                end
            end
            
            SelectBundlePop(obj,PopBundle,[],patchH) ;
            if evn.Button==1
                popupBuni.Value=BundleInd ;
                QueryConn(obj,popupBuni,[], popupBunj,textConnQS,PopBundle,patchH) ;
                PopBundle.UserData.FirstOne = BundleInd;
            elseif evn.Button==3 && checkH_Rcenter.Value==0 
                popupBunj.Value=BundleInd ;
                QueryConn(obj,popupBunj,[], popupBuni,textConnQS,PopBundle,patchH) ;
                PopBundle.UserData.FirstOne = BundleInd;
            end
            
        end
        
        
        function lineselect(obj,src,evn,lineIndex,patchH)
            %            fH= src.Parent.Parent;
            fH= findobj(gcf,'Tag','ss_Assembly') ;
            AdjM=obj.BundleAdjM;
            
            B=triu(AdjM) ;
            C=tril(AdjM) ;
            II= (B)~=transpose(C) ;
%             if sum(sum(II))>0
%                 fsdfsf=1213
%             end
%             
            AdjM=B+transpose(B) ;
            obj.BundleAdjM=AdjM ;
            
            Combination= obj.LineCombination;
            fH.UserData.tb3.Data=fH.UserData.tb2.Data ;
            if  evn.Button==1    %Add connect
                %            strx = inputdlg('Enter Number of Connection:','Sample', [1 30]);
                AssignConn(Combination(lineIndex,:),obj,fH,patchH,lineIndex);
                %                                obj.choice
                if strcmp(obj.choice.UseMethod,'Manual')
                    obj.TakeSeveralV3{lineIndex}=  obj.choice.ManualRes ;
                    obj.TakeSeveralV3Type{lineIndex} ='manual';
                else
                    obj.TakeSeveralV3Type{lineIndex} ='auto';
                end
                %              end
                x=round(obj.choice.NumberOfConnect);   % this specific pair choice
                
                
                obj.LineOption(lineIndex,1)=obj.choice.Bundle1;  % real data, 1->Both, 2->Side , 3-> End  ,4->Manual
                obj.LineOption(lineIndex,2)=obj.choice.Bundle2;
                switch obj.choice.Bundle1   % just for showing, not data
                    case 1
                        fH.UserData.tb2.Data{ lineIndex,2}='Both';       fH.UserData.tb3.Data{ lineIndex,2}='Both';
                    case 2
                        fH.UserData.tb2.Data{ lineIndex,2}='Side only';  fH.UserData.tb3.Data{ lineIndex,2}='Side only';
                    case 3
                        fH.UserData.tb2.Data{ lineIndex,2}='End only';  fH.UserData.tb3.Data{ lineIndex,2}='End only';
                    case 4
                        fH.UserData.tb2.Data{ lineIndex,2}='Manual';  fH.UserData.tb3.Data{ lineIndex,2}='Manual';
                end
                switch obj.choice.Bundle2 % just for showing, not data
                    case 1
                        fH.UserData.tb2.Data{ lineIndex,4}='Both'; fH.UserData.tb3.Data{ lineIndex,4}='Both';
                    case 2
                        fH.UserData.tb2.Data{ lineIndex,4}='Side only';  fH.UserData.tb3.Data{ lineIndex,4}='Side only';
                    case 3
                        fH.UserData.tb2.Data{ lineIndex,4}='End only';  fH.UserData.t3.Data{ lineIndex,4}='End only';
                    case 4
                        fH.UserData.tb2.Data{ lineIndex,4}='Manual'; fH.UserData.tb3.Data{ lineIndex,4}='Manual';
                end
                fH.UserData.tb2.Data{ lineIndex,5} = obj.choice.NumberOfConnect ;   % just for showing, not data
                fH.UserData.tb3.Data{ lineIndex,5} = obj.choice.NumberOfConnect ;   % just for showing, not data
                if   isempty(x)
                    warndlg('Please enter numerical value','!! Warning !!'); % a little make no sense, since using popupmenu to limit choice
                    return;
                end
                if    x<=0  || x>30
                    warndlg('Please Appropriate Integer','!! Warning !!') % a little make no sense, since using popupmenu to limit choice
                    return;
                end
                
                obj.BundleAdjM(Combination(lineIndex,1),Combination(lineIndex,2))=x;  % real data
                obj.BundleAdjM(Combination(lineIndex,2),Combination(lineIndex,1))=x;
                
                %----- graphical showing update
                if strcmp(obj.choice.UseMethod,'Auto')
                    fH.UserData.BundleConnectivity{lineIndex}.Color=[1 0 0] ; %turn into red
                elseif strcmp(obj.choice.UseMethod,'Manual')
                    fH.UserData.BundleConnectivity{lineIndex}.Color=[1 0 1] ; %turn into pink
                end
                fH.UserData.BundleConnectivity{lineIndex}.LineStyle='-';  % graph representation
                fH.UserData.textNumberHandle{lineIndex}.String=num2str(x);
                %-------
                
                textConnQS=findobj(gcf,'Tag','textConnQS') ;  % uitext to show # of connections
                textConnQS.String=strcat('Conn =',num2str(x));
                
                
            else                 %remove connect
                AdjM(Combination(lineIndex,1),Combination(lineIndex,2))=0;
                AdjM(Combination(lineIndex,2),Combination(lineIndex,1))=0;
                AdjM;
                GBundle=graph(AdjM);bins = conncomp(GBundle);
                %                 if max(bins)==1   %All bundle is connected, no along bundle
                obj.TakeSeveralV3{lineIndex}=[];  % clear forcedConnection in this edge
                obj.BundleAdjM=AdjM;
                fH.UserData.BundleConnectivity{lineIndex}.Color=[0 0 0] ; %turn into black
                fH.UserData.BundleConnectivity{lineIndex}.LineStyle=':';
                fH.UserData.textNumberHandle{lineIndex}.String=num2str(0);
                textConnQS=findobj(gcf,'Tag','textConnQS') ;
                textConnQS.String=strcat('Conn =',num2str(0));
                fH.UserData.tb2.Data{ lineIndex,5} = 0;  fH.UserData.tb3.Data{ lineIndex,5} = 0;
                fH.UserData.tb2.Data{ lineIndex,2}='Both'; fH.UserData.tb3.Data{ lineIndex,2}='Both';
                fH.UserData.tb2.Data{ lineIndex,4}='Both'; fH.UserData.tb3.Data{ lineIndex,4}='Both';
                obj.TakeSeveralV3Type{lineIndex} ='auto';
                delete(fH.UserData.pH{lineIndex}) ;
                %                 else
                %                     fprintf('Create isolated bundles. Nothing happen. \n')
                %                 end
            end
            obj.BundleAdjM;
            IndsRemain= zeros(size(fH.UserData.tb2.Data ,1),1) ;
            for k=1:size(fH.UserData.tb2.Data ,1)
                if fH.UserData.tb2.Data{k,5}~=0
                    IndsRemain(k)=1;
                end
            end
            fH.UserData.tb3.Data=fH.UserData.tb2.Data(IndsRemain==1,:) ;
            
        end % end of lineselect
        
        function tH2Change(obj,src,evn,patchH,TheOtherTable)  % if users change the connection setting by table
            %                             obj.LineOption(lineIndex,1)=obj.choice.Bundle1;  % real data, 1->Both, 2->Side , 3-> End  ,4->Manual
            %                 obj.LineOption(lineIndex,2)=obj.choice.Bundle2;
            % % % evn.Indices(1)== lineIndex
            %             lineIndex=evn.Indices(1);
            ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
            [~,lineIndex]=ismember([src.Data{evn.Indices(1),1} ,src.Data{evn.Indices(1),3}] , obj.LineCombination ,'rows') ;
            ss_Assembly.UserData.tb3.Data=ss_Assembly.UserData.tb2.Data ;
            if evn.Indices(2)==5
                ss_Assembly.UserData.tb3.Data{lineIndex,5} =evn.NewData ;            ss_Assembly.UserData.tb2.Data{lineIndex,5} =evn.NewData ;
            elseif evn.Indices(2)==4
                ss_Assembly.UserData.tb3.Data{lineIndex,4} =evn.NewData ;            ss_Assembly.UserData.tb2.Data{lineIndex,4} =evn.NewData ;
            elseif evn.Indices(2)==2
                ss_Assembly.UserData.tb3.Data{lineIndex,2} =evn.NewData ;            ss_Assembly.UserData.tb2.Data{lineIndex,2} =evn.NewData ;
            end
            
            if evn.Indices(2) ~=5   % just number of connections
                [~,b]=ismember( evn.NewData, src.ColumnFormat{evn.Indices(2)} ) ;
                if b<4
                    if evn.Indices(2) ==2   % change Bi type
                        obj.LineOption(lineIndex ,1)=b;
                        if  strcmp(evn.PreviousData,'Manual')
                            src.Data{lineIndex,4} =  src.Data{lineIndex,2} ;
                        end
                    elseif  evn.Indices(2) ==4   % change Bi type
                        obj.LineOption(lineIndex ,2)=b;
                        if  strcmp(evn.PreviousData,'Manual')
                            src.Data{lineIndex,2} =  src.Data{lineIndex,4} ;
                        end
                    end
                    if src.Data{lineIndex,5}==0
                        %----- graphical showing update
                        ss_Assembly.UserData.BundleConnectivity{lineIndex}.Color=[0 0 0] ; %turn into pink
                        ss_Assembly.UserData.BundleConnectivity{lineIndex}.LineStyle='--';  % graph representation
                        ss_Assembly.UserData.textNumberHandle{lineIndex}.String=num2str(0);
                        %-------
                        textConnQS=findobj(gcf,'Tag','textConnQS') ;  % uitext to show # of connections
                        textConnQS.String=strcat('Conn =',num2str(0));
                    elseif src.Data{lineIndex,5}>0
                        %----- graphical showing update
                        ss_Assembly.UserData.BundleConnectivity{lineIndex}.Color=[1 0 0] ; %turn into pink
                        ss_Assembly.UserData.BundleConnectivity{lineIndex}.LineStyle='-';  % graph representation
                        ss_Assembly.UserData.textNumberHandle{lineIndex}.String=num2str(src.Data{lineIndex,5});
                        %-------
                        textConnQS=findobj(gcf,'Tag','textConnQS') ;  % uitext to show # of connections
                        textConnQS.String=strcat('Conn =',num2str(src.Data{lineIndex,5}));
                        
                    end
                    
                elseif b==4 && or(evn.Indices(2) ==2 ,evn.Indices(2) ==4) % use manually assign
                    AssignConn_OnlyManual(obj.LineCombination(lineIndex,:),obj,ss_Assembly,patchH,lineIndex);  % special dialog excluding searching options, only for manual
                    obj.TakeSeveralV3{lineIndex}=  obj.choice.ManualRes ;
                    obj.TakeSeveralV3Type{lineIndex} ='manual';
                    
                    x=round(obj.choice.NumberOfConnect);   % this specific pair choice
                    obj.LineOption(lineIndex,1)=obj.choice.Bundle1;  % real data, 1->Both, 2->Side , 3-> End  ,4->Manual
                    obj.LineOption(lineIndex,2)=obj.choice.Bundle2;
                    ss_Assembly.UserData.tb2.Data{ lineIndex,5} = x ;   % just for showing, not data
                    obj.BundleAdjM(obj.LineCombination(lineIndex,1),obj.LineCombination(lineIndex,2))=x;  % real data
                    obj.BundleAdjM(obj.LineCombination(lineIndex,2),obj.LineCombination(lineIndex,1))=x;
                    src.Data{ lineIndex,2}='Manual';  src.Data{ lineIndex,4}='Manual';
                    %----- graphical showing update
                    ss_Assembly.UserData.BundleConnectivity{lineIndex}.Color=[1 0 1] ; %turn into pink
                    ss_Assembly.UserData.BundleConnectivity{lineIndex}.LineStyle='-';  % graph representation
                    ss_Assembly.UserData.textNumberHandle{lineIndex}.String=num2str(x);
                    %-------
                    textConnQS=findobj(gcf,'Tag','textConnQS') ;  % uitext to show # of connections
                    textConnQS.String=strcat('Conn =',num2str(x));
                end
            else  % assign number of connection
                if isnan(evn.NewData) ||  (evn.NewData)<0  || ( strcmp(src.Data{ lineIndex,4},'Manual') && evn.NewData~=0)  % prevent unreasonable input
                    fw=msgbox('Please enter reasonable number. If opts are Manual, please cancel connection and do manually assign again. !! ');
                    waitfor(fw);
                    src.Data{lineIndex,evn.Indices(2)}=evn.PreviousData ;
                    return ;
                end
                
                if evn.NewData>0  % change # connection to positive numbers
                    ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
                    MatchCase(obj,src.Data);
                    obj.BundleAdjM(obj.LineCombination(lineIndex,1),obj.LineCombination(lineIndex,2))=evn.NewData;  % real data
                    obj.BundleAdjM(obj.LineCombination(lineIndex,2),obj.LineCombination(lineIndex,1))=evn.NewData;  % real data
                    %----- graphical showing update
                    ss_Assembly.UserData.BundleConnectivity{lineIndex}.Color=[1 0 0] ; %turn into pink
                    ss_Assembly.UserData.BundleConnectivity{lineIndex}.LineStyle='-';  % graph representation
                    ss_Assembly.UserData.textNumberHandle{lineIndex}.String=num2str(evn.NewData);
                    %-------
                    textConnQS=findobj(gcf,'Tag','textConnQS') ;  % uitext to show # of connections
                    textConnQS.String=strcat('Conn =',num2str(evn.NewData));
                elseif evn.NewData==0  % change # connection to positive numbers
                    AdjM=obj.BundleAdjM;
                    AdjM(obj.LineCombination(lineIndex,1),obj.LineCombination(lineIndex,2))=0;
                    AdjM(obj.LineCombination(lineIndex,2),obj.LineCombination(lineIndex,1))=0;
                    GBundle=graph(AdjM);bins = conncomp(GBundle);
                    if max(bins)==1   %All bundle is connected, no along bundle
                        obj.TakeSeveralV3{lineIndex}=[];  % clear forcedConnection in this edge
                        obj.BundleAdjM=AdjM;
                        ss_Assembly.UserData.BundleConnectivity{lineIndex}.Color=[0 0 0] ; %turn into black
                        ss_Assembly.UserData.BundleConnectivity{lineIndex}.LineStyle=':';
                        ss_Assembly.UserData.textNumberHandle{lineIndex}.String=num2str(0);
                        textConnQS=findobj(gcf,'Tag','textConnQS') ;
                        textConnQS.String=strcat('Conn =',num2str(0));
                        ss_Assembly.UserData.tb2.Data{ lineIndex,5} = 0;
                        %
                        obj.TakeSeveralV3Type{lineIndex} ='auto';
                        delete(ss_Assembly.UserData.pH{lineIndex}) ;
                    else
                        src.Data{lineIndex,evn.Indices(2)}=evn.PreviousData ;
                        ss_Assembly.UserData.tb2.Data{lineIndex,evn.Indices(2)}=evn.PreviousData ;
                    end
                end
                
            end
            
            %             if isequal(ss_Assembly.UserData.tb2,src)
            IndsRemain= zeros(size(ss_Assembly.UserData.tb2.Data ,1),1) ;
            for k=1:size(ss_Assembly.UserData.tb2.Data ,1)
                if ss_Assembly.UserData.tb2.Data{k,5}~=0
                    IndsRemain(k)=1;
                end
            end
            ss_Assembly.UserData.tb3.Data=src.Data(IndsRemain==1,:) ;
            %                 sychro=1
            %             end
            
        end % end of tH2Change
        
        function MatchCase(obj,TableData)
            Mat= zeros(size(obj.LineOption)) ;
            for k=1: size(Mat,1)
                switch TableData{k,2}
                    case 'Both'
                        Mat(k,1) =1 ;
                    case 'Side only'
                        Mat(k,1) =2 ;
                    case 'End only'
                        Mat(k,1) =3 ;
                    case 'Manual'
                        Mat(k,1) =4 ;
                end
                
                switch TableData{k,4}
                    case 'Both'
                        Mat(k,2) =1 ;
                    case 'Side only'
                        Mat(k,2) =2 ;
                    case 'End only'
                        Mat(k,2) =3 ;
                    case 'Manual'
                        Mat(k,2) =4 ;
                end
            end
            obj.LineOption=Mat ;
        end
        
        function obj=SimulateAssembly(obj,src,evn)
            %             fH=src.Parent  ;% figure handle
            NewConnectList=zeros( sum(sum(obj.BundleAdjM))/2,4) ;  %[ Bi, BiXpIndex, Bj ,BjXpindex]
            nC=1;
            for edgei=1:size( obj.LineCombination,1)
                Bij=obj.LineCombination(edgei,:) ;
                for nConn=1: size( obj.TakeSeveralV3{edgei} ,1)
                    XX=   obj.TakeSeveralV3{edgei};
                    NewConnectList(nC,:)=[ Bij(1),XX(nConn,1) ,  Bij(2),XX(nConn,2) ];
                    nC=nC+1;
                end
            end
            
            %             obj.containBundle{1}.ExternalXoverAsFB;    old          % [cyl, min(base), max(base),X,Y,Z]  or [base, -Cyl1,-Cyl2,X,Y,Z]
            %  obj.containBundle{1}.ExternalXoverAsFB;              % [cyl, Ba1, C2,ba2,X,Y,Z]  or [-cyl1, Ba1, -C2,ba2,X,Y,Z]
            %Convert Xover pos index to two base representaton with 5'& 3'
            ForcedConnecteList=zeros(2*size(NewConnectList,1) ,6);
            for k=1:size(NewConnectList,1)
                Link=NewConnectList(k,:) ;   %[ Bi, BiXpIndex, Bj ,BjXpindex]
                %                 if k==55
                %                     sdfsf=3
                %                 end
                %                 k
                MA=obj.containBundle{Link(1)}.ExternalXoverAsFB(Link(2)  ,1:4);
                MB=obj.containBundle{Link(3)}.ExternalXoverAsFB(Link(4)  ,1:4);
                
                % MA2 | MB2 ==[ 's' ; 'd' ]  follow this form
                if sum(MA<0)~=2    %normal outer xover
                    if xor( ismember(MA(1),obj.containBundle{Link(1)}.AGroup),obj.containBundle{Link(1)}.AGroupGoUp)
                        MA2=[Link(1),MA(1),min(MA([2,4])) ; Link(1),MA(1),max(MA([2,4])) ];
                    else
                        MA2=[Link(1),MA(1),max(MA([2,4])) ; Link(1),MA(1),min(MA([2,4])) ];
                    end
                else              %side xover
                    MA([1,3])=-MA([1,3]);
                    Crit1=~xor( ismember(MA(1),obj.containBundle{Link(1)}.AGroup),obj.containBundle{Link(1)}.AGroupGoUp) ; % Cyl1 MA2 move up
                    if xor(Crit1, MA(2) <obj.containBundle{Link(1)}.Zbase1(MA(1)))
                        MA2=[Link(1),MA(3),MA(4) ; Link(1),MA(1),MA(2)];
                    else
                        MA2=[Link(1),MA(1),MA(2) ; Link(1),MA(3),MA(4)];
                    end
                end
                
                if sum(MB<0)~=2   %normal outer xover
                    if xor( ismember(MB(1),obj.containBundle{Link(3)}.AGroup),obj.containBundle{Link(3)}.AGroupGoUp)
                        MB2=[Link(3),MB(1),min(MB([2,4])) ; Link(3),MB(1),max(MB([2,4])) ];
                    else
                        MB2=[Link(3),MB(1),max(MB([2,4])) ; Link(3),MB(1),min(MB([2,4]))];
                    end
                else              %side xover
                    MB([1,3])=-MB([1,3]);
                    Crit2=~xor( ismember(MB(1),obj.containBundle{Link(3)}.AGroup),obj.containBundle{Link(3)}.AGroupGoUp); % Cyl MB2 move up
                    
                    if  xor(Crit2, MB(2) <obj.containBundle{Link(3)}.Zbase1(MB(1))  )
                        MB2=[Link(3),MB(3),MB(4) ; Link(3),MB(1),MB(2)];
                    else
                        MB2=[Link(3),MB(1),MB(2) ; Link(3),MB(3),MB(4)];
                    end
                end
                QQQ=[MB2(2,:), MA2(1,:) ; MA2(2,:), MB2(1,:)];
                ForcedConnecteList(2*k-1,:)=[MB2(1,:), MA2(2,:)];
                ForcedConnecteList(2*k,:)=[MA2(1,:), MB2(2,:)];
            end
            obj.ForcedConnectList=ForcedConnecteList;
            %--------------
            FCListRelative=zeros(size(ForcedConnecteList));  %more convenient for extract XYZ
            for iFC=1:size(FCListRelative,1)
                VV=ForcedConnecteList(iFC,:);
                FCListRelative(iFC,1:2)= VV(1:2);
                QQA=VV(3)-  obj.containBundle{VV(1)}.Zbase1(VV(2))+1;
                MAxBase=obj.containBundle{VV(1)}.Zbase2(VV(2))-obj.containBundle{VV(1)}.Zbase1(VV(2))+1;
                if QQA>0 && QQA<=MAxBase
                    FCListRelative(iFC,3)= VV(3)-  obj.containBundle{VV(1)}.Zbase1(VV(2))+1;
                elseif QQA<=0   %Bottom end Xover
                    FCListRelative(iFC,3)=1;
                elseif QQA>MAxBase %Top  end Xover
                    FCListRelative(iFC,3)=MAxBase;
                end
                
                FCListRelative(iFC,4:5)= VV(4:5);
                QQB=VV(6)-  obj.containBundle{VV(4)}.Zbase1(VV(5))+1;
                MAxBase2=obj.containBundle{VV(4)}.Zbase2(VV(5))-obj.containBundle{VV(4)}.Zbase1(VV(5))+1;
                
                if QQB>0 && QQB<=MAxBase2
                    FCListRelative(iFC,6)= VV(6)-  obj.containBundle{VV(4)}.Zbase1(VV(5))+1;
                elseif QQB<=0   %Bottom end Xover
                    FCListRelative(iFC,6)=1;
                elseif  QQB>MAxBase2 %Top  end Xover
                    FCListRelative(iFC,6)=MAxBase2;
                end
            end
            obj.FCListRelativeO=FCListRelative;
            
            alpha=-0.1;Amp=0.1;
            alpha2=-0.1;Amp2=0.5;
            %----------------
            f444=figure(444);cla;
            ax1=subplot(2,1,1);
            %              xlabel('Current distance  (nm)') ;ylabel('Attactive Displacement (nm)');
            %              xlabel(ax1,'X axis #1', 'FontSize', 20);
            
            ylabel(ax1,'Population')
            x=0:0.1:10;
            
            y=Amp*(1-exp(alpha.*x));
            plot(x,y)
            xlabel(ax1,'Current distance  (nm)') ;ylabel(ax1,'Attactive Displacement (nm)');
            ax2=subplot(2,1,2);
            %              xlabel('Current Sum( r *distance) (nm2)') ;ylabel('Attactive Angular Displacement  (deg) ');
            x=0:0.1:50;
            
            y=Amp2*(1-exp(alpha2.*x));
            plot(x,y)
            xlabel(ax2,'Current Sum( r *distance) (nm2)') ;ylabel(ax2,'Attactive Angular Displacement  (deg)');
            
            %--------find Gcenter in terms of cylinder and base which is
            SaveGHelix=cell(1,length(obj.containBundle));
            for Bundlei=1:length(obj.containBundle)
                QQWithPM10Bases =obj.containBundle{Bundlei}.HelixXYZG;
                for k=1:length(QQWithPM10Bases)   %excluding +- 10 bases, which is extra for in case.
                    QQWithPM10Bases{k}(1:10,:)=[];
                    QQWithPM10Bases{k}(end-9:end,:)=[];
                end
                SaveGHelix{Bundlei}= QQWithPM10Bases;
            end
            
            f22=figure(22);clf;hold on; axis equal;
            xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
            
            plotH=cell(length(SaveGHelix),20);
            for Bundlej=1:length(SaveGHelix)
                for Cylinderk=1:length(SaveGHelix{Bundlej})
                    PXYZ=SaveGHelix{Bundlej}{Cylinderk};
                    plotH{Bundlej,Cylinderk}= plot3(PXYZ(:,1),PXYZ(:,2),PXYZ(:,3),'.:'); % plot each helix
                end
            end
            
            pHofFC=cell(size(FCListRelative,1),1);
            for FCi=1:size(FCListRelative,1)
                LineAB=FCListRelative(FCi,:);
                indA=LineAB(3);
                PA=SaveGHelix{LineAB(1)}{LineAB(2)}(indA,:) ;
                indB=LineAB(6);
                %                  FCi
                PB=SaveGHelix{LineAB(4)}{LineAB(5)}(indB,:) ;
                PP=[PA;PB];
                pHofFC{FCi}=plot3(PP(:,1),PP(:,2),PP(:,3),'-r','LineWidth',0.1);
                %                  pHofFC{FCi}.Color='none';
            end
            
            %most closest
            PseudoGCenter=zeros(length(SaveGHelix),2); % closest to [Cyl, Base];  %Base Relative representation
            for Bi=1:length(SaveGHelix)
                BigMatrix=zeros(10000,5);kB=1;
                for cylj=1:  length(SaveGHelix{Bi})
                    cc=size(SaveGHelix{Bi}{cylj},1);
                    BigMatrix(kB:kB-1+cc,: )=[SaveGHelix{Bi}{cylj} cylj*ones(cc,1) (1:cc)']  ;
                    kB=kB+size(SaveGHelix{Bi}{cylj},1);
                end
                BigMatrix(sum(BigMatrix,2)==0,:)=[];
                GXYZ= mean(BigMatrix(:,1:3));
                BigMatrix(:,1:3)=BigMatrix(:,1:3)- ones(size(BigMatrix,1) ,1)*GXYZ;
                dVec=BigMatrix(:,1).^2+BigMatrix(:,2).^2+BigMatrix(:,3).^2 ;
                ind=find(dVec==min(dVec)); ind=ind(1);   %----------4/20
                PseudoGCenter(Bi,:)= BigMatrix(ind,4:5);
            end
            
            pHCenter=cell(length(SaveGHelix),1);
            for check=1:length(SaveGHelix)
                XYZ=SaveGHelix{check}{PseudoGCenter(check,1)}(PseudoGCenter(check,2),:);
                pHCenter{check}=scatter3(XYZ(1),XYZ(2),XYZ(3),'O');pHCenter{check}.SizeData=1;
                pHCenter{check}.MarkerFaceColor=[1,0,0];
            end
            
            FCListRelative;
            plotH; % plotH{Bundlej,Cylinderk}= plot3(PXYZ(:,1),PXYZ(:,2),PXYZ(:,3),'.:');
            %              PseudoGCenter
            nloop=1;
            nBundle=length(SaveGHelix);
            %             view([16.4000 30.0000])
            %              pause
            
            tic
            MaxnLoop=2;
            
            ForceConnSD=zeros(size(FCListRelative,1),1);
            for connecti=1:length(ForceConnSD)
                Vec=FCListRelative(connecti,:);
                PLeftConn=SaveGHelix{Vec(1)}{Vec(2)}(Vec(3),:);
                PRightConn=SaveGHelix{Vec(4)}{Vec(5)}(Vec(6),:);
                ForceConnSD(connecti)=norm(PLeftConn-PRightConn);
            end
            SaveTrajSD=zeros(MaxnLoop+1,1);
            SaveTrajSD(1)=norm(ForceConnSD);
            %              F(MaxnLoop) = struct('cdata',[],'colormap',[]);
            %              axis off;
            %              title('');
            %             cc=0;
            while nloop<MaxnLoop
                drawnow
                %                        cc=cc+1;
                %                       F(cc) = getframe(gcf);
                if nloop<20
                    Amp2=Amp2*(nloop+54)/(nloop+55);
                    Amp=Amp*(nloop+81)/(nloop+82);
                end
                dX=zeros(nBundle,3); %3D translation of each iteration, like summation of all force one each bundle
                dM=zeros(nBundle,3); %3D Rotation of each iteration
                for nFC=1:size(FCListRelative)%-------Calcualte small transformation
                    Vec=FCListRelative(nFC,:);
                    PLeft=[ plotH{Vec(1),Vec(2)}.XData(Vec(3)),plotH{Vec(1),Vec(2)}.YData(Vec(3)),plotH{Vec(1),Vec(2)}.ZData(Vec(3))];
                    PRight=[ plotH{Vec(4),Vec(5)}.XData(Vec(6)),plotH{Vec(4),Vec(5)}.YData(Vec(6)),plotH{Vec(4),Vec(5)}.ZData(Vec(6))];
                    distance=norm(PLeft-PRight)-0.9;
                    Ydist=Amp*(1-exp(alpha*distance));
                    VV=(PLeft-PRight)/norm((PLeft-PRight));
                    XX=Ydist*VV*distance/3;   %  internal force
                    %------Left Bundle
                    dX(Vec(1),:)= dX(Vec(1),:) -XX;
                    %-----Right Bundle
                    dX(Vec(4),:)= dX(Vec(4),:) +XX;
                    
                    %-----------------------------Calculate rotaion part
                    %                       %%Right
                    %                     GCcenterLeft=SaveGHelix{Vec(4)}{PseudoGCenter(Vec(4),1)}(PseudoGCenter(Vec(4),2),:);
                    %                     GCcenterRight=SaveGHelix{Vec(1)}{PseudoGCenter(Vec(1),1)}(PseudoGCenter(Vec(1),2),:);
                    
                    GCcenterRight(3)= plotH{Vec(4),PseudoGCenter(Vec(4),1)}.ZData(PseudoGCenter(Vec(4),2));
                    GCcenterRight(2)= plotH{Vec(4),PseudoGCenter(Vec(4),1)}.YData(PseudoGCenter(Vec(4),2)) ;
                    GCcenterRight(1)= plotH{Vec(4),PseudoGCenter(Vec(4),1)}.XData(PseudoGCenter(Vec(4),2)) ;
                    
                    GCcenterLeft(3)= plotH{Vec(1),PseudoGCenter(Vec(1),1)}.ZData(PseudoGCenter(Vec(1),2));
                    GCcenterLeft(2)= plotH{Vec(1),PseudoGCenter(Vec(1),1)}.YData(PseudoGCenter(Vec(1),2)) ;
                    GCcenterLeft(1)= plotH{Vec(1),PseudoGCenter(Vec(1),1)}.XData(PseudoGCenter(Vec(1),2)) ;
                    RVecRight=PRight-GCcenterRight;
                    RVecLeft=PLeft-GCcenterLeft;
                    OneMommentOnRight=cross(RVecRight,XX);
                    OneMommentOnLeft=cross(RVecLeft,-XX);
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
                    dM(Vec(1),:)= dM(Vec(1),:) +SmallRotationLeft;
                    dM(Vec(4),:)= dM(Vec(4),:) +SmallRotationRight;
                end
                
                for Bundlej=1:length(SaveGHelix)
                    for Cylinderk=1:length(SaveGHelix{Bundlej})
                        AllXYZinThisCylinder=[ plotH{Bundlej,Cylinderk}.XData', plotH{Bundlej,Cylinderk}.YData', plotH{Bundlej,Cylinderk}.ZData'];
                        Rx=[1,0,0;0 cosd(dM(Bundlej,1)), sind(dM(Bundlej,1)); 0 -sind(dM(Bundlej,1)), cosd(dM(Bundlej,1))];
                        Ry=[cosd(dM(Bundlej,2)),0,-sind(dM(Bundlej,2)) ;0 ,1 ,0 ; sind(dM(Bundlej,2)), 0 ,cosd(dM(Bundlej,2))];
                        Rz=[cosd(dM(Bundlej,3)),sind(dM(Bundlej,3)),0 ;  -sind(dM(Bundlej,3)) cosd(dM(Bundlej,3))  0 ; 0 ,0 ,1];
                        
                        NewXYZAll=AllXYZinThisCylinder*Rx*Ry*Rz ;
                        
                        plotH{Bundlej,Cylinderk}.XData= NewXYZAll(:,1)' + dX(Bundlej,1);
                        plotH{Bundlej,Cylinderk}.YData= NewXYZAll(:,2)' + dX(Bundlej,2);
                        plotH{Bundlej,Cylinderk}.ZData= NewXYZAll(:,3)' + dX(Bundlej,3);
                    end
                end
                %                  pause(0.1)
                nloop=nloop+1 ;
                ForceConnSD=zeros(size(FCListRelative,1),1);
                for connecti=1:length(ForceConnSD)
                    Vec=FCListRelative(connecti,:);
                    PRight(3)= plotH{Vec(4),Vec(5)}.ZData(Vec(6));
                    PRight(2)= plotH{Vec(4),Vec(5)}.YData(Vec(6)) ;
                    PRight(1)= plotH{Vec(4),Vec(5)}.XData(Vec(6)) ;
                    
                    PLeft(3)= plotH{Vec(1),Vec(2)}.ZData(Vec(3));
                    PLeft(2)= plotH{Vec(1),Vec(2)}.YData(Vec(3)) ;
                    PLeft(1)= plotH{Vec(1),Vec(2)}.XData(Vec(3)) ;
                    ForceConnSD(connecti)=norm(PLeft-PRight);
                    pHofFC{connecti}.XData=[ PLeft(1), PRight(1)];
                    pHofFC{connecti}.YData=[ PLeft(2), PRight(2)];
                    pHofFC{connecti}.ZData=[ PLeft(3), PRight(3)];
                end
                SaveTrajSD(nloop)=norm(ForceConnSD);
                %                      pHofFC
                %                    GCcenterRight(3)= plotH{Vec(4),PseudoGCenter(Vec(4),1)}.ZData(PseudoGCenter(Vec(4),2));
                %                    GCcenterRight(2)= plotH{Vec(4),PseudoGCenter(Vec(4),1)}.YData(PseudoGCenter(Vec(4),2)) ;
                %                    GCcenterRight(1)= plotH{Vec(4),PseudoGCenter(Vec(4),1)}.XData(PseudoGCenter(Vec(4),2)) ;
                
                for check=1:size(PseudoGCenter,1)
                    pHCenter{check}.ZData=plotH{check,PseudoGCenter(check,1)}.ZData(PseudoGCenter(check,2));
                    pHCenter{check}.YData=plotH{check,PseudoGCenter(check,1)}.YData(PseudoGCenter(check,2));
                    pHCenter{check}.XData=plotH{check,PseudoGCenter(check,1)}.XData(PseudoGCenter(check,2));
                end
                title(strcat('Iteration = [',num2str(nloop),' , ', num2str(SaveTrajSD(nloop))  ));
            end  %end of while loop
            f30=figure(30);clf;hold on;
            subplot(2,1,1)
            plot(1:nloop,SaveTrajSD(1:nloop))
            subplot(2,1,2);hold on;
            ForceConnSD=zeros(size(FCListRelative,1),1);
            for connecti=1:length(ForceConnSD)
                Vec=FCListRelative(connecti,:);
                PRight(3)= plotH{Vec(4),Vec(5)}.ZData(Vec(6));
                PRight(2)= plotH{Vec(4),Vec(5)}.YData(Vec(6)) ;
                PRight(1)= plotH{Vec(4),Vec(5)}.XData(Vec(6)) ;
                
                PLeft(3)= plotH{Vec(1),Vec(2)}.ZData(Vec(3));
                PLeft(2)= plotH{Vec(1),Vec(2)}.YData(Vec(3)) ;
                PLeft(1)= plotH{Vec(1),Vec(2)}.XData(Vec(3)) ;
                plot3([PRight(1),PLeft(1)],  [PRight(2),PLeft(2)],[PRight(3),PLeft(3)])
                scatter3(PRight(1),PRight(2),PRight(3),'s');
                scatter3(PLeft(1),PLeft(2),PLeft(3),'d');
            end
            for makevideo=1:1
                %                     v = VideoWriter('AssembleSteward.avi');
                %                     v.FrameRate = 5;
                %                     open(v);
                %                     Ft=F(2:cc);  Lastframe=Ft(end);%the one with fg123
                % %                     Ft(end:end+10)=Ft(end-1);
                %                     Ft(end:end+10)=Lastframe;
                %                     writeVideo(v,Ft)
                %                     close(v)
            end
            %              Amp2
            toc
            for bundi=1:length(obj.containBundle)
                FinalHelixPostionM=[plotH{bundi,1}.XData',plotH{bundi,1}.YData',plotH{bundi,1}.ZData'];
                InitailHelixPostionM=SaveGHelix{bundi}{1};
                %              [regParams,~,~]=absor(FinalHelixPostionM',InitailHelixPostionM');
                [regParams,~,~]=absor(InitailHelixPostionM',FinalHelixPostionM');
                obj.containBundle{bundi}.SimulateTransMFromTM2=regParams.M;
            end
            
            close(f22);
            close(f30);             close(f444);
        end  %end of function simulation
        
        
        function  obj=GetOrthogonal(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ,PremPair )
            clc;
            
            iBundle=popupH.Value;
            
            if length(iBundle)>1 ;return ; end
            
            scatterXYZ=[ patchH{3,iBundle}.XData', patchH{3,iBundle}.YData',patchH{3,iBundle}.ZData' ];
            
            T=obj.containBundle{iBundle}.TransformMatrix2;
            transPart=T(1:3,4)';
            RCurrent=T(1:3,1:3);
            OrthRotationM=cell(24,1);
            List=[];
            for rx= 0:90:270
                Rx=[1,0 , 0; 0 ,cosd(rx), -sind(rx) ;0, sind(rx), cosd(rx)];
                for ry= 0:90:270
                    Ry=[cosd(ry),0 , sind(ry); 0 ,1, 0 ; -sind(ry), 0,cosd(ry)];
                    for rz= 0:90:270
                        Rz=[cosd(rz), -sind(rz), 0 ;sind(rz),cosd(rz),0;0 0 ,1];
                        TotR=Rx*Ry*Rz;
                        VecR=reshape(TotR,[1,9]);
                        List=[List;VecR];
                    end
                end
            end
            List=unique(List,'rows');
            for k=1:size(List,1)
                OrthRotationM{k}=reshape(List(k,:),[3,3]);
            end
            Closeness=zeros(24,1);
            for k=1:24
                MatCompared=OrthRotationM{k};
                Closeness(k)=sum(sum( abs( MatCompared-RCurrent)));
            end
            %              ClosestCase=find();
            Rtarget=OrthRotationM{Closeness==min(Closeness)};
            QQ= patchH{1,iBundle}.Vertices;
            Gcenter=mean(QQ ); transPart2=-Gcenter;
            %---------------
            if nargin ~=7   % not trigger by undo but by button
                TransR =RCurrent/Rtarget ;
                ThisCommand =AssemblyCommand('snapXYZ', TransR,popupH.Value)  ;
                obj.TransformHistroy{end+1}= ThisCommand  ;
            else % undo
                TransR=     inv(obj.TransformHistroy{end}.Trans_Rotate) ;
            end
            
            PlotXYZ=zeros(length( patchH{2,iBundle}),6);
            for gatherXYZ=1:length( patchH{2,iBundle})
                PlotXYZ(gatherXYZ,:)= [patchH{2,iBundle}(gatherXYZ).XData , patchH{2,iBundle}(gatherXYZ).YData,patchH{2,iBundle}(gatherXYZ).ZData];
            end
            PlotXYZV2=[PlotXYZ(:,[1 3 5]) ; PlotXYZ(:,[2 4 6])];
            PlotXYZV2Int=PlotXYZV2;
            
            NewPXYZV2= (PlotXYZV2+ ones(size(PlotXYZV2,1),1)*transPart2)*TransR - ones(size(PlotXYZV2,1),1)*transPart2 ;
            ntop=1:size(NewPXYZV2,1)/2;  nBot=size(NewPXYZV2,1)/2+1:size(NewPXYZV2,1);
            NewPXYZ=[NewPXYZV2(ntop',1) , NewPXYZV2(nBot',1),NewPXYZV2(ntop',2) , NewPXYZV2(nBot',2),NewPXYZV2(ntop',3) , NewPXYZV2(nBot',3)];
            
            for dispatch=1:length( patchH{2,iBundle})
                patchH{2,iBundle}(dispatch).XData=  NewPXYZ(dispatch,1:2);
                patchH{2,iBundle}(dispatch).YData=  NewPXYZ(dispatch,3:4);
                patchH{2,iBundle}(dispatch).ZData=  NewPXYZ(dispatch,5:6);
            end
            NewScaXYZ= (scatterXYZ+ ones(size(scatterXYZ,1),1)*transPart2)*TransR - ones(size(scatterXYZ,1),1)*transPart2 ;
            patchH{3,iBundle}.XData=NewScaXYZ(:,1)' ; patchH{3,iBundle}.YData=NewScaXYZ(:,2)' ; patchH{3,iBundle}.ZData=NewScaXYZ(:,3)' ;
            
            VerticesP1=patchH{1,iBundle}.Vertices;
            NewVerticesP1= (VerticesP1+ ones(size(VerticesP1,1),1)*transPart2)*TransR - ones(size(VerticesP1,1),1)*transPart2 ;
            patchH{1,iBundle}.Vertices=NewVerticesP1;
            %                 OP6=V0CylinderXYZ{iBundle};   %old version , always
            %                 compar to the very beginning
            %                 OP3=[[OP6(:,1);OP6(:,2)], [OP6(:,3);OP6(:,4)], [OP6(:,5);OP6(:,6)]] ;
            OP3=PlotXYZV2Int;
            [regParams,~,~]=absor(OP3',NewPXYZV2');
            tH.Data=regParams.M;
            obj.containBundle{iBundle}.TransformMatrix2=[regParams.R, regParams.t;0,0,0,1]*obj.containBundle{iBundle}.TransformMatrix2 ;
            
            
        end
        %          function getForcedConnecteList(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ,PremPair)
        %
        %          end
        function getForcedConnecteList(obj)
            fH=findobj(gcf,'Tag','ss_Assembly') ; % sub tab handle
            NewConnectList=zeros( sum(sum(obj.BundleAdjM))/2,4) ;  %[ Bi, BiXpIndex, Bj ,BjXpindex]
            nC=1;
            for edgei=1:size( obj.LineCombination,1)
                Bij=obj.LineCombination(edgei,:) ;
                for nConn=1: size( obj.TakeSeveralV3{edgei} ,1)
                    XX=   obj.TakeSeveralV3{edgei};
                    NewConnectList(nC,:)=[ Bij(1),XX(nConn,1) ,  Bij(2),XX(nConn,2) ];
                    nC=nC+1;
                end
            end
            ForcedConnecteList=zeros(2*size(NewConnectList,1) ,6);
            for k=1:size(NewConnectList,1)
                Link=NewConnectList(k,:) ;   %[ Bi, BiXpIndex, Bj ,BjXpindex]
                
                MA=obj.containBundle{Link(1)}.ExternalXoverAsFB(Link(2)  ,1:4);
                MB=obj.containBundle{Link(3)}.ExternalXoverAsFB(Link(4)  ,1:4);
                
                % MA2 | MB2 ==[ 's' ; 'd' ]  follow this form
                if sum(MA<0)~=2    %normal outer xover
                    if xor( ismember(MA(1),obj.containBundle{Link(1)}.AGroup),obj.containBundle{Link(1)}.AGroupGoUp)
                        MA2=[Link(1),MA(1),min(MA([2,4])) ; Link(1),MA(1),max(MA([2,4])) ];
                    else
                        MA2=[Link(1),MA(1),max(MA([2,4])) ; Link(1),MA(1),min(MA([2,4])) ];
                    end
                else              %side xover
                    MA([1,3])=-MA([1,3]);
                    Crit1=~xor( ismember(MA(1),obj.containBundle{Link(1)}.AGroup),obj.containBundle{Link(1)}.AGroupGoUp) ; % Cyl1 MA2 move up
                    if xor(Crit1, MA(2) <obj.containBundle{Link(1)}.Zbase1(MA(1)))
                        MA2=[Link(1),MA(3),MA(4) ; Link(1),MA(1),MA(2)];
                    else
                        MA2=[Link(1),MA(1),MA(2) ; Link(1),MA(3),MA(4)];
                    end
                    
                end
                
                if sum(MB<0)~=2   %normal outer xover
                    if xor( ismember(MB(1),obj.containBundle{Link(3)}.AGroup),obj.containBundle{Link(3)}.AGroupGoUp)
                        MB2=[Link(3),MB(1),min(MB([2,4])) ; Link(3),MB(1),max(MB([2,4])) ];
                    else
                        MB2=[Link(3),MB(1),max(MB([2,4])) ; Link(3),MB(1),min(MB([2,4]))];
                    end
                else              %side xover
                    MB([1,3])=-MB([1,3]);
                    Crit2=~xor( ismember(MB(1),obj.containBundle{Link(3)}.AGroup),obj.containBundle{Link(3)}.AGroupGoUp); % Cyl MB2 move up
                    
                    if  xor(Crit2, MB(2) <obj.containBundle{Link(3)}.Zbase1(MB(1))  )
                        MB2=[Link(3),MB(3),MB(4) ; Link(3),MB(1),MB(2)];
                    else
                        MB2=[Link(3),MB(1),MB(2) ; Link(3),MB(3),MB(4)];
                    end
                    
                end
                ForcedConnecteList(2*k-1,:)=[MB2(1,:), MA2(2,:)];
                ForcedConnecteList(2*k,:)=[MA2(1,:), MB2(2,:)];
                
            end
            obj.ForcedConnectList=ForcedConnecteList;
        end
        
        
        function  obj=ConnectFC(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ,PremPair)
            fH=findobj(gcf,'Tag','ss_Assembly') ; % sub tab handle
            assemblyAxs=findobj(gcf,'Tag','AssemblyMain') ;
            
            obj.getForcedConnecteList;
            obj=obj.findRT;
            obj.ssOption ;
            
            
            
            if ~isempty(obj.ssOption)
                obj.ssOption =GUIAsignSSInfo3( obj,fH, obj.ssOption ,assemblyAxs.View );   % call ssDNA asking Info GUI
            else
                obj.ssOption =GUIAsignSSInfo3( obj,fH,[] ,assemblyAxs.View);   % call ssDNA asking Info GUI
            end
            
            
            %             obj.findCycleList([],PremPair,2);   %get property ScafRouting-------
            %
            %             ConvertScafG(obj);
            %             obj=ConvertScafSQ(obj);      % Get properties:ScafdigitSQ
            %             obj=ConvertScafHC(obj);      % Get properties:ScafdigitHC
            %
            %             fprintf(' found scaffold routing \n')
            %             plotScaf;
            % %             return
            %
            %             type=2;
            %             obj=obj.FindStap(type );    % Get properties: StapList ......., stapBP ,stapBPinCell
            %             obj=obj.FindStapStep2;           %Get property:StapList2
            %             obj=obj.ApartStaples;     %Get property:StapList3
            % %             obj.StapList3=obj.StapList2;
            %
            %             obj=obj.ConvertStap('Square');    % Get properties: DigitStapSQ,   HeadOfStep
            %             obj=obj.ConvertStap('Honeycomb');    % Get properties: DigitStapHC
            %
            %
            %             obj=obj.ExportJSON;
            %              scriptOfPlot
            
            %-------------------------end of main process
        end  % end of function ConnectFC
        
        
        
        function  SortDist(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ)
            %             helloSortDistD=1;
            fH=findobj(gcf,'Tag','ss_Assembly') ;  % figure handle
            aH=findobj(fH.Children,'Tag','AssemblyMain');  %axe handle
            nBundle=size(patchH,2);
            
            GlobalXYZofPossibleXover=cell(1,nBundle);
            for iBun=1:nBundle
                Mat=zeros(length(patchH{3,iBun}.XData),3);  %preallocate for XYZ data
                for jdataPoint=1:size(Mat,1)
                    Mat(jdataPoint,:)=[ patchH{3,iBun}.XData(jdataPoint),patchH{3,iBun}.YData(jdataPoint),patchH{3,iBun}.ZData(jdataPoint) ];
                end
                GlobalXYZofPossibleXover{iBun}=Mat;
            end
            %---------------- ^^ pre-process
            BundleAdj= obj.BundleAdjM   ;
            Preallocateisze=0;
            NumberOfFB=zeros(size(obj.LineCombination,1),1);
            for BunEdgeCheck=1:size(obj.LineCombination,1)
                BundleA=obj.LineCombination(BunEdgeCheck,1);
                BundleB=obj.LineCombination(BunEdgeCheck,2);
                if        BundleAdj(BundleA,BundleB)>0  %need to calculate all combination between BunA and BunB
                    Preallocateisze=Preallocateisze+ size(GlobalXYZofPossibleXover{BundleA} ,1)*size(GlobalXYZofPossibleXover{BundleB} ,1);
                end
                NumberOfFB(BunEdgeCheck)=BundleAdj(BundleA,BundleB);
            end
            tic
            AllCombinationDistTableV2=zeros(Preallocateisze,5); k=1; %[dist, Bundlei, Xoveri,Bundlej, Xoverj]
            for BunEdgeCheck=1:size(obj.LineCombination,1)
                BundleA=obj.LineCombination(BunEdgeCheck,1);
                BundleB=obj.LineCombination(BunEdgeCheck,2);
                if  BundleAdj(BundleA,BundleB)>0  %need to calculate all combination between BunA and BunB
                    for iPoint=1:size(GlobalXYZofPossibleXover{BundleA},1)
                        for jPoint=1:size(GlobalXYZofPossibleXover{BundleB},1)
                            dd=norm([GlobalXYZofPossibleXover{BundleA}(iPoint,:)-GlobalXYZofPossibleXover{BundleB}(jPoint,:)]);
                            AllCombinationDistTableV2(k,:)=[dd,BundleA,iPoint,BundleB,jPoint];
                            k=k+1;
                        end
                    end
                end
            end
            toc
            %--------------
            SortRes=sortrows(AllCombinationDistTableV2,1); %[dist, Bundlei, Xoverindex, Bundlej,Xoverindex]
            
            %              TakeSeveralV3=cell(size(obj.LineCombination,1),1);
            
            
            %               TakeSeveralV3L=cell(size(NumberOfFB));   % manual result put here----------
            TakeSeveralV3L= obj.TakeSeveralV3;    % load all, This is the index of FB between bundles. If use different Tol, it will lead to wrong index.
            AlreadyUseInManually=[-1,-1] ;
            for k=1:size(obj.LineCombination,1)  % clear "auto" line, leave it able to search later
                if strcmp(obj.TakeSeveralV3Type(k) ,"auto")
                    TakeSeveralV3L{k}=[];
                else
                    AlreadyUseInManually=union(AlreadyUseInManually, [obj.LineCombination(k,1)*ones(size( TakeSeveralV3L{k}(:,1))) , TakeSeveralV3L{k}(:,1)],'rows') ;
                    AlreadyUseInManually=union(AlreadyUseInManually, [obj.LineCombination(k,2)*ones(size( TakeSeveralV3L{k}(:,2))) , TakeSeveralV3L{k}(:,2)],'rows') ;
                    
                end
            end
            
            UsedNodeinAll=[-1,-1,-1];   % [Bundleb, Bundlea, XoverBPb] & [Bundlea, Bundleb, XoverBPa]
            %             UsedNode=[-1,-1];
            UsedNode=AlreadyUseInManually ;
            searchRow=1;
            FClinkrecord=[-1,-1,-1,-1];
            while searchRow < size(SortRes,1)-1
                AsOneConsider=SortRes(searchRow,:);
                Cons1=[AsOneConsider(2),AsOneConsider(4),AsOneConsider(3)]; Cons1S=[AsOneConsider(2),AsOneConsider(3)];
                Cons2=[AsOneConsider(4),AsOneConsider(2),AsOneConsider(5)]; Cons2S=[AsOneConsider(4),AsOneConsider(5)];
                
                %                 if ~ismember(Cons1,UsedNodeinAll,'rows') && ~ismember(Cons2,UsedNodeinAll,'rows')
                if ~ismember(Cons1S,UsedNode,'rows') && ~ismember(Cons2S,UsedNode,'rows')
                    [~,indE2]=ismember(AsOneConsider([2,4]),obj.LineCombination,'rows');
                    
                    if  size(TakeSeveralV3L{indE2},1)< NumberOfFB(indE2)  % check if has room
                        switch obj.LineOption(indE2,1)   % real data for searching, no in the table
                            case 1
                                Avail1=(1:size(obj.containBundle{AsOneConsider(2)}.ExternalXoverAsFB,1))';
                            case 2
                                Avail1=find( obj.containBundle{AsOneConsider(2)}.ExternalXoverAsFB (:,1)>0) ;
                            case 3
                                Avail1=find( obj.containBundle{AsOneConsider(2)}.ExternalXoverAsFB (:,1)<0) ;
                        end
                        
                        switch obj.LineOption(indE2,2)
                            case 1
                                Avail2=(1:size(obj.containBundle{AsOneConsider(4)}.ExternalXoverAsFB,1))';
                            case 2
                                Avail2=find( obj.containBundle{AsOneConsider(4)}.ExternalXoverAsFB (:,1)>0) ;
                            case 3
                                Avail2=find( obj.containBundle{AsOneConsider(4)}.ExternalXoverAsFB (:,1)<0) ;
                        end
                        
                        if ismember(AsOneConsider(3),Avail1) && ismember(AsOneConsider(5),Avail2)
                            
                            TakeSeveralV3L{indE2}=[ TakeSeveralV3L{indE2} ; AsOneConsider([3,5])];
                            UsedNodeinAll=[UsedNodeinAll;Cons1; Cons2];
                            UsedNode=[UsedNode;AsOneConsider(2),AsOneConsider(3);AsOneConsider(4),AsOneConsider(5)];
                            FClinkrecord=[FClinkrecord;AsOneConsider(2:5)];
                        end
                        
                    end
                end
                CurrentNumFB=zeros(size(NumberOfFB));  %update
                for updatei=1:length(CurrentNumFB)
                    CurrentNumFB(  updatei)= size( TakeSeveralV3L{updatei},1) ;
                end
                
                if sum(CurrentNumFB==NumberOfFB)==length(CurrentNumFB)
                    break;
                end
                searchRow=searchRow+1 ;
            end   %end of while loop
            %            FClinkrecord
            FClinkrecord=setdiff(FClinkrecord,[-1,-1,-1,-1],'rows');  CompJoint.FClinkrecord=FClinkrecord;
            involveFC=unique([FClinkrecord(:,1:2) ;FClinkrecord(:,3:4)],'rows'); CompJoint.involveFC=involveFC;
            NumberofSharing=size(involveFC,1)-size(FClinkrecord,1);
            IndexRep=zeros(size(FClinkrecord,1),2);
            for kFC=1:size(FClinkrecord,1)
                [~,  IndexRep(kFC,1)]= ismember(FClinkrecord(kFC,1:2),involveFC,'rows');
                [~,  IndexRep(kFC,2)]= ismember(FClinkrecord(kFC,3:4),involveFC,'rows') ;
            end
            GraphFC=graph(IndexRep(:,1),IndexRep(:,2) ); CompJoint.GraphFC=GraphFC;
            bins = conncomp(GraphFC); CompJoint.bins=bins;
            obj.CompJointInfo=CompJoint;
            
            
            for eJ=1:size(TakeSeveralV3L,1)
                if ~isempty(fH.UserData.pH{eJ})
                    delete(fH.UserData.pH{eJ})
                end
            end
            TakeSeveralV3L;
            for edgeJ=1:length(TakeSeveralV3L)
                if ~isempty(TakeSeveralV3L{edgeJ})
                    
                    plotXYZ=zeros(NumberOfFB(edgeJ),6);
                    for putin=1:NumberOfFB(edgeJ)
                        BundleA=obj.LineCombination(edgeJ,1);
                        BundleB=obj.LineCombination(edgeJ,2);
                        % try
                        P1=GlobalXYZofPossibleXover{BundleA}(TakeSeveralV3L{edgeJ}(putin,1),:);
                        % catch
                        %     sdf=3
                        % end
                        P2=GlobalXYZofPossibleXover{BundleB}(TakeSeveralV3L{edgeJ}(putin,2),:)   ;
                        plotXYZ(putin,:)=[P1(1),P2(1) , P1(2),P2(2) ,P1(3),P2(3) ];
                    end
                    if ~isempty(fH.UserData.pH{edgeJ})
                        delete(fH.UserData.pH{edgeJ})
                    end
                    fH.UserData.pH{edgeJ}=plot3(aH,plotXYZ(:,1:2)' ,plotXYZ(:,3:4)' ,plotXYZ(:,5:6)','LineWidth',2);
                    obj.TakeSeveralV3=TakeSeveralV3L;
                end
            end
            obj.TakeSeveralV3=TakeSeveralV3L;
            obj.ssOption=[];
            %             TakeSeveralV3L
            fprintf('finish searching connection. \n ')
        end   %end of function SortDist/FindXover
        
        function TotalL = EstimateDsScafL(obj)
            TotalL=0;
            for k=1: length(obj.containBundle )
                dsInThis =  sum(obj.containBundle{k}.Zbase2-obj.containBundle{k}.Zbase1+1)  ;
                TotalL=TotalL+dsInThis;
            end
            fprintf('The estimate scaffold lengths = %s \n', num2str(TotalL));
        end
        
        function TotalL = EstimateDsScafL_part(obj,defaultBundles)
            
            prompt = {'Enter Bundle indexes to calculate approximate scaffold size:'};
            dlgtitle = 'Input';
            dims = [1 35];
            definput = {num2str(defaultBundles)};
            answer = inputdlg(prompt,dlgtitle,dims,definput) ;
            
            QuerryBundles = str2num(answer{1}) ;
            TotalL=0;
            for k=1: length(QuerryBundles )
                dsInThis =  sum(obj.containBundle{QuerryBundles(k)}.Zbase2-obj.containBundle{QuerryBundles(k)}.Zbase1+1)  ;
                TotalL=TotalL+dsInThis;
            end
            fprintf('The estimate scaffold lengths for bundles [%s] = %s \n',answer{1} ,num2str(TotalL));
        end
        
        
        
        
        function obj=keyMove(obj,src,evn,patchH,popupH,checkH,tH,V0CylinderXYZ,editH,textBundle ,checkH2_GL,txtH2,checkH_Rcenter )
%                         df=111
            %             evn
            switch evn.Key
                case 'rightarrow'  % for cadnano tab, mapping 2D base with 3D position
                    sH3D=findobj(gcf,'Tag','sH3D')  ;
                    sH=findobj(gcf,'Tag','sH')  ;
                    if isempty(sH);return;end
                    stp=sH.UserData.Highlight.UserData.Step ;
                    NewInd = sH.UserData.ind +stp;
                    if NewInd<=length(sH.XData)
                        NewXY =[sH.XData(NewInd) , sH.YData(NewInd)  ] ;
                        evnSH.IntersectionPoint=NewXY;
                        HighLightBase(sH,evnSH,sH3D) ;
                    end
                    return
                case 'leftarrow'
                    sH3D=findobj(gcf,'Tag','sH3D')  ;
                    sH=findobj(gcf,'Tag','sH')  ;
                    if isempty(sH);return;end
                    stp=sH.UserData.Highlight.UserData.Step ;
                    NewInd = sH.UserData.ind -stp;
                    if NewInd>0
                        NewXY =[sH.XData(NewInd) , sH.YData(NewInd)  ] ;
                        evnSH.IntersectionPoint=NewXY;
                        HighLightBase(sH,evnSH,sH3D) ;
                    end
                    return
                    
                case 'h'
                    implay('AsssemblyTest.mp4')
                case 'l'
                    obj.EstimateDsScafL ;    return;
                case 'k'
                    defaultBundles= popupH.Value ;
                    obj.EstimateDsScafL_part(defaultBundles) ;    return;
                    
                    %                 case 'm'
                    %                     helpCADDOM(2) ;          return;
                case 'c'
                    fprintf('printing cylinder models for Chimera.\n')
                    file_name='CylinderModel' ;   Radius =1.26;
                    fileID = fopen([pwd filesep file_name '.bild'],'w');
                    for Bi = 1:size(patchH ,2)
                        fprintf(fileID ,'\n' );
                        %                         fprintf(fileID ,'.transparency 0.5\n' );
                        fprintf(fileID ,'.comment (If need colors for bundle %i) .color  r g b\n',Bi );
                        for Cj= 1:length(patchH{2,Bi})
                            XYZ= [patchH{2,Bi}(Cj).XData ;patchH{2,Bi}(Cj).YData ;patchH{2,Bi}(Cj).ZData] ;
                            fprintf(fileID , '.cylinder %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n',XYZ(:,1)',XYZ(:,2)',Radius )    ;
                        end
                    end
                    fclose(fileID);
                    return;
                    
                case 'u'
                    if ~isempty(obj.TransformHistroy)
                        LastCommend=  obj.TransformHistroy{end} ;
                        if strcmp(LastCommend.type,'XYZ_input')
                            checkH.Value=LastCommend.Trans_Rotate ;
                            checkFcn(obj,checkH,[],editH,txtH2)
                            checkH2_GL.Value=LastCommend.Global_Local ;
                            editH.String=num2str( LastCommend.Increment );
                            popupH.Value=LastCommend.Bundle  ;
                            if length(LastCommend.Bundle)>1
                                popupH.UserData.FirstOne= LastCommend.Bundle(1) ;
                            end
                            evn2.Character= LastCommend.Direction ;
                            evn2.Key ='?';
                            keyMove(obj,src,evn2,patchH,popupH,checkH,tH,V0CylinderXYZ,editH,textBundle ,checkH2_GL,txtH2,checkH_Rcenter) ;
                            obj.TransformHistroy(end)=[];
                        elseif strcmp(LastCommend.type,'snapXYZ')
                            popupH.Value=LastCommend.Bundle ;
                            
                            GetOrthogonal(obj,src,evn,patchH,popupH,checkH,tH);
                            obj.TransformHistroy(end)=[];
                            
                        end
                        SelectBundlePop(obj,popupH,'',patchH) ;
                    else
                        fprintf('No transform command history or have all used up Undos \n' );
                    end
                    
                case 'i'
                    ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
                    ReOrderBundles(obj,ss_Assembly )
                    
                    
                    %                 ThisCommand =AssemblyCommand(checkH.Value, checkH2_GL.Value, RotaionIncrement, evn.Character )   ;
                    %                 case 'z'
                    %                     ggtab =gctab  ; %
                    %                     ss_json=findall(gcf,'Tag','ss_json') ;
                    %                     if isequal(ggtab,ss_json)             % curent tab on cadnano
                    %                         sH=findobj(ss_json,'Tag','sH') ;
                    %                         if ~isempty(sH)   % have initial the routing
                    %                             if isfield(sH.UserData,'ScafBCB') % have create red dot in cadnano tab
                    %                                 if sH.UserData.SFC(sH.UserData.ind,1)==sH.UserData.SFC(sH.UserData.ind-1,1)  % select (kth,k+1) base, avoid routing on single base
                    %                                 XY2d= [ sH.XData( sH.UserData.ind) , sH.YData( sH.UserData.ind) ] ;
                    %
                    %                                 sH.UserData.Highlight.XData(2)=XY2d(1) ;
                    %                                 sH.UserData.Highlight.YData(2)=XY2d(2) ;
                    %                                 sH.UserData.Highlight.CData(1:2,1:3) =[1,0,0 ;  0 , 1 , 0.3];
                    %                                 sH.UserData.ScafBCBz=  sH.UserData.ScafBCB ;
                    %                                 %                                 sH.UserData
                    %                                 %                                 fprintf(' key z \n')
                    %                                 else
                    %                                   f = msgbox('Not valid base.');
                    %                                 end
                    %                             end
                    %                         end
                    %                     end
                    %                 case 'x'
                    %                     ggtab =gctab  ; %
                    %                     ss_json=findall(gcf,'Tag','ss_json') ;
                    %                     if isequal(ggtab,ss_json)             % curent tab on cadnano
                    %                         sH=findobj(ss_json,'Tag','sH') ;
                    %                         if ~isempty(sH)   % have initial the routing
                    %                             if isfield(sH.UserData,'ScafBCB') % have create red dot in cadnano tab
                    %                                 if sH.UserData.SFC(sH.UserData.ind,1)==sH.UserData.SFC(sH.UserData.ind-1,1)  % select (kth,k+1) base, avoid routing on single base
                    %                                     XY2d= [ sH.XData( sH.UserData.ind) , sH.YData( sH.UserData.ind) ] ;
                    %                                     if length(sH.UserData.Highlight.XData) ~= 1 % make sure use z key before x key
                    %                                         sH.UserData.Highlight.XData(3)=XY2d(1) ;
                    %                                         sH.UserData.Highlight.YData(3)=XY2d(2) ;
                    %                                         sH.UserData.Highlight.CData =[1,0,0;  0 , 1 , 0.3 ;  0 , 1 , 0.8];
                    %                                         sH.UserData.ScafBCBx=  sH.UserData.ScafBCB ;
                    %                                     end
                    %                                 else                                    f = msgbox('Not valid base.');
                    %
                    %                                 end
                    %                                 %                                 fprintf(' key x \n')
                    %                             end
                    %                         end
                    %                     end
                    %                 case 'v'
                    %                     ggtab =gctab  ; %
                    %                     ss_json=findall(gcf,'Tag','ss_json') ;
                    %
                    %                     if isequal(ggtab,ss_json)             % curent tab on cadnano
                    %                         sH=findobj(ss_json,'Tag','sH') ;
                    %                         if ~isempty(sH)   % have initial the routing
                    %                             if isfield(sH.UserData,'ScafBCB') % have create red dot in cadnano tab
                    %                                 XY2d= [ sH.XData( sH.UserData.ind) , sH.YData( sH.UserData.ind) ] ;
                    %                                 if length(sH.UserData.Highlight.XData) == 3 % make sure have used 'z' and 'x'
                    %                                     if ~isfield(sH.UserData,'ExtraForcedScafXover')
                    %                                         axJson2D=findobj(gcf,'Tag','json2D') ; axes(axJson2D) ;
                    %
                    %                                         MM =[sH.UserData.ScafBCBz ;  sH.UserData.ScafBCBx  ] ;
                    %                                         MM= [ MM(2,4:6), MM(1,1:3) ; MM(1,4:6), MM(2,1:3) ];
                    % %                                                     MM= [ MM(1,1:3), MM(2,4:6) ; MM(2,1:3), MM(1,4:6) ];
                    %                                         sH.UserData.ExtraForcedScafXover = MM ;
                    %                                     else
                    %                                         axJson2D=findobj(gcf,'Tag','json2D') ; axes(axJson2D) ;
                    %                                         MM =[sH.UserData.ScafBCBz ;  sH.UserData.ScafBCBx  ] ;
                    %                                         MM= [ MM(2,4:6), MM(1,1:3) ; MM(1,4:6), MM(2,1:3) ];
                    % %                                                     MM= [ MM(1,1:3), MM(2,4:6) ; MM(2,1:3), MM(1,4:6) ];
                    %                                         sH.UserData.ExtraForcedScafXover = [sH.UserData.ExtraForcedScafXover ; MM] ;
                    %                                     end
                    %                                     plot(sH.UserData.Highlight.XData(2:3),sH.UserData.Highlight.YData(2:3),'-o','Color','k' ,'MarkerFaceColor',ones(1,3)) ;
                    %                                     sH.UserData.Highlight.XData=sH.UserData.Highlight.XData(1) ;
                    %                                     sH.UserData.Highlight.YData=sH.UserData.Highlight.YData(1) ;
                    %                                     sH.UserData.Highlight.CData=[1,0 ,0];
                    %                                     sH.UserData = rmfield(sH.UserData,'ScafBCBz') ;
                    %                                     sH.UserData = rmfield(sH.UserData,'ScafBCBx') ;
                    %
                    %                                 end
                    %                             end
                    %                         end
                    %                     end
                    
            end
            %             gctab
            iBundle=popupH.Value;
            if length(iBundle)>1
                iBundle= [ popupH.UserData.FirstOne, setdiff(iBundle,popupH.UserData.FirstOne)]  ;
            end
            
            
            %-------gather data
            nCyl2=0 ;
            for k =1: length(iBundle)
                nCyl2=nCyl2+length( patchH{2,iBundle(k)}) ;
            end
            
            PlotXYZ=zeros(nCyl2,6);   %[x1 x2 y1 y2 z1 z2]
            PlotXYZV2Ini=cell(length(iBundle),1) ; % for absor use
            k=1; Indp2 = 0;
            for bbi= 1 :length(iBundle)
                ki=k ;
                for gatherXYZ=1:length( patchH{2,iBundle(bbi)})
                    PlotXYZ(k,:)= [patchH{2,iBundle(bbi)}(gatherXYZ).XData , patchH{2,iBundle(bbi)}(gatherXYZ).YData,patchH{2,iBundle(bbi)}(gatherXYZ).ZData];
                    k=k +1 ;
                end
                Indp2=[ Indp2 ,repelem(bbi , gatherXYZ) ];
                PlotXYZV2Ini{bbi} =[PlotXYZ(ki:k-1,[1 3 5]) ; PlotXYZ(ki:k-1,[2 4 6])] ;
            end
            Indp2=Indp2(2:end)' ;
            PlotXYZV2=[PlotXYZ(:,[1 3 5]) ; PlotXYZ(:,[2 4 6])];
            
            scatterXYZ=[nan,nan,nan] ; Indp3 = 0;
            for bbi= 1 :length(iBundle)
                scatterXYZ=[ scatterXYZ ;patchH{3,iBundle(bbi)}.XData', patchH{3,iBundle(bbi)}.YData',patchH{3,iBundle(bbi)}.ZData' ];
                Indp3=[ Indp3 ,repelem(bbi , length(patchH{3,iBundle(bbi)}.XData')) ];
            end
            scatterXYZ=scatterXYZ(2:end, :) ;  Indp3=Indp3(2:end)' ;
            
            EffectiveInput= true;
            %             OriP1= patchH{1,iBundle}.Vertices;
            OriP1=[nan,nan,nan] ;  Indp1 = 0;
            for bbi= 1 :length(iBundle)
                OriP1=[ OriP1 ; patchH{1,iBundle(bbi)}.Vertices ];
                Indp1=[ Indp1 ,repelem(bbi , size( patchH{1,iBundle(bbi)}.Vertices ,1)) ];
            end
            OriP1=OriP1(2:end, :) ;  Indp1=Indp1(2:end)' ;
            
            
            Gcenter= mean(reshape(PlotXYZ, size(PlotXYZ,1)*2,3)) ;
            if checkH_Rcenter.Value == 1 && ~isempty(obj.Rcenter_scatterH)
                Gcenter = [obj.Rcenter_scatterH.XData,obj.Rcenter_scatterH.YData,obj.Rcenter_scatterH.ZData ] ;
            end
            %--------
            %             return
            if checkH.Value==0   %translation
                LinearIncrement= str2double(editH.String)       ;
                if isempty(LinearIncrement)
                    LinearIncrement=1;
                end
                if  checkH2_GL.Value ==0
                    switch evn.Character
                        case 'a'
                            Vect = [-LinearIncrement 0 0]  ;
                        case 'q'
                            Vect = [LinearIncrement 0 0]  ;
                        case 's'
                            Vect = [0 -LinearIncrement 0]  ;
                        case 'w'
                            Vect =  [0 LinearIncrement 0]  ;
                        case 'd'
                            Vect = [0 0 -LinearIncrement]  ;
                        case 'e'
                            Vect = [0 0 LinearIncrement] ;
                        otherwise
                            EffectiveInput=false ;         return;
                    end
                    Vect=Vect' ;
                else  % if checkH2_GL.Value ==0  % local coordinate
                    switch evn.Character
                        case 'a'
                            Vect =  -obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 1)  ;
                        case 'q'
                            Vect =  obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 1)  ;
                        case 's'
                            Vect =  -obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 2)  ;
                        case 'w'
                            Vect =  obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 2)  ;
                        case 'd'
                            Vect =  -obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 3)  ;
                        case 'e'
                            Vect =  obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 3)  ;
                        otherwise
                            EffectiveInput=false ;         return;
                    end
                    Vect=Vect*LinearIncrement ;
                    
                    
                end  % G/L
                OriP1= OriP1 +  ones(size(OriP1,1),1)*Vect'  ;
                PlotXYZV2=PlotXYZV2 + ones(size(PlotXYZV2,1),1)*Vect' ;
                scatterXYZ=scatterXYZ+ ones(size(scatterXYZ,1),1)*Vect' ;
            else   %rotation
                RotaionIncrement= str2double(editH.String);
                RotaionIncrement=RotaionIncrement*pi/180;
                if isempty(RotaionIncrement)
                    RotaionIncrement=0.05;
                end
                
                theta=RotaionIncrement;
                switch evn.Character
                    case 'a'
                        if  checkH2_GL.Value ==0
                            RMat=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
                        else
                            LocalAxis= obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 1) ;
                            RMat = RotationAxisToRotaionMatrix( LocalAxis,theta*180/pi )     ;
                        end
                        OriP1=OriP1 - ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 - ones(size(PlotXYZV2,1),1)*Gcenter ;
                        OriP1=OriP1* RMat ;                            PlotXYZV2=PlotXYZV2*RMat;
                        OriP1=OriP1 + ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 + ones(size(PlotXYZV2,1),1)*Gcenter ;
                    case 'q'
                        if  checkH2_GL.Value ==0
                            RMat=[1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
                        else
                            LocalAxis= -obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 1) ;
                            RMat = RotationAxisToRotaionMatrix( LocalAxis,theta*180/pi )     ;
                        end
                        OriP1=OriP1 - ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 - ones(size(PlotXYZV2,1),1)*Gcenter ;
                        OriP1=OriP1* RMat ;                            PlotXYZV2=PlotXYZV2*RMat;
                        OriP1=OriP1 + ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 + ones(size(PlotXYZV2,1),1)*Gcenter ;
                    case 's'
                        if  checkH2_GL.Value ==0
                            RMat=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
                        else
                            LocalAxis= obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 2) ;
                            RMat = RotationAxisToRotaionMatrix( LocalAxis,theta*180/pi )     ;
                        end
                        
                        OriP1=OriP1 - ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 - ones(size(PlotXYZV2,1),1)*Gcenter ;
                        OriP1=OriP1* RMat ;                            PlotXYZV2=PlotXYZV2*RMat;
                        OriP1=OriP1 + ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 + ones(size(PlotXYZV2,1),1)*Gcenter ;
                    case 'w'
                        if  checkH2_GL.Value ==0
                            RMat=[cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
                        else
                            LocalAxis= -obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 2) ;
                            RMat = RotationAxisToRotaionMatrix( LocalAxis,theta*180/pi )     ;
                        end
                        OriP1=OriP1 - ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 - ones(size(PlotXYZV2,1),1)*Gcenter ;
                        OriP1=OriP1* RMat ;                            PlotXYZV2=PlotXYZV2*RMat;
                        OriP1=OriP1 + ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 + ones(size(PlotXYZV2,1),1)*Gcenter ;
                    case 'd'
                        if checkH2_GL.Value ==0
                            RMat=[cos(theta) -sin(theta) 0; sin(theta)  cos(theta)  0; 0 0 1];
                        else
                            LocalAxis= obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 3) ;
                            RMat = RotationAxisToRotaionMatrix( LocalAxis,theta*180/pi )     ;
                        end
                        
                        OriP1=OriP1 - ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 - ones(size(PlotXYZV2,1),1)*Gcenter ;
                        OriP1=OriP1* RMat ;                            PlotXYZV2=PlotXYZV2*RMat;
                        OriP1=OriP1 + ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 + ones(size(PlotXYZV2,1),1)*Gcenter ;
                    case 'e'
                        if checkH2_GL.Value ==0
                            RMat=[cos(theta) sin(theta) 0; -sin(theta)  cos(theta)  0; 0 0 1];
                        else
                            LocalAxis= -obj.containBundle{iBundle(1)}.TransformMatrix2(1:3, 3) ;
                            RMat = RotationAxisToRotaionMatrix( LocalAxis,theta*180/pi )     ;
                        end
                        OriP1=OriP1 - ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 - ones(size(PlotXYZV2,1),1)*Gcenter ;
                        OriP1=OriP1* RMat ;                            PlotXYZV2=PlotXYZV2*RMat;
                        OriP1=OriP1 + ones(size(OriP1,1),1)*Gcenter ;  PlotXYZV2=PlotXYZV2 + ones(size(PlotXYZV2,1),1)*Gcenter ;
                        %                    case 'r'
                        %                         rotate3d(src)
                    otherwise
                        EffectiveInput=false ;   return;
                end
                
                if checkH_Rcenter.Value == 1 && ~isempty(obj.Rcenter_scatterH)
                   Gcenter = [obj.Rcenter_scatterH.XData,obj.Rcenter_scatterH.YData,obj.Rcenter_scatterH.ZData ] ;
                end
                
                scatterXYZ=scatterXYZ- ones(size(scatterXYZ,1),1)*Gcenter ;
                scatterXYZ=scatterXYZ*RMat;
                scatterXYZ=scatterXYZ+ ones(size(scatterXYZ,1),1)*Gcenter ;
            end  % end of trans/rotate
            
            %%%% dispatch to graphical objects  Indp1, Indp2,Indp3
            for bbi= 1 :length(iBundle)
                patchH{1,iBundle(bbi)}.Vertices=OriP1(Indp1==bbi , :);
                
                nCol=size(PlotXYZV2,1)/2;
                
                PlotXYZV3= [PlotXYZV2(1:nCol,1), PlotXYZV2(nCol+1:end,1) , PlotXYZV2(1:nCol,2), PlotXYZV2(nCol+1:end,2),PlotXYZV2(1:nCol,3), PlotXYZV2(nCol+1:end,3)];
                PlotXYZV3TT= PlotXYZV3(Indp2==bbi ,: ) ;
                for spread=1: length( patchH{2,iBundle(bbi)})
                    [bbi,spread];
                    patchH{2,iBundle(bbi)}(spread).XData =   PlotXYZV3TT(spread,1:2);
                    patchH{2,iBundle(bbi)}(spread).YData =   PlotXYZV3TT(spread,3:4);
                    patchH{2,iBundle(bbi)}(spread).ZData =   PlotXYZV3TT(spread,5:6);
                end
                
                patchH{3,iBundle(bbi)}.XData=scatterXYZ(Indp3==bbi,1)';
                patchH{3,iBundle(bbi)}.YData=scatterXYZ(Indp3==bbi,2)';
                patchH{3,iBundle(bbi)}.ZData=scatterXYZ(Indp3==bbi,3)';
            end
            %----------------------
            %             return
            
            %Calculate Overall TF to V0, Update to table
            FinalCylXYZ=cell(size(V0CylinderXYZ));
            for gatherc=1:length( FinalCylXYZ)
                MM=zeros(length(patchH{2,gatherc}),6);
                for  nLine=1:length(patchH{2,gatherc})
                    MM(nLine,:)= [patchH{2,gatherc}(nLine).XData , patchH{2,gatherc}(nLine).YData,patchH{2,gatherc}(nLine).ZData];
                end
                FinalCylXYZ{gatherc}=MM;
            end
            bundlej =iBundle;
            %                 OP6=V0CylinderXYZ{bundlej};   %alway compare to the very beginning
            %                 OP3=[[OP6(:,1);OP6(:,2)], [OP6(:,3);OP6(:,4)], [OP6(:,5);OP6(:,6)]] ;
            
            for Bi= 1 :length(bundlej)
                OP3=PlotXYZV2Ini{Bi}  ;
                
                NewP6=FinalCylXYZ{bundlej(Bi)};
                NewP3=[[NewP6(:,1);NewP6(:,2)], [NewP6(:,3);NewP6(:,4)], [NewP6(:,5);NewP6(:,6)]] ;
                [regParams,~,~]=absor(OP3',NewP3');
                %
                tH.Data=regParams.M;
                [regParams.R, regParams.t;0,0,0,1];
%                QQ= obj.containBundle{iBundle(Bi)}.TransformMatrix2 
                obj.containBundle{iBundle(Bi)}.TransformMatrix2= [regParams.R, regParams.t;0,0,0,1]*obj.containBundle{iBundle(Bi)}.TransformMatrix2 ;
%                     QQ2=   obj.containBundle{iBundle(Bi)}.TransformMatrix2
                XYZALL= [  patchH{3,iBundle(Bi)}.XData; patchH{3,iBundle(Bi)}.YData; patchH{3,iBundle(Bi)}.ZData];
                xyz=mean(XYZALL,2);
                textBundle{iBundle(Bi)}.Position=xyz';
            end
            
            
            %             EffectiveInput
            if EffectiveInput==1 && ~strcmp(evn.Key,'?')
                if checkH.Value==0
                    ThisCommand =AssemblyCommand('XYZ_input', checkH.Value, checkH2_GL.Value, LinearIncrement, evn.Character ,bundlej)  ;
                else
                    ThisCommand =AssemblyCommand('XYZ_input', checkH.Value, checkH2_GL.Value, RotaionIncrement/pi*180, evn.Character,bundlej )   ;
                end
                obj.TransformHistroy{end+1}= ThisCommand  ;
            end
            
            
        end
        %
        
        
        function  plotXYZ2=plotScafR_BaseByBase(obj,ScafForScatter)
            
            %             AllGXYZ=cell(1,length(obj.containBundle));
            %             for k=1:length(obj.containBundle)
            %             AllGXYZ{k}=obj.containBundle{k}.CylinderXYZGlobal ;
            %             end
            %
            %             plotXYZ=zeros(size(ScafForScatter,1) ,3);  % C5
            %             for k=1:size(SacfR,1)
            %                [~, ind] = ismember(
            %
            %                bundle=SacfR(k,1);  Cyl=SacfR(k,2);
            %                alpha=SacfR(k,3)- obj.containBundle{bundle}.Zbase1(Cyl);
            %                beta=obj.containBundle{bundle}.Zbase2(Cyl)-SacfR(k,3);
            %
            %
            %                P= AllGXYZ{bundle}(Cyl,1:3);
            %                Q=AllGXYZ{bundle}(Cyl,4:6);
            %                XYZ=(beta*P + alpha*Q )/(alpha+beta);
            %                plotXYZ(k,:)=XYZ;
            %             end
            %
            %             tic
            %             tic
            %             plotXYZ=zeros(size(ScafForScatter,1) ,3);  % C5   % old way
            %             query base by base , not efficient 08/26/2019
            %             for k=1:size(ScafForScatter,1)
            %                 ColRow =  obj.RelateTable(obj.RelateTable(:,5)==ScafForScatter(k,1),6:7);  ColRow=ColRow(1,:) ;
            %                 indB =find(obj.RelateTable(:,5)==ScafForScatter(k,1)) ; indB=indB(1) ;
            %                 bundle= obj.RelateTable(indB,1) ;
            %                 Bundle= obj.containBundle{bundle };
            %                 if mod(obj.RelateTable(indB,4) ,2) ==0
            %                     SimilarDirCylders = and( mod(obj.RelateTable(:,4) ,2)==0 , obj.RelateTable(:,1)==bundle) ;
            %                 else
            %                     SimilarDirCylders = and( mod(obj.RelateTable(:,4) ,2)==1 , obj.RelateTable(:,1)==bundle) ;
            %                 end
            %                 fSimilarDirCylders = find(SimilarDirCylders) ; fSimilarDirCylders=fSimilarDirCylders(1) ;
            %                 ThefSimilarDirCylder= obj.RelateTable(fSimilarDirCylders,2) ;
            %
            %                 RefCyl = and(obj.RelateTable(:,1)==bundle, obj.RelateTable(:,2)~=-1  ) ;
            %                 fRefCyl= find(RefCyl) ;  fRefCyl=fRefCyl(1) ;
            %                 BasesArr=ScafForScatter(k,2) ;
            %
            %                 InplaneXY=Bundle.findExtraCylInplanePosition(  obj.RelateTable, ColRow,fRefCyl)   ;%
            %                 Global_XYZ=Bundle.HelixRegardlessCylinder(0,1,InplaneXY,BasesArr,ThefSimilarDirCylder) ;
            %                 plotXYZ(k,:) = Global_XYZ(1:3) ;
            %
            %             end
            %             toc
            %
            
            %-------improve efficiency with vectorize functions, 08/26/2019
            BCB =ConverC5Routing_BCB(obj,ScafForScatter) ;
            [~,b] =   ismember( ScafForScatter(:,1) , obj.RelateTable(:,5) ) ;
            ColRows  =  obj.RelateTable(b,6:7) ;   % vectorized
            plotXYZ2=zeros(size(ScafForScatter,1) ,3);  % C5
            for Bi = 1 : max(BCB(:,1))                 % query bundle by bundle, instead of base by base
%                 Bi
                BaseIndInThisBundle = BCB(:,1)==Bi ;  ffBaseIndInThisBundle = find(BaseIndInThisBundle);
                bundle= Bi ;
                Bundle= obj.containBundle{bundle };
                RefCyl = and(obj.RelateTable(:,1)==bundle, obj.RelateTable(:,2)~=-1  ) ;
                fRefCyl= find(RefCyl) ;  fRefCyl=fRefCyl(1) ;
                
                InplaneXY=Bundle.findExtraCylInplanePosition(  obj.RelateTable, ColRows(BaseIndInThisBundle,: ),fRefCyl) ;
                
                [~,bCR] =   ismember(  ColRows(BaseIndInThisBundle,: ) , obj.RelateTable(:,6:7) ,'rows' ) ;
                EvenOrOldCylinder  =  obj.RelateTable(bCR,4) ;   % vectorized
                EvenOrOldCylinder= mod(EvenOrOldCylinder,2)==0 ;
                
                ThefSimilarDirCylder_A= unique( obj.RelateTable(bCR(EvenOrOldCylinder) ,2)')  ; % query Agroup and BGroup
                %
                %                 if isempty(ThefSimilarDirCylder_A)
                %                 sdsf=3
                %                 end
                if ~isempty(ThefSimilarDirCylder_A)
                    ThefSimilarDirCylder_A=ThefSimilarDirCylder_A(1) ;
                else
                    ThefSimilarDirCylder_A=Bundle.AGroup(1) ;
                end
                ThefSimilarDirCylder_B= unique( obj.RelateTable(bCR(~EvenOrOldCylinder) ,2)') ;
                if ~isempty(ThefSimilarDirCylder_B)
                    ThefSimilarDirCylder_B=ThefSimilarDirCylder_B(1) ;
                else
                    ThefSimilarDirCylder_B=Bundle.BGroup(1) ;
                end
                BCB_inthisBundle = BCB(BaseIndInThisBundle,3) ;
                
                Global_XYZ_A=Bundle.HelixRegardlessCylinder(0,1,InplaneXY(EvenOrOldCylinder,:), BCB_inthisBundle(EvenOrOldCylinder), ThefSimilarDirCylder_A) ;
                Global_XYZ_B=Bundle.HelixRegardlessCylinder(0,1,InplaneXY(~EvenOrOldCylinder,:),  BCB_inthisBundle(~EvenOrOldCylinder), ThefSimilarDirCylder_B) ;
                
                plotXYZ2(ffBaseIndInThisBundle(EvenOrOldCylinder) ,:) = Global_XYZ_A(:, 1:3) ;
                plotXYZ2(ffBaseIndInThisBundle(~EvenOrOldCylinder) ,:) = Global_XYZ_B(:, 1:3) ;
                
                %                 figure; scatter3(Global_XYZ_A(:,1),Global_XYZ_A(:,2),Global_XYZ_A(:,3)); hold on ; scatter3(Global_XYZ_B(:,1),Global_XYZ_B(:,2),Global_XYZ_B(:,3));
            end
            %                         toc
            %                         dsfsf=3;
        end
        
        function transXYZ = Get_transScafXYZ(obj, SaveXYZ)
            % use for obtain scaffold routing in extrusion model. Use for
            % extrusion only. SaveXYZ: unlock codes in GUI_conti.. to get.
            
            BCB_scaf = obj.ConverC5Routing_BCB(obj.ScafAllBase{1}) ;
            transXYZ = zeros( size(BCB_scaf )) ; % Nx3
            for k =1: size(transXYZ ,1 )
                BCBi  = BCB_scaf(k,:) ;
                AllCyl_xyz =  SaveXYZ{ BCBi(2)} ;
                
                AllbaseOnCyl = BCB_scaf(BCB_scaf(:,2)==BCBi(2) ,: ) ;
                RefXYZ = sortrows(AllbaseOnCyl,[1 3]) ;
                [ ~ ,s]= ismember(BCBi , RefXYZ ,'rows' ) ;
                s=s/size(RefXYZ,1) ;
                RefX=linspace(0 , 1 , size(AllCyl_xyz ,1)) ;
                NewXYZ = interp1(RefX  ,AllCyl_xyz, s) ;
                transXYZ(k,:) = NewXYZ ;
            end
            
            for k =1:3                
                dXYZ = diff(transXYZ) ;
                ds = sqrt(dXYZ(:,1).^2 + dXYZ(:,2).^2 +dXYZ(:,3).^2 );
                Inds = find(ds > 1) ; RemoveInd =[Inds; Inds+1] ;
                transXYZ(RemoveInd ,:  ) =[];
            end                        
        end
        
        function transXYZ_cell = Get_transStapXYZ(obj, SaveXYZ ,h_staple)
            BCB_stap = cell(length(obj.StapAllBase) ,1) ;
            for k =1:length(obj.StapAllBase)
            BCB_stap{k}=obj.ConverC5Routing_BCB(obj.StapAllBase{k}) ;            
            end
            
            StackAll = cell2mat( BCB_stap) ;
            
            transXYZ_cell=cell(length(obj.StapAllBase) ,1) ;
          
            for k =1:length(obj.StapAllBase)
                
                BCB_oneStrand  = BCB_stap{k} ;
                transXYZ = zeros( size(BCB_oneStrand )) ; % Nx3
                for j = 1: size(BCB_oneStrand ,1 )
                BCBi  = BCB_oneStrand(j,:) ;
                AllCyl_xyz =  SaveXYZ{ BCBi(2)} ;
            
                
                AllbaseOnCyl = StackAll(StackAll(:,2)==BCBi(2) ,: ) ;
                RefXYZ = sortrows(AllbaseOnCyl,[1 3]) ;
                [ ~ ,s]= ismember(BCBi , RefXYZ ,'rows' ) ;
                s=s/size(RefXYZ,1) ;
                RefX=linspace(0 , 1 , size(AllCyl_xyz ,1)) ;
                NewXYZ = interp1(RefX  ,AllCyl_xyz, s) ;
                transXYZ(j,:) = NewXYZ ;            
               
                end
                transXYZ_cell{k}=transXYZ ;
                %---------
                for dsk =1:1
                    dXYZ = diff(transXYZ) ;
                    ds = sqrt(dXYZ(:,1).^2 + dXYZ(:,2).^2 +dXYZ(:,3).^2 );
                    Inds = find(ds > 1) ; RemoveInd =[Inds; Inds+1] ;
                    transXYZ(RemoveInd ,:  ) =[];
                end
                
                transXYZ(1,:)=[];
                transXYZ(end,:)=[];
                %
                
                
                h_staple(end-k+1).XData = transXYZ(:,1)' ;
                h_staple(end-k+1).YData = transXYZ(:,2)' ;
                h_staple(end-k+1).ZData = transXYZ(:,3)' ;
                

            end
            
        end
        
        
        
        function BCB =ConverC5Routing_BCB(obj,C5Routing)
            
            BCB=zeros( size(C5Routing,1) ,3) ;
            BCB(:, end) =     C5Routing(:, end)  ;
            [~,b] =   ismember( C5Routing(:,1) , obj.RelateTable(:,5) ) ;
            BC =   obj.RelateTable(b,1:2) ;
            BCB(:, 1:2) = BC ;
        end
        
        function surfH=plotScafR_cylindermodelAllBase(obj,whichScafSource, scafInd )
            if nargin==1
                whichScafSource=1;
            end
            AllGXYZ=cell(1,length(obj.containBundle));
            for k=1:length(obj.containBundle)
                AllGXYZ{k}=obj.containBundle{k}.CylinderXYZGlobal ;
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
            if whichScafSource==1
                SacfR=obj.ScafRouting{scafInd} ;
            elseif  whichScafSource==2
                SacfR=obj.Scaf_fromCadDOM ;
            elseif whichScafSource==3
                SacfR=obj.Scaf_fromJSON ;
            end
            SacfR = interpolateBase_ThreeCol( SacfR ) ;   % BCB notation
            
            
            [ScafForScatter, ia]= setdiff(interpolateBase( obj.scafC5{scafInd} ) , obj.skipBase,'rows' ,'stable') ;  %C5 notation
            %  consider skip in SQ
            SacfR=SacfR(ia,:) ;
            %             SacfR=ScafForScatter ;
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
            
            % SSplotXYZ=size(plotXYZ) ;
            % plot3(plotXYZ(:,1), plotXYZ(:,2), plotXYZ(:,3) )
            x = plotXYZ(:,1)';
            y = plotXYZ(:,2)';
            z = plotXYZ(:,3)';
            %             col = (1:size(plotXYZ,1))*1000;  % This is the color
            col = ((1:size(plotXYZ,1))-1)/( size(plotXYZ,1)-1);  % This is the color
            surfH=surface([x;x],[y;y],[z;z],[col;col], 'facecol','no', 'edgecol','interp', 'linew',2);
            
            %     DSh=plot3(XYZdata(:,1),-XYZdata(:,2),XYZdata(:,3),'Linewidth',2);   %draw route
            scatter3(plotXYZ(1,1),plotXYZ(1,2),plotXYZ(1,3),'s','filled','SizeData',100);   %mark head
            scatter3(plotXYZ(end,1),plotXYZ(end,2),plotXYZ(end,3),'d','filled','SizeData',100);  %mark tail
            axis equal;grid on;  xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
            surfH.UserData.SacfR=SacfR ;
        end % end of plotScafR_cylindermodelAllBase
        
        function surfH=plotScafR_cylindermodel(obj,whichScaf)
            if nargin==1
                whichScaf=1;
            end
            AllGXYZ=cell(1,length(obj.containBundle));
            for k=1:length(obj.containBundle)
                AllGXYZ{k}=obj.containBundle{k}.CylinderXYZGlobal ;
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
            if whichScaf==1
                SacfR=obj.ScafRouting ;
            elseif  whichScaf==2
                SacfR=obj.Scaf_fromCadDOM ;
            elseif whichScaf==3
                SacfR=obj.Scaf_fromJSON ;
            elseif whichScaf==4
                SacfR=obj.ScafRouting ;
            end
            %             SacfR=ScafForScatter ;
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
            
            % SSplotXYZ=size(plotXYZ) ;
            % plot3(plotXYZ(:,1), plotXYZ(:,2), plotXYZ(:,3) )
            x = plotXYZ(:,1)';
            y = plotXYZ(:,2)';
            z = plotXYZ(:,3)';
            col = (1:size(plotXYZ,1))*1000;  % This is the color
            surfH=surface([x;x],[y;y],[z;z],[col;col], 'facecol','no', 'edgecol','interp', 'linew',2);
            
            %     DSh=plot3(XYZdata(:,1),-XYZdata(:,2),XYZdata(:,3),'Linewidth',2);   %draw route
            scatter3(plotXYZ(1,1),plotXYZ(1,2),plotXYZ(1,3),'s','filled','SizeData',100);   %mark head
            scatter3(plotXYZ(end,1),plotXYZ(end,2),plotXYZ(end,3),'d','filled','SizeData',100);  %mark tail
            axis equal;grid on;  xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
            surfH.UserData.SacfR=SacfR ;
        end % end of plotScafR_cylindermodel
        
        function surfH=plotScafR_cylindermodelMulti(obj,whichScaf, GradOrIsoColor,GivenScafR)
            if nargin==1
                whichScaf=1;
                GradOrIsoColor='Gradient' ;
            end
            AllGXYZ=cell(1,length(obj.containBundle));
            for k=1:length(obj.containBundle)
                AllGXYZ{k}=obj.containBundle{k}.CylinderXYZGlobal ;
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
            if whichScaf==1
                SacfC=obj.ScafRouting ;
            elseif  whichScaf==2
                SacfC=obj.Scaf_fromCadDOM ;
            elseif whichScaf==3
                SacfC=obj.Scaf_fromJSON ;
            elseif whichScaf==4
                SacfC=obj.ScafRouting ;
            elseif whichScaf==5 
                SacfC= GivenScafR ;
            end
            %             SacfR=ScafForScatter ;
            CC= get(gca,'defaultAxesColorOrder') ; QCC=CC(4,:) ; 
            %             TempCC = [ 1,0,0  ; 0 ,1,0 ; 0,0 1 ; 1, 0.5,0 ;  0,1,1; 1,0 ,1];
            CC(2,:)=[ 0.7500  0  0.7500] ;CC(3,:)=[]; CC(3,:)=[];  
            CC(1,:)= [0.3010    0.7450    0.9330];         CC(3,:)=[0    0.5000         0];CC;
            CC=[CC;QCC];CC(4,:)=[];
            set(gca,'ColorOrder',CC) ;
            %             set(0,'defaultAxesColorOrder',TempCC*0.95) ;
            SaveXYZ= zeros( length(SacfC) ,6 ) ;
            cla ;hold on;
            surfH= cell(length(SacfC),1) ;
            for kc= 1 :length(SacfC)
                SacfR= SacfC{kc} ;
                plotXYZ=zeros(size(SacfR));
                for k=1:size(SacfR,1)
                    bundle=SacfR(k,1);  Cyl=SacfR(k,2);
                    alpha=SacfR(k,3)- obj.containBundle{bundle}.Zbase1(Cyl);
                    
                    beta=obj.containBundle{bundle}.Zbase2(Cyl)-SacfR(k,3);
                    P= AllGXYZ{bundle}(Cyl,1:3);
                    Q= AllGXYZ{bundle}(Cyl,4:6);
                    XYZ=(beta*P + alpha*Q )/(alpha+beta);
                    plotXYZ(k,:)=XYZ;
                end
                
                % SSplotXYZ=size(plotXYZ) ;
                % plot3(plotXYZ(:,1), plotXYZ(:,2), plotXYZ(:,3) )
                x = plotXYZ(:,1)';
                y = plotXYZ(:,2)';
                z = plotXYZ(:,3)';
                col = ((1:size(plotXYZ,1))-1)/( size(plotXYZ,1)-1);  % This is the color
                %              col = (1:size(plotXYZ,1))*1000;  % This is the color
                %             ssdfsf=3 ;
                % 'Gradient','IsoColor'
                switch GradOrIsoColor
                    case 'Gradient'
                        surfH{kc}=surface([x;x],[y;y],[z;z],[col;col], 'facecol','no', 'edgecol','interp', 'linew',2);
                    case 'IsoColor'
                        surfH{kc}=plot3(x,y,z ,'LineWidth',2);
                        %                         surfH.Color = [surfH.Color ,0.7];
                end
                %                 ax = gca;
                % ax.ColorOrderIndex = 2;
                SaveXYZ(kc,:) = [plotXYZ(1,:),plotXYZ(end,:)] ;
                %     DSh=plot3(XYZdata(:,1),-XYZdata(:,2),XYZdata(:,3),'Linewidth',2);   %draw route
                %                 scatter3(plotXYZ(1,1),plotXYZ(1,2),plotXYZ(1,3),'sr','filled','SizeData',100);   %mark head
                %                 scatter3(plotXYZ(end,1),plotXYZ(end,2),plotXYZ(end,3),'db','filled','SizeData',100);  %mark tail
                
                surfH{kc}.UserData.SacfR=SacfR ;
            end
            scatter3(SaveXYZ(:,1),SaveXYZ(:,2),SaveXYZ(:,3),'sr','filled','SizeData',100);   %mark head
            scatter3(SaveXYZ(:,4),SaveXYZ(:,5),SaveXYZ(:,6),'db','filled','SizeData',100);  %mark tail
            
            
            %             set(0, 'DefaultAxesColorOrder', 'factory');
            axis equal;grid on;  xlabel('X-axis') ;ylabel('Y-axis');zlabel('Z-axis') ;
            title( sprintf(' # of scaffold = %i ',kc ))
%             surfH.UserData.SacfR=SacfR ;
        end % end of plotScafR_cylindermodelMulti
        
        function pH=plotScaf3Dposition(obj,Xovers)
            figure(3252);clf ; hold on ;
            obj.plotScafR_cylindermodelMulti(2, 'Gradient') ;
            AllGXYZ=cell(1,length(obj.containBundle));
            for k=1:length(obj.containBundle)
                AllGXYZ{k}=obj.containBundle{k}.CylinderXYZGlobal ;
            end
            
            
            pH=cell(size(Xovers,1) ,1) ;
            for kc= 1 :size(Xovers,1)
                SacfR= Xovers(kc,:) ; SacfR=reshape(SacfR',3,4)' ;
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
                
                plotXYZ=[plotXYZ(1:2,:) ; [nan nan nan] ;plotXYZ(3:4,:) ];
                % SSplotXYZ=size(plotXYZ) ;
                % plot3(plotXYZ(:,1), plotXYZ(:,2), plotXYZ(:,3) )
                x = plotXYZ(:,1)';
                y = plotXYZ(:,2)';
                z = plotXYZ(:,3)';
                pH{kc}=plot3(x,y,z, 'k','LineWidth',1.5);
            end
            
            
            
        end % end of plotScaf3Dposition
        
        
        
        
        function obj=ApartStaplesSimple(obj, varargin)
            %                         skipBase= GetHyperB.skipBase ;
            %             if strcmp(obj.StapOption.halfXover,'yes')  && strcmp(obj.StapOption.type,'straightcut')
            DoHalf=1;
            if nargin==1
                StapleCell= obj.StapList2 ;
            elseif nargin==2  % debug use
                StapleCell=    varargin{1} ;
            elseif nargin==3  % Don't do halxover after first cutting
                StapleCell= obj.StapList2 ;
                DoHalf=0;
            end
            if strcmp(obj.StapOption.halfXover,'yes') && DoHalf==1
                scriptcase=1;%----------------halfCrossover
                AddHalfXoverV3LargeCrossSec;   %input: StapleCell
            end
            obj.StapList2 =  StapleCell;
            %             end
            
            
            StapleCell;
            NewStapleCell =cell(1000,1) ; iC=1;
%             Bound=[60,90];  mB =75;
            
            for k=1: length(StapleCell)
                Bound=[35,60];  mB =50;
%                  Bound=[60,90];  mB =75;
                BaseRout = setdiff( interpolateBase( StapleCell{k} ) , obj.skipBase ,'rows','stable') ;
                %                 [WW,AlwaysNotCut]=intersect( BaseRout, StapleCell{k},'rows' ) ;
                [WW2,AlwaysNotCut]=intersect( BaseRout, [StapleCell{k}; obj.scafC5All ],'rows' ) ;
                
                AlwaysNotCut= repmat(AlwaysNotCut, 1 , 5) + repmat(-2:1:2 ,   size(AlwaysNotCut,1) ,1) ;  % nicks and thier neighbors
                %-------- +-2nt from xover is the limit, since
                % 5(+-2)+2(xover itself) =7 is HC domain. Increasing the
                % tolerance result in no position to cut for inner HC
                % cylinder.
                AlwaysNotCut=AlwaysNotCut(:); AlwaysNotCut(AlwaysNotCut<=0)=[] ;AlwaysNotCut(AlwaysNotCut>size(BaseRout,1))=[];
                OkeyToCutInd =setdiff(1:size(BaseRout,1) , AlwaysNotCut) ;
                nw=1 ;
                
                %------debug use. Monitor where can be cut.
                %                 NoCut =zeros(size(CutPlace)) ;
                %                 NoCut(AlwaysNotCut)=-1 ;
                %                 [BaseRout , CutPlace, NoCut]
                %------
                
                if size(BaseRout,1)<= Bound(2)
                    BaseArray= zeros(size(BaseRout,1) ,1 ) ;
                else
                    while 1
                        CutPlace = zeros(size(BaseRout,1) ,1 ) ;
                        N_CutPoint =  round(size(BaseRout,1)/mB ) + randi(2 +round(size(BaseRout,1)/50)  ) -2  ;   % flexible +-1
                        ToCut = randi(numel(OkeyToCutInd),N_CutPoint ,1 );
                        CutPlace(OkeyToCutInd(ToCut)  ) = 1 ; CutPlace(1) =0 ;
                        CutPlace(AlwaysNotCut)=0;
                        for adjustI=1:3
                            BaseArray = cumsum(CutPlace) ;
                            [a,~]=hist(BaseArray,unique(BaseArray))  ;
                            TooShort = find(a < 20 ) ; TooShort=setdiff(TooShort,1);
                            Inds =find(CutPlace)   ;
                            CutPlace( Inds(TooShort-1)) = 0 ;
                        end
                        
                        if  and(sum( and( a>=Bound(1) ,a<=Bound(2)) ) ==length(a),~isempty(a))  ||nw>=2000*3
                            break
                        end
                        nw=nw+1 ;
                        if nw==1000*3
                            Bound=[25,60];   mB=35;
                            fprintf('Cutting loops up to half. Loose Criterion. \n')
                            %                             StapleCell{k}
                        end
                        if nw==1800*3
                            Bound=[20,150];   mB=45;
                            fprintf('Cutting loops close to max. Loose Criterion. \n')
                        end
                    end
                end
                
                for j= 0:max(BaseArray)
                    NewOneStap = BaseRout( BaseArray==j ,:) ;
                    CornerRep = ConvertCornerRep( NewOneStap ) ;
                    NewStapleCell{iC} =CornerRep ;iC=iC+1 ;
                end
            end
            NewStapleCell=NewStapleCell(1:iC-1) ;
            obj.StapList3=NewStapleCell;
            fprintf('Finished ApartStaplesSimple\n')
        end % end of ApartStaplesSimple
        
        
        function obj=ApartStaples(obj)
            
            if strcmp(obj.StapOption.halfXover,'yes')  && strcmp(obj.StapOption.type,'straightcut')
                scriptcase=1;%----------------halfCrossover
                StapleCell= obj.StapList2 ;
                AddHalfXoverV3LargeCrossSec;   %input: StapleCell
                obj.StapList2 =  StapleCell;
            end
            
            %         obj.StapList2{42}=obj.StapList2{42}(1:152,:);
            [ Res ] = showstapins( obj.StapList2 );
            Bound=[30,50];    DefLongSegThan=8;  %value, beyond which is considered as long edge
            meanbound=mean(Bound);
            %               DrawStapp3(obj.StapList2,obj.RelateVec,1)
            newcell=cell(300,1);  NewCell_count=1;
            for iStap=1:    size(Res.IniSegLength,1)
                %                 iStap
                if Res.InitialLength(iStap)< 60  || length( Res.IniSegLength{iStap})==2     %if too short, remain old staples
                    newcell(NewCell_count)=obj.StapList2(iStap);
                    NewCell_count=NewCell_count+1;
                    continue
                end
                if Res.InitialLength(iStap)>1000
                    loop2Reset=200;
                else
                    loop2Reset=50;
                end
                Ex=Res.IniSegLength{iStap} ;
                %                if   iStap==4
                %                    sdfsf=3
                %                end
                Ex2=Ex;
                DigOneForLong=Ex>DefLongSegThan+4;
                %------
                CountIntervals=0;trigger=0;
                for CI=1:length(DigOneForLong)  %Count Number of intervals
                    if DigOneForLong(CI)==0 && trigger==0
                        CountIntervals=CountIntervals+1;
                        trigger=1;
                    end
                    if DigOneForLong(CI)==1
                        trigger=0 ;%reset
                    end
                end
                %                 Qsum=cumsum(~DigOneForLong);
                
                %                 IntervalLengthAs7s=Qsum(Qsum(2:end)==Qsum(1:end-1)); TT=IntervalLengthAs7s;
                numberGroup=  bwlabel(~DigOneForLong)    ;
                IntervalLengthAs7s=zeros(1,max(numberGroup));
                for kG=1:max(numberGroup)
                    IntervalLengthAs7s(kG)=sum(numberGroup==kG);
                end
                Intervalstart= zeros(1,max(numberGroup));
                for kG=1:max(numberGroup)
                    AskGindex=find(numberGroup==kG);AskGindex=AskGindex(1);
                    Intervalstart(kG)=AskGindex;
                end
                NotDigOneForLong=~DigOneForLong;
                
                Interval.L=IntervalLengthAs7s;
                Interval.SP=Intervalstart;
                Interval.N=length(IntervalLengthAs7s);
                
                if  Interval.N==0   %|| iStap==2
                    sdfsf=324;
                    
                end
                
                LongEdges=  bwlabel(DigOneForLong)  ;
                LongGInfo=zeros(max(LongEdges),3);  %[ L (unit:nt),L (eges), SP(unit:Edge in Ex)]
                for k=1:size(LongGInfo,1)
                    LongGInfo(k,1)=  sum(Ex(find(LongEdges==k)));
                    LongGInfo(k,2)=sum(LongEdges==k);
                    QW1=find(LongEdges==k);
                    LongGInfo(k,3)=QW1(1);
                end
                LongInterval=find(LongGInfo(:,1)>meanbound); OLI=LongInterval;
                
                %                 if ~isempty(LongInterval)
                %                     sfg=3
                %                 end
                %
                nloop=0;
                Saved.Lengths=[];
                Saved.GlobalCut=[];
                Saved.Val=1000000;
                
                HeadTailisLong=~[NotDigOneForLong(1),NotDigOneForLong(end)];
                CompL=[0,0];
                if HeadTailisLong(1)==1;  CompL(1)=Ex(1);end
                if HeadTailisLong(2)==1;  CompL(2)=Ex(end);end
                
                HeadAndTail=[1,length(Interval.L)];  OHT=HeadAndTail;
                MiddleIntervals=[2:1:length(Interval.L)-1]; OMI=MiddleIntervals;
                SaveCase=ones(size(DigOneForLong));
                cc=1;
                
                Tesin=~DigOneForLong;
                TestLforCutAll=CalculateSegLengths(obj,Ex,Tesin);
                BestTargetVal=sum(TestLforCutAll>Bound(2));
                Ex;
                
                Intervalselect=cell(size(Interval.L));
                GlobalCut=-(int8(DigOneForLong));
                reset=1;
                saveplot=zeros(5000,4);
                
                %                 ddhandle=gobjects(100,3);
                %                 F(1500) = struct('cdata',[],'colormap',[]);
                %                 cutmoreLong=0;
                %                 pause
                if sum(Ex<=DefLongSegThan)==0  || Interval.N==0
                    HeadAndTail=[];
                    OHT=HeadAndTail;
                end
                %               HeadAndTail
                while 1
                    cc=cc+1;
                    %                     HeadAndTail
                    for SPCareofHead=1:length(HeadAndTail)   %special care of head and tai;
                        
                        %                         iStap
                        %                         HeadAndTail
                        %                         SPCareofHead
                        if  Interval.SP(HeadAndTail(SPCareofHead))==1 || -1+Interval.SP(HeadAndTail(SPCareofHead))+ Interval.L(HeadAndTail(SPCareofHead))==length(Ex)
                            NNa=Interval.L(HeadAndTail(SPCareofHead));
                            M = (dec2bin(0:(2^NNa)-1)=='1') ;
                            NumberofInsertBP=0: ceil(Interval.L(HeadAndTail(SPCareofHead))*8/meanbound);
                            RandPickCandidatesHT= M(ismember(sum(M,2),NumberofInsertBP),:);  %setup number of true( break points)
                            QWE= Interval.SP(HeadAndTail(SPCareofHead)):Interval.SP(HeadAndTail(SPCareofHead))+Interval.L(HeadAndTail(SPCareofHead))-1;
                            %                            if NNa*DefLongSegThan+CompL(SPCareofHead)+DefLongSegThan+Ex(QWE(end))*rand()<Bound(1)
                            %                            Pick=  RandPickCandidatesHT(1,:);
                            %                            else
                            Pick=  RandPickCandidatesHT(randi(size(RandPickCandidatesHT,1)),:);
                            %                            end
                            
                            %                            dfsf=234
                        else
                            %start with long edge or end with
                            NNa=Interval.L(HeadAndTail(SPCareofHead));
                            M = (dec2bin(0:(2^NNa)-1)=='1') ;
                            NumberofInsertBP=1: ceil(Interval.L(HeadAndTail(SPCareofHead))*8/meanbound);   % at least cut one
                            RandPickCandidatesHT= M(ismember(sum(M,2),NumberofInsertBP),:);  %setup number of true( break points)
                            QWE= Interval.SP(HeadAndTail(SPCareofHead)):Interval.SP(HeadAndTail(SPCareofHead))+Interval.L(HeadAndTail(SPCareofHead))-1;
                            %                            if NNa*DefLongSegThan+CompL(SPCareofHead)+DefLongSegThan+Ex(QWE(end))*rand()<Bound(1)
                            %                            Pick=  RandPickCandidatesHT(1,:);
                            %                            else
                            Pick=  RandPickCandidatesHT(randi(size(RandPickCandidatesHT,1)),:);
                            %                            end
                        end
                        
                        GlobalCut(QWE)=Pick;
                        Intervalselect{HeadAndTail(SPCareofHead)}=Pick;
                    end
                    
                    for N_int2=1:length(MiddleIntervals)   %for middle interval
                        N_int=MiddleIntervals(N_int2);
                        NN=Interval.L(N_int);         %number of digits in this interval
                        seeds=NN-1:-1:5;
                        seedinDec=2.^seeds+1; seedinDec=[seedinDec 1];
                        for k=1:length(seedinDec)
                            QQseed=seedinDec(k);
                            while  QQseed*2<(2^NN)
                                seedinDec=[seedinDec QQseed*2];
                                QQseed=QQseed*2 ;
                            end
                        end
                        RandPickCandidates=dec2bin(seedinDec,NN)=='1';
                        Pick2=  RandPickCandidates(randi(size(RandPickCandidates,1)),:);
                        GlobalCut( Interval.SP(N_int):Interval.SP(N_int)+Interval.L(N_int)-1)=Pick2;
                        Intervalselect{N_int}=Pick2;
                    end
                    
                    %                     LongGInfo;  %[ L (unit:nt),L (deges), SP(unit:Edge in Ex)]
                    %                     LongEdges;
                    for LongIi=1:length(LongInterval)
                        NN= LongGInfo(LongInterval(LongIi),2);   %number of digits
                        ran20=20*(0.5-rand);
                        %                       NumberofCut= floor((LongGInfo(LongInterval(LongIi),1)+DefLongSegThan)/(meanbound-cutmoreLong));
                        NumberofCut= floor((LongGInfo(LongInterval(LongIi),1)+DefLongSegThan) /(meanbound-ran20));
                        %                          NumberofCut=0: ceil(Interval.L(HeadAndTail(SPCareofHead))*8/meanbound);
                        %                       NN
                        M = (dec2bin(0:(2^NN)-1)=='1') ;
                        QuF=sum(M,2)<=NumberofCut;
                        RandPickCandidatesL=M(QuF,:) ;
                        %                       if  isempty(RandPickCandidatesL)
                        %                           sdfsf=1
                        %                       end
                        %                        size(RandPickCandidatesL,1);
                        PickL=RandPickCandidatesL(randi(size(RandPickCandidatesL,1)) ,:);
                        SP=LongGInfo(LongInterval(LongIi),3);
                        PutInPosition=SP:SP+LongGInfo(LongInterval(LongIi),2)-1;
                        GlobalCut(PutInPosition)=PickL;
                    end
                    if sum(GlobalCut==1)==0
                        nloop= nloop+1   ;
                        if  nloop<500
                            continue;
                        else
                            sdfsf= 2111;
                        end
                    else
                        Lengths=CalculateSegLengths(obj,Ex,GlobalCut);   %CALCULATE LENGTH
                        BadLengh=or(Lengths>Bound(2),Lengths<Bound(1))';
                        
                        SaveCase=union(SaveCase,GlobalCut,'rows');   %
                        OutOfBound=Lengths>Bound(2); twoshort=Lengths<Bound(1);
                        
                        LengHasTrouble=or(OutOfBound,twoshort);
                        dLofTooShort=Lengths(twoshort)-Bound(1);
                        %                     TheLessTheBetter=sum(OutOfBound)+0.2*sum(twoshort) ;
                        %                     TheLessTheBetter=10*sum(OutOfBound)+sum(abs(dLofTooShort))/100 ;
                        %                     TheLessTheBetter= sum(abs(dLofTooShort))/100 ;
                        
                        a1=1 ; a2=4; n1=2;n2=4;
                        TheLessTheBetter=a1*sum( abs(Lengths(twoshort)-Bound(1) ).^n1)+ a2*sum( abs(Lengths(OutOfBound)-Bound(2) ).^n2) ;
                        
                        
                    end
                    %--------------
                    if TheLessTheBetter<Saved.Val
                        Saved.Lengths=Lengths; Saved.GlobalCut= GlobalCut;
                        Saved.Val=TheLessTheBetter; Saved.TroubleLeng=LengHasTrouble;
                        nloop;
                        nloop=0;
                        cutmoreLong=0;
                    elseif TheLessTheBetter==Saved.Val || cc==2
                        Saved.Lengths=Lengths; Saved.GlobalCut= GlobalCut;
                        Saved.Val=TheLessTheBetter; Saved.TroubleLeng=LengHasTrouble;
                        nloop= nloop+1 ;
                    else
                        if ~isempty(Saved.GlobalCut) && mod(nloop,loop2Reset)>= 1
                            GlobalCut=Saved.GlobalCut;
                        end
                        nloop= nloop+1 ;
                        
                    end
                    %                   [iStap , nloop]
                    relax=rand(1);
                    %                                             Lengths
                    %                                             Saved.Lengths
                    %                                             if ismember(0,Lengths)
                    %                                                 sdfsf=3
                    %                                             end
                    
                    [IncorrectInterval,LongwrongInt]=SPbreakfindleng2Interval(obj, Saved.Lengths,Interval,Ex,relax,LongInterval,LongGInfo,Bound);
                    MiddleIntervals=intersect(OMI,IncorrectInterval);
                    LongInterval=intersect(OLI,LongwrongInt);
                    
                    if mod(nloop,loop2Reset)==loop2Reset-1
                        %                     GlobalCut=Saved.GlobalCut;
                        GlobalCut=-(int8(DigOneForLong));
                        %                      Saved.GlobalCut=GlobalCut;
                        MiddleIntervals=OMI; HeadAndTail=OHT; LongInterval=OLI;
                        reset=1+reset;
                        saveplot(cc,[1 3])=[cc,1];
                    end
                    if nloop == 1500 && Saved.Val>10   %still has one too long seg
                        cutmoreLong=15;
                    end
                    
                    saveplot(cc,1:2)=[cc,Saved.Val];
                    saveplot(cc,4)=TheLessTheBetter;
                    if TheLessTheBetter<=0  || nloop> 500   %or use BestTargetVal
                        break
                    end
                    nloop;
                end  %end of while loop
                
                Lengths=Saved.Lengths;
                GlobalCut=Saved.GlobalCut;  % [ 1 2 3 ...] indicate which edge to cut, -1 means long edge(no cut)
                GlobalCut(GlobalCut==1)=1:length(find(GlobalCut==1));
                WW=Saved.GlobalCut;
                ActualCutedge=find(WW==1);
                LengOfCutEdge=Ex(ActualCutedge);
                %                  if sum(LengOfCutEdge>7)>0
                HaveToCutOnLeng=sum(LengOfCutEdge>DefLongSegThan);
                HaveToCutOnLeng= [ HaveToCutOnLeng ,  iStap];
                
                
                %                  end
                Apart=round(LengOfCutEdge/2);
                Bpart=LengOfCutEdge-Apart;
                
                OldStp= obj.StapList2{iStap};
                %                 StpStart=OldStp(1,:);
                AddCell=cell(length(ActualCutedge)+1,1);
                
                % GlobalCut
                % ActualCutedge
                % OldStp
                if OldStp(2*ActualCutedge(1),2)> OldStp(2*ActualCutedge(1)-1,2)
                    AddCell{1}=[OldStp(1:2*ActualCutedge(1)-1,:) ;   [OldStp(2*ActualCutedge(1)-1,1),OldStp(2*ActualCutedge(1)-1,2)+Apart(1)-1] ] ;
                else
                    AddCell{1}=[OldStp(1:2*ActualCutedge(1)-1,:) ;   [OldStp(2*ActualCutedge(1)-1,1),OldStp(2*ActualCutedge(1)-1,2)-Apart(1)+1] ]  ;
                end
                
                if OldStp(2*ActualCutedge(end)-1,2)<OldStp(2*ActualCutedge(end),2)
                    AddCell{end}=[ [OldStp(2*ActualCutedge(end),1),OldStp(2*ActualCutedge(end),2)-(Bpart(end)-1)]   ;OldStp(2*ActualCutedge(end):end,:)];
                else
                    AddCell{end}=[ [OldStp(2*ActualCutedge(end),1),OldStp(2*ActualCutedge(end),2)+(Bpart(end)-1)]   ;OldStp(2*ActualCutedge(end):end,:)];
                end
                
                for AddNewCell=2:length(ActualCutedge)
                    Vec1=OldStp(2*ActualCutedge(AddNewCell-1),:) -OldStp(2*ActualCutedge(AddNewCell-1)-1,:);
                    %                                         Vec1=OldStp(2*ActualCutedge(AddNewCell-1)+1,:)-OldStp(2*ActualCutedge(AddNewCell-1)+2,:);
                    
                    Vec1=Vec1/norm(Vec1);
                    BeforeLine=OldStp(2*ActualCutedge(AddNewCell-1),:)-(Bpart(AddNewCell-1)-1)*Vec1;
                    Vec2= OldStp(2*ActualCutedge(AddNewCell),:) -OldStp(2*ActualCutedge(AddNewCell)-1,:);
                    %                                                         Vec2=OldStp(2*ActualCutedge(AddNewCell)-3,:)-OldStp(2*ActualCutedge(AddNewCell)-2,:);
                    
                    Vec2=Vec2/norm(Vec2)  ;
                    AfterLine=OldStp(2*ActualCutedge(AddNewCell)-1,:)+(Apart(AddNewCell)-1)*Vec2;
                    
                    Temp=[BeforeLine; OldStp(2*ActualCutedge(AddNewCell-1):2*ActualCutedge(AddNewCell)-1,:) ;AfterLine] ; % temperory for straight cut case.
                    if length(unique(Temp(:,1)))==1
                        Temp=sortrows(Temp) ;
                    end
                    AddCell{AddNewCell}=Temp;
                end
                
                newcell(NewCell_count:NewCell_count+length(ActualCutedge))=AddCell(1:end);
                NewCell_count=NewCell_count+length(ActualCutedge)+1;
                
            end   %end if for all iStap
            
            OutPutStapList=newcell(1:NewCell_count-1);
            [ Res2 ] = showstapins( OutPutStapList );
            figure(21);cla;  subplot(2,1,1 );cla;
            histogram(Res.InitialLength);  title(strcat('Before Cutting, stap Lengths ,  N= ',num2str(length(Res.InitialLength)) ));xlabel('L (nt)'); ylabel('N');
            subplot(2,1, 2);cla;
            histogram(Res2.InitialLength);
            title(strcat('After Cutting, stap Lengths,  N= ',num2str(length(OutPutStapList)) ));xlabel('L (nt)'); ylabel('N');
            flagisEdgestp=zeros(size(OutPutStapList));
            %              for kstp=1:length(flagisEdgestp)
            %                 CC=OutPutStapList{kstp};
            %
            %
            %                  HeadisOnZ1Surface=  obj.Z1(CC(1,1))== CC(1,2);
            %                  tailisOnZ1Surface=  obj.Z1(CC(end,1))== CC(end,2);
            %                  HeadisOnZ2Surface=  obj.Z2(CC(1,1))== CC(1,2);
            %                  tailisOnZ2Surface=  obj.Z2(CC(end,1))== CC(end,2);
            %                  IsEdgeStp=HeadisOnZ1Surface||HeadisOnZ1Surface||tailisOnZ1Surface||HeadisOnZ2Surface||tailisOnZ2Surface;
            %                  if IsEdgeStp==1;flagisEdgestp(kstp)=1;end
            %              end
            %             subplot(2,2,3); cla;
            %              histogram(Res2.InitialLength(flagisEdgestp==1));
            %              title(  strcat('Edge Stap Lengths,  N= ',num2str( sum(flagisEdgestp))) );xlabel('L (nt)'); ylabel('N');
            %              subplot(2,2,4); cla;
            %               histogram(Res2.InitialLength(flagisEdgestp==0));
            %              title(  strcat('Non-Edge Stap Lengths,  N= ',num2str(sum(flagisEdgestp==0))));xlabel('L (nt)'); ylabel('N');
            obj.StapList3=OutPutStapList;
        end
        
        function extendOverhang(obj,ax  )
            axes(ax) ;
            table_OH=findobj(gcf,'Tag','OHTable') ;
            Info=  table_OH.Data;
            IndEff= cellfun(@isequal,Info(:,3) , num2cell(ones(size(Info(:,3) ) ) ) ) ;
            InfoHide= table_OH.UserData ;   InfoHide=InfoHide(IndEff,:)  ;
            Nick= InfoHide(:,[1 4]);   % Nick= InfoHide(:,[1 4]);
            Nick=Nick(:);
            Info=Info(IndEff,:)  ;%
            ConvertC3= zeros(length(Nick)*2 , 2) ;
            for k=1:length(Nick)
                Arr= str2num(Nick{k} ) ; %#ok<ST2NM>
                [~,Ind]= ismember(  Arr(1:2), obj.RelateTable(:,1:2),'rows') ;
                C3_ind=  obj.RelateTable(Ind,3) ;
                ConvertC3(2*k-1,:) =[C3_ind, Arr(3)] ;
                ConvertC3(2*k,:) =[C3_ind, Arr(4)] ;
            end
            
            HeadAndTail =zeros( 2*length(obj.StapList3),2) ;
            for k=1:  length(obj.StapList3)
                HeadAndTail(2*k-1,:) =  obj.StapList3{k}(1,:) ;
                HeadAndTail(2*k,:) =  obj.StapList3{k}(end,:) ;
            end
            ConvertC3   ;
            [CheckNick,inds] = ismember(ConvertC3, HeadAndTail,'rows') ;
            nCylOri= size(obj.RelateTable,1) ;
            AppRT = zeros(length(Nick),6) ;
            for k=1:length(Nick)
                comefrom =  str2num(Nick{k} ) ;
                SourceXY=  obj.RelateTable( ConvertC3(2*k-1,1) ,6:7) ;
                PotentialXY = [SourceXY+[1,0] ;SourceXY+[-1,0] ; SourceXY+[0,1] ;SourceXY+[0,-1]] ;   %
                PotentialXY=setdiff(PotentialXY ,obj.RelateTable(:,6:7),'rows' ) ;
                Bundle=obj.containBundle{comefrom(1)} ;
                if size( PotentialXY,1) ~=1
                    BaseL= comefrom(3)- obj.containBundle{comefrom(1)}.Zbase1(comefrom(2))+1   ;
                    BaseR= comefrom(4)- obj.containBundle{comefrom(1)}.Zbase1(comefrom(2))+1   ;
                    
                    HCoor1=HelixXYZGStapNoPMLocal(obj.containBundle{comefrom(1)},2,1)  ;
                    HCoor2=HelixXYZGStapNoPMLocal(obj.containBundle{comefrom(1)},0,1)  ;
                    %                  TMatrix= obj.containBundle{comefrom(1)}.TransformMatrix2 ;
                    %                  XYZ_3D=HCoor{comefrom(2)}(Base,:)  ;
                    VecL=HCoor1{comefrom(2)}(BaseL,:) - HCoor2{comefrom(2)}(BaseL,:);
                    VecR=HCoor1{comefrom(2)}(BaseR,:) - HCoor2{comefrom(2)}(BaseR,:);
                    Vec=VecL+VecR;
                    
                    Vec=Vec(1:2); Vec=Vec/norm(Vec) ;
                    TargetXY =  SourceXY+ Vec ;
                    d=  PotentialXY - ones(size( PotentialXY,1),1)* TargetXY ;    d=d(:,1).^2 +d(:,2).^2 ;
                    Select = find(d==min(d)) ;
                    PotentialXY=PotentialXY( Select,:) ;
                    
                end
                
                AddRow= [comefrom PotentialXY] ;
                AppRT(k,:)=AddRow ;
            end
            AppRT=union(AppRT,AppRT,'rows','stable') ;  % user given overhangs to new cylinder
            OverhangToNewCylinder =  AppRT ;   % user given overhangs to new cylinder
            
            AppRT=AppRT(:,[1,5:6]) ;[AppRT,ia,ib]=union(AppRT,AppRT,'rows','stable') ;
            cut=AppRT(:,1) ;  AppRT=AppRT(:,2:3) ;
            AppRT=[ -1*ones(size(AppRT)) , AppRT] ;
            maxOdd= obj.RelateTable(:,4);  maxOdd=maxOdd(mod(maxOdd,2)==1); maxOdd=max(maxOdd) ;
            maxEvv= obj.RelateTable(:,4);  maxEvv=maxEvv(mod(maxEvv,2)==0); maxEvv=max(maxEvv) ;
            for k=1:size(AppRT,1)
                AppRT(k,2) = max(obj.RelateTable(:,5))+k ;
                %                 AppRT(k,2) = size(obj.RelateTable,1)+k ;
                if  mod(AppRT(k,3)+AppRT(k,4) ,2)==0
                    AppEvn=AppRT(:,1); AppEvn=AppEvn(mod(AppEvn,2)==0) ;
                    CurrentEvenMax = max([maxEvv;AppEvn]) ;
                    AppRT(k,1)=CurrentEvenMax+2 ;
                else
                    AppOdd=AppRT(:,1); AppOdd=AppOdd(mod(AppOdd,2)==1) ;
                    CurrentOddMax = max([maxOdd;AppOdd]) ;
                    AppRT(k,1)=CurrentOddMax+2 ;
                end
            end
            AppRT= [ cut,-1*ones(size(AppRT,1) , 2) , AppRT ] ;   % columns 1-3 hopefully don't matter
            
            OverhangToNewCylinder ;   % user given overhangs to new cylinder
            Info;
            InfoHide ;
            
            %---------Staple assembly, staple crossovers, Dec 22 2020 , 
% % % % % %-----------    Jan 19 2021
            Extras = sum( or(contains( Info(:,8) , 'Double crossover'), contains( Info(:,8) , 'Double overhang1')) ) ;
            PairMappForOverhangCases =zeros(2*(size(Info,1)+Extras) , 1) ;    % due to double Xovers
            currentCount = 1;
            for c= 1:size(Info,1)
                if  contains( Info(c,8) , 'Double crossover')||contains( Info(c,8) , 'Double overhang1')||contains( Info(c,8) , 'Double overhang2')                                                 
                    % contains( Info(c,8) , 'Double crossover')
                PairMappForOverhangCases(currentCount:currentCount+3) =  c ;
                currentCount=currentCount+4;
                else
                PairMappForOverhangCases(currentCount:currentCount+1) =  c ;
                currentCount=currentCount+2;
                end
            end
            
            DecidePrims=zeros(2*(size(Info,1)+Extras) , 3) ;  % root(C5 + Base) + lengths(4) 
%             DecidePrims=zeros(2*size(Info,1) , 3) ;  % root(C5 + Base) + lengths(4)
            
            
%             for extraExtend = 1 :size(Info,1) ;
%             
%             end
                
            ck=1 ; Saveck =zeros(1000,2); cSk =1 ;
            for k=1:size(Info,1)
                comefrom =  str2num(InfoHide{k,1} ) ;   %
                [~,C3Cylinder]= ismember(  comefrom(1:2), obj.RelateTable(:,1:2),'rows') ;
                if xor(strcmp(Info{k,2},'3'''),  mod(obj.RelateTable(C3Cylinder,4) ,2) == 0)   %
                    root= [obj.RelateTable(C3Cylinder,5) ,comefrom(3) ] ;  Dir = - Info{k,6}+1 ;
                else
                    root= [obj.RelateTable(C3Cylinder,5) ,comefrom(4) ] ;  Dir =  Info{k,6}-1 ;
                end
                if Info{k,6}==0 ; Dir=0;end   % for single extension
                
                DecidePrims(ck,:) =[root ,Dir ] ; ck=ck+1;
                %  ------------
                comefrom =  str2num(InfoHide{k,4} ) ;   %                C3Cylinder =
                [~,C3Cylinder]= ismember(  comefrom(1:2), obj.RelateTable(:,1:2),'rows') ;
                if xor(strcmp(Info{k,5},'3'''),  mod(obj.RelateTable(C3Cylinder,4) ,2) == 0)   %
                    root= [obj.RelateTable(C3Cylinder,5) ,comefrom(3) ] ; Dir = - Info{k,7}+1 ;
                else
                    root= [obj.RelateTable(C3Cylinder,5) ,comefrom(4) ] ; Dir = Info{k,7}-1 ;
                end
                if Info{k,7}==0 ; Dir=0;end % for single extension
                DecidePrims(ck,:) =[root , Dir ] ; ck=ck+1;
                
                if  contains( Info(k,8), 'Double crossover') ||  contains( Info(k,8), 'Single crossover') || contains( Info(k,8), 'Double overhang1')
                    Saveck(cSk,: ) =[ck-2 , ck-1 ] ;  cSk=cSk+1 ;
                end
                %-------------
                
                if  contains( Info(k,8), 'Double crossover') || contains( Info(k,8), 'Double overhang1') || contains( Info(k,8), 'Double overhang2')
                    comefrom =  str2num(InfoHide{k,1} ) ;   %
                    [~,C3Cylinder]= ismember(  comefrom(1:2), obj.RelateTable(:,1:2),'rows') ;
                    if ~xor(strcmp(Info{k,2},'3'''),  mod(obj.RelateTable(C3Cylinder,4) ,2) == 0)   %
                        root= [obj.RelateTable(C3Cylinder,5) ,comefrom(3) ] ;  Dir = - Info{k,6}+1 ;
                    else
                        root= [obj.RelateTable(C3Cylinder,5) ,comefrom(4) ] ;  Dir =  Info{k,6}-1 ;
                    end
                    if Info{k,6}==0 ; Dir=0;end   % for single extension
                    
                    DecidePrims(ck,:) =[root ,Dir ] ; ck=ck+1;
                    %  ------------
                    comefrom =  str2num(InfoHide{k,4} ) ;   %                C3Cylinder =
                    [~,C3Cylinder]= ismember(  comefrom(1:2), obj.RelateTable(:,1:2),'rows') ;
                    if ~xor(strcmp(Info{k,5},'3'''),  mod(obj.RelateTable(C3Cylinder,4) ,2) == 0)   %
                        root= [obj.RelateTable(C3Cylinder,5) ,comefrom(3) ] ; Dir = - Info{k,7}+1 ;
                    else
                        root= [obj.RelateTable(C3Cylinder,5) ,comefrom(4) ] ; Dir = Info{k,7}-1 ;
                    end
                    if Info{k,7}==0 ; Dir=0;end % for single extension
                    DecidePrims(ck,:) =[root , Dir ] ; ck=ck+1;
                    
                    DecidePrims([ck-1 , ck-3] ,:) =    DecidePrims([ck-3 , ck-1] ,:) ;
                    %-------------
                    Saveck(cSk,: ) =[ck-2 , ck-1 ] ;  cSk=cSk+1 ;
                end
                
            
                
            end
            Saveck=Saveck(1:cSk-1 ,: )  ;% use to connect single/double crossover, row index in DecidePrims
            % DecidePrims column 3 == means 0-nt extension            
            count=0;
            DecidePrims ;
            ExtendScaf = zeros(2*size(DecidePrims,1) ,2 ) ; nn=1;
            TargetRCylandDecidePrims=[-1*ones(size(DecidePrims,1),1) , DecidePrims ] ;
            % later used for linking pseudo scaffold section
            for stpi= 1:length( obj.StapList3)
                stpR= obj.StapList3{stpi} ;
                %                 AA=ismember(   stpR(1,:), DecidePrims(:,1:2) , 'rows') ;
                %                 BB=ismember(   stpR(end,:), DecidePrims(:,1:2) , 'rows')  ;
                %                 [stpi,AA,BB,and(AA,BB) ]
                if ismember(   stpR(1,:), DecidePrims(:,1:2) , 'rows')
                    [~, whichOverhangs] =   ismember(   stpR(1,:), DecidePrims(:,1:2) , 'rows') ;
                    %                     RootXY=obj.RelateTable( DecidePrims(whichOverhangs,1),6:7) ;% will have error for discontinuous cylinders
                    indsC55 = find(obj.RelateTable(:,5)==DecidePrims(whichOverhangs,1));indsC55=indsC55(1) ;
                    RootXY=obj.RelateTable( indsC55,6:7) ;
                    d= OverhangToNewCylinder(:,5:6)- ones(size(OverhangToNewCylinder,1),1)*RootXY;
                    d=d(:,1).^2 +d(:,2).^2 ;    d=d==1;
                    check2= DecidePrims(whichOverhangs,2)==OverhangToNewCylinder(:,3) ;
                    check3= DecidePrims(whichOverhangs,2)==OverhangToNewCylinder(:,4) ;
                    targetCylC5ind=find(and(or(check2,check3),d))   ;
                    TargetXY= OverhangToNewCylinder(targetCylC5ind,5:6) ;
                    [~,TargetC5ind ]=   ismember(   TargetXY, AppRT(:,6:7) , 'rows');
                    TargetC5=AppRT(TargetC5ind,5);
                    [~,qqind] =ismember(   stpR(1,:), DecidePrims(:,1:2) , 'rows') ;
                    TargetRCylandDecidePrims(qqind,1)=TargetC5 ;
                    if DecidePrims(whichOverhangs,3)~=0
                        NewR = [[TargetC5, stpR(1,2)+DecidePrims(whichOverhangs,3)];[TargetC5, stpR(1,2)];stpR    ];
                    else % single extension
                        NewR=   stpR;
                    end
                    obj.StapList3{stpi}= NewR;
                    ExtendScaf(2*nn-1:2*nn, : ) =[[TargetC5, stpR(1,2)];[TargetC5, stpR(1,2)+DecidePrims(whichOverhangs,3)]   ] ;
                    nn=nn+1 ;
                end
            end
            for stpi= 1:length( obj.StapList3)
                stpR= obj.StapList3{stpi} ;
                if  ismember(   stpR(end,:), DecidePrims(:,1:2) , 'rows')
                    [~, whichOverhangs] =   ismember(   stpR(end,:), DecidePrims(:,1:2) , 'rows') ;
                    %                     RootXY=obj.RelateTable( DecidePrims(whichOverhangs,1),6:7) ;  % will have error for discontinuous cylinders
                    indsC55 = find(obj.RelateTable(:,5)==DecidePrims(whichOverhangs,1));indsC55=indsC55(1) ;
                    RootXY=obj.RelateTable( indsC55,6:7) ;
                    d= OverhangToNewCylinder(:,5:6)- ones(size(OverhangToNewCylinder,1),1)*RootXY;
                    d=d(:,1).^2 +d(:,2).^2 ;    d=d==1;
                    check2= DecidePrims(whichOverhangs,2)==OverhangToNewCylinder(:,3) ;
                    check3= DecidePrims(whichOverhangs,2)==OverhangToNewCylinder(:,4) ;
                    targetCylC5ind=find(and(or(check2,check3),d))   ;
                    TargetXY= OverhangToNewCylinder(targetCylC5ind,5:6) ;
                    [~,TargetC5ind ]=   ismember(   TargetXY, AppRT(:,6:7) , 'rows');
                    TargetC5=AppRT(TargetC5ind,5);
                    [~,qqind] =ismember(    stpR(end,:), DecidePrims(:,1:2) , 'rows');
                    TargetRCylandDecidePrims(qqind,1)=TargetC5 ;
                    if DecidePrims(whichOverhangs,3)~=0
                        NewR = [stpR ; [TargetC5, stpR(end,2)] ; [TargetC5, stpR(end,2)+DecidePrims(whichOverhangs,3)]];
                    else   % single extension
                        NewR = stpR ;
                    end
                    obj.StapList3{stpi}= NewR;
                    ExtendScaf(2*nn-1:2*nn, : ) =[[TargetC5, stpR(end,2)+DecidePrims(whichOverhangs,3)];[TargetC5, stpR(end,2)]   ] ;
                    nn=nn+1 ;
                end
            end
            
            ExtendScaf;
            %------------connect nick or overhang for single/double Xover,
            %Dec 22 2020
            DecidePrims ; TargetRCylandDecidePrims;
            for k = 1: size(Saveck ,1)
                if  DecidePrims( Saveck(k,1) ,3) ==0
                    ConnectEndA  =     DecidePrims( Saveck(k,1) ,1:2) ;
                else
%                     TargetRCy = TargetRCylandDecidePrims(TargetRCylandDecidePrims(:,2)== DecidePrims( Saveck(k,1) ,1),1) ;
%                     TargetRCy=unique(TargetRCy)  ;                   
                    TargetRCy = TargetRCylandDecidePrims(Saveck(k,1),1 ) ;
                    
                    ConnectEndA  = [TargetRCy,  DecidePrims( Saveck(k,1) ,2)+ DecidePrims( Saveck(k,1) ,3) ] ;
                end
                if  DecidePrims( Saveck(k,2) ,3) ==0
                    ConnectEndB  =     DecidePrims( Saveck(k,2) ,1:2) ;
                else
%                     TargetRCy = TargetRCylandDecidePrims(TargetRCylandDecidePrims(:,2)== DecidePrims( Saveck(k,2) ,1), 1) ;
%                     TargetRCy=unique(TargetRCy) ;
                    
                    TargetRCy = TargetRCylandDecidePrims(Saveck(k,2),1 ) ;
                    ConnectEndB  = [TargetRCy,  DecidePrims( Saveck(k,2) ,2)+ DecidePrims( Saveck(k,2) ,3) ] ;
                end
                
                %                 StpA =[];
%                 length(obj.StapList3)
                for stpj= 1: length(obj.StapList3)
                    if ismember( ConnectEndA ,obj.StapList3{stpj} ,'rows')
                        StpA=stpj;
                        [~, StpA_ind]  = ismember( ConnectEndA ,obj.StapList3{stpj} ,'rows') ;
                    end
                    if ismember( ConnectEndB ,obj.StapList3{stpj},'rows' )
                        StpB=stpj;
                        [~, StpB_ind]  = ismember( ConnectEndB ,obj.StapList3{stpj} ,'rows') ;
                    end
                end
                
                if StpB_ind==1 && StpA_ind~=1
                obj.StapList3{StpA} = [obj.StapList3{StpA} ; obj.StapList3{StpB}] ;
                obj.StapList3{StpB} =[];
                elseif StpA_ind==1 && StpB_ind~=1
                obj.StapList3{StpA} = [obj.StapList3{StpB} ; obj.StapList3{StpA}] ;
                obj.StapList3{StpB} =[];
                end
                [~,bb] = cellfun( @size ,obj.StapList3 ) ;
                obj.StapList3  =obj.StapList3(bb~=0) ;
            end
            %-----------
            
            obj.RelateTable=[obj.RelateTableOrr ; AppRT] ;  %update RTable
            
            [~,index]=sortrows(obj.RelateTable,5);    %update RVec
            Vec=obj.RelateTable(:,4);
            TTRelateVec=Vec(index);
            WQER =union(TTRelateVec,[],'stable');
            obj.RelateVec=WQER';
            %--------------
            
            obj.ConvertScafG;
            obj.ConvertScafSQ;      % Get properties:ScafdigitSQ
            obj.ConvertScafHC;      % Get properties:ScafdigitHC
            %             dfgdfg=4 ; ExtendScaf ;
            %---------------update  %obj.ScafdigitHC    obj.ScafdigitSQ
            for k = 1 :  nn-1
                FromTo = ExtendScaf(2*k-1 :2*k,:) ;
                DDor= FromTo(2,2)-FromTo(1,2) ;  DDor=DDor/norm(DDor) ;
                Arr= (FromTo(1,2):DDor:FromTo(2,2))' ; Arr=[ FromTo(1,1)*ones(length(Arr),1)  ,Arr ];
                Arr2 = [circshift(Arr,1) , circshift(Arr,-1) ] ;
                Arr2(1,1:2)=-1 ; Arr2(end, 3:4)=-1 ;
                if  DDor==1
                    obj.ScafdigitSQ{FromTo(1,1)}(FromTo(1,2):FromTo(2,2),: ) =Arr2 ;
                    obj.ScafdigitHC{FromTo(1,1)}(FromTo(1,2):FromTo(2,2),: ) =Arr2 ;
                else
                    obj.ScafdigitSQ{FromTo(1,1)}(FromTo(2,2):1:FromTo(1,2),: ) =flip(Arr2 );
                    obj.ScafdigitHC{FromTo(1,1)}(FromTo(2,2):1:FromTo(1,2),: ) =flip(Arr2 );
                end
            end
            %----------------------||||| discuss case
            %             dfsaf = 23 ;  TargetRCylandDecidePrims;
            PairMappForOverhangCases ;
            ForEasy =TargetRCylandDecidePrims ;
            ConnScaf2= cell(size(PairMappForOverhangCases ,1)/2,1) ; extra_by_doubleXover_n = 0 ;
            for k= 1:size(Info ,1)
                %                 k
                %                 switch Info{k,8}
                %                     case 'r_r'
                %                          MM =  [ForEasy(2*k-1:2*k ,1 ) , ForEasy(2*k-1:2*k ,3 )];
                %
                %                     case 'e_e'
                %                          MM =  [ForEasy(2*k-1:2*k ,1 ) , ForEasy(2*k-1:2*k ,3 )+ForEasy(2*k-1:2*k ,4 )];
                %
                %                     case 'Ar_Be'
                %                           MM =  [ForEasy(2*k-1:2*k ,1 ) , [ForEasy(2*k-1 ,3 );  ForEasy(2*k ,3 )+ForEasy(2*k ,4 )    ]   ];
                %
                %                     case 'Ae_Br'
                %                           MM =  [ForEasy(2*k-1:2*k ,1 ) , [ForEasy(2*k-1 ,3 )+ForEasy(2*k-1 ,4 );  ForEasy(2*k ,3 )    ]   ];
                %
                %                 end
                if strcmp(Info{k,8},'Connected')  &&  ~strcmp(Info{k,2},Info{k,5})
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , ForEasy(2*k-1:2*k ,3 )];
                elseif  strcmp(Info{k,8},'Free end')  &&  ~strcmp(Info{k,2},Info{k,5})
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , ForEasy(2*k-1:2*k ,3 )+ForEasy(2*k-1:2*k ,4 )];
                elseif strcmp(Info{k,8},'Connected')  &&  strcmp(Info{k,2},Info{k,5})
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , [ForEasy(2*k-1 ,3 );  ForEasy(2*k ,3 )+ForEasy(2*k ,4 )    ]   ];
                elseif strcmp(Info{k,8},'Free end')  &&  strcmp(Info{k,2},Info{k,5})
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , [ForEasy(2*k-1 ,3 )+ForEasy(2*k-1 ,4 );  ForEasy(2*k ,3 )    ]   ];
                    
                    %                 elseif  strcmp(Info{k,8},'Double overhang2')
                    
                    
                elseif strcmp(Info{k,8},'Double crossover') ||  strcmp(Info{k,8},'Double overhang1') ||strcmp(Info{k,8},'Double overhang2') % similar to case 2
                    % %                     MM =  [ForEasy(2*k-1:2*k ,1 ) , [ForEasy(2*k-1 ,3 )+ForEasy(2*k-1 ,4 );  ForEasy(2*k ,3 )    ]   ];
                    IndsMapping  = find(PairMappForOverhangCases == k) ;
                    IndsMapping = [IndsMapping(3:4 );IndsMapping(1:2 )];
                    
                    
                    if ~strcmp(Info{k,8},'Double overhang2')
                        MM =  [ForEasy(IndsMapping(1:2) ,1 ) , ForEasy(IndsMapping(1:2) ,3 )+ForEasy(IndsMapping(1:2) ,4 )];
                        MM2 =[ForEasy(IndsMapping(3:4) ,1 ) , ForEasy(IndsMapping(3:4) ,3 )+ForEasy(IndsMapping(3:4) ,4 )];
                    else
                        MM =  [ForEasy(IndsMapping([2 3]) ,1 ) , ForEasy(IndsMapping([2 3]) ,3 )+0*ForEasy(IndsMapping([2 3]) ,4 )];
                        MM2 =[ForEasy(IndsMapping([4 1]) ,1 ) , ForEasy(IndsMapping([4 1]) ,3 )+ForEasy(IndsMapping([4 1]) ,4 )];
                        %                          MM =  [ForEasy(IndsMapping([2 3]) ,1 ) , ForEasy(IndsMapping([2 3]) ,3 )+ForEasy(IndsMapping([2 3]) ,4 )];
                        %                          MM2 =[ForEasy(IndsMapping([4 1]) ,1 ) , ForEasy(IndsMapping([4 1]) ,3 )+ForEasy(IndsMapping([4 1]) ,4 )];
                        
                        
                    end
                    
                    
                    MM2=flip(MM2)   ; % first row 5' end , second row 3' end ,otherwise cadnano not show
                    if MM2(1,1)~=-1
                        FirstRow= obj.ScafdigitSQ{MM2(1,1)}(MM2(1,2),:);
                        Inds=find(FirstRow==-1) ;
                        FirstRow(Inds(1)) = MM2(2,1) ;   FirstRow(Inds(2)) = MM2(2,2) ;
                        obj.ScafdigitSQ{MM2(1,1)}(MM2(1,2),:)=FirstRow ;
                        obj.ScafdigitHC{MM2(1,1)}(MM2(1,2),:)=FirstRow ;
                    end
                    
                    if MM2(2,1)~= -1
                        SecondRow= obj.ScafdigitSQ{MM2(2,1)}(MM2(2,2),:) ;
                        Inds2=find(SecondRow==-1)  ;
                        SecondRow(Inds2(1)) = MM2(1,1) ;   SecondRow(Inds2(2)) = MM2(1,2) ;
                        obj.ScafdigitSQ{MM2(2,1)}(MM2(2,2),:)= SecondRow ;
                        obj.ScafdigitHC{MM2(2,1)}(MM2(2,2),:)= SecondRow ;
                    end
                    extra_by_doubleXover_n=extra_by_doubleXover_n +1 ;
                    ConnScaf2{k+extra_by_doubleXover_n-1} =MM2 ;
                    %                     MM = [-1 , -1;-1 -1 ];
                elseif strcmp(Info{k,8},'Single crossover')
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , ForEasy(2*k-1:2*k ,3 )+ForEasy(2*k-1:2*k ,4 )];
                    %                      MM = [-1 , -1;-1 -1 ];
                end
                
                MM=flip(MM);
                % first row 5' end , second row 3' end ,otherwise cadnano not show
                if MM(1,1)~=-1
                    FirstRow= obj.ScafdigitSQ{MM(1,1)}(MM(1,2),:);
                    Inds=find(FirstRow==-1);

                    FirstRow(Inds(1)) = MM(2,1) ;   FirstRow(Inds(2)) = MM(2,2) ;
                    obj.ScafdigitSQ{MM(1,1)}(MM(1,2),:)=FirstRow ;
                    obj.ScafdigitHC{MM(1,1)}(MM(1,2),:)=FirstRow ;
                end
                
                if MM(2,1)~= -1
                    SecondRow= obj.ScafdigitSQ{MM(2,1)}(MM(2,2),:) ;
                    Inds2=find(SecondRow==-1)  ;
                    SecondRow(Inds2(1)) = MM(1,1) ;   SecondRow(Inds2(2)) = MM(1,2) ;
                    obj.ScafdigitSQ{MM(2,1)}(MM(2,2),:)= SecondRow ;
                    obj.ScafdigitHC{MM(2,1)}(MM(2,2),:)= SecondRow ;
                end
                
                ConnScaf2{k+extra_by_doubleXover_n} =MM ;
                %                 MM
            end
                %                 k
            obj.ClosingStrand.ExtendScaf=ExtendScaf;
            obj.ClosingStrand.ConnScaf2=ConnScaf2;
            obj.ClosingStrand.RootAndExtC5 = TargetRCylandDecidePrims ;
        end  % end of  extendOverhang
        
        function LoadJsonRecoverClosingStrand(obj  )
            %-------only use for feeding json file with OH to generate
            %closing strands on json files. Do not update other properties
            table_OH=findobj(gcf,'Tag','OHTable') ;
            Info=  table_OH.Data;
            IndEff= cellfun(@isequal,Info(:,3) , num2cell(ones(size(Info(:,3) ) ) ) ) ;
            InfoHide= table_OH.UserData ;   InfoHide=InfoHide(IndEff,:)  ;
            Nick= InfoHide(:,[1 4]);   % Nick= InfoHide(:,[1 4]);
            Nick=Nick(:);
            Info=Info(IndEff,:)  ;%
            ConvertC3= zeros(length(Nick)*2 , 2) ;
            for k=1:length(Nick)
                Arr= str2num(Nick{k} ) ; %#ok<ST2NM>
                [~,Ind]= ismember(  Arr(1:2), obj.RelateTableOrr(:,1:2),'rows') ;
                C3_ind=  obj.RelateTableOrr(Ind,3) ;
                ConvertC3(2*k-1,:) =[C3_ind, Arr(3)] ;
                ConvertC3(2*k,:) =[C3_ind, Arr(4)] ;
            end
            
            HeadAndTail =zeros( 2*length(obj.StapList3),2) ;
            for k=1:  length(obj.StapList3)
                HeadAndTail(2*k-1,:) =  obj.StapList3{k}(1,:) ;
                HeadAndTail(2*k,:) =  obj.StapList3{k}(end,:) ;
            end
            [CheckNick,inds] = ismember(ConvertC3, HeadAndTail,'rows') ;
            nCylOri= size(obj.RelateTableOrr,1) ;
            AppRT = zeros(length(Nick),6) ;
            for k=1:length(Nick)
                
                comefrom =  str2num(Nick{k} ) ;
                SourceXY=  obj.RelateTableOrr( ConvertC3(2*k-1,1) ,6:7) ;
                PotentialXY = [SourceXY+[1,0] ;SourceXY+[-1,0] ; SourceXY+[0,1] ;SourceXY+[0,-1]] ;   %
                PotentialXY=setdiff(PotentialXY ,obj.RelateTableOrr(:,6:7),'rows' ) ;
                Bundle=obj.containBundle{comefrom(1)} ;
                if size( PotentialXY,1) ~=1
                    BaseL= comefrom(3)- obj.containBundle{comefrom(1)}.Zbase1(comefrom(2))+1   ;
                    BaseR= comefrom(4)- obj.containBundle{comefrom(1)}.Zbase1(comefrom(2))+1   ;
                    
                    HCoor1=HelixXYZGStapNoPMLocal(obj.containBundle{comefrom(1)},2,1)  ;
                    HCoor2=HelixXYZGStapNoPMLocal(obj.containBundle{comefrom(1)},0,1)  ;
                    %                  TMatrix= obj.containBundle{comefrom(1)}.TransformMatrix2 ;
                    %                  XYZ_3D=HCoor{comefrom(2)}(Base,:)  ;
                    VecL=HCoor1{comefrom(2)}(BaseL,:) - HCoor2{comefrom(2)}(BaseL,:);
                    VecR=HCoor1{comefrom(2)}(BaseR,:) - HCoor2{comefrom(2)}(BaseR,:);
                    Vec=VecL+VecR;
                    
                    Vec=Vec(1:2); Vec=Vec/norm(Vec) ;
                    TargetXY =  SourceXY+ Vec ;
                    d=  PotentialXY - ones(size( PotentialXY,1),1)* TargetXY ;    d=d(:,1).^2 +d(:,2).^2 ;
                    Select = find(d==min(d)) ;
                    PotentialXY=PotentialXY( Select,:) ;
                end
                AddRow= [comefrom PotentialXY] ;
                AppRT(k,:)=AddRow ;
            end
            AppRT=union(AppRT,AppRT,'rows','stable') ;  % user given overhangs to new cylinder
            OverhangToNewCylinder = AppRT;   % user given overhangs to new cylinder
            
            AppRT=AppRT(:,[1,5:6]) ;[AppRT,ia,ib]=union(AppRT,AppRT,'rows','stable') ;
            cut=AppRT(:,1) ;  AppRT=AppRT(:,2:3) ;
            AppRT=[ -1*ones(size(AppRT)) , AppRT] ;
            maxOdd= obj.RelateTableOrr(:,4);  maxOdd=maxOdd(mod(maxOdd,2)==1); maxOdd=max(maxOdd) ;
            maxEvv= obj.RelateTableOrr(:,4);  maxEvv=maxEvv(mod(maxEvv,2)==0); maxEvv=max(maxEvv) ;
            for k=1:size(AppRT,1)
                AppRT(k,2) = max(obj.RelateTableOrr(:,5))+k ;
                %                 AppRT(k,2) = size(obj.RelateTable,1)+k ;
                if  mod(AppRT(k,3)+AppRT(k,4) ,2)==0
                    AppEvn=AppRT(:,1); AppEvn=AppEvn(mod(AppEvn,2)==0) ;
                    CurrentEvenMax = max([maxEvv;AppEvn]) ;
                    AppRT(k,1)=CurrentEvenMax+2 ;
                else
                    AppOdd=AppRT(:,1); AppOdd=AppOdd(mod(AppOdd,2)==1) ;
                    CurrentOddMax = max([maxOdd;AppOdd]) ;
                    AppRT(k,1)=CurrentOddMax+2 ;
                end
            end
            AppRT= [ cut,-1*ones(size(AppRT,1) , 2) , AppRT ] ;   % columns 1-3 hopefully don't matter
            DecidePrims=zeros(2* size(Info,1) , 3) ;  % root(C5 + Base) + lengths(4)
            for k=1:size(Info,1)
                comefrom =  str2num(InfoHide{k,1} ) ;   %
                [~,C3Cylinder]= ismember(  comefrom(1:2), obj.RelateTableOrr(:,1:2),'rows') ;
                if xor(strcmp(Info{k,2},'3'''),  mod(obj.RelateTableOrr(C3Cylinder,4) ,2) == 0)   %
                    root= [obj.RelateTableOrr(C3Cylinder,5) ,comefrom(3) ] ;  Dir = - Info{k,6}+1 ;
                else
                    root= [obj.RelateTableOrr(C3Cylinder,5) ,comefrom(4) ] ;  Dir =  Info{k,6}-1 ;
                end
                DecidePrims(2*k-1,:) =[root ,Dir ] ;
                %  ------------
                comefrom =  str2num(InfoHide{k,4} ) ;   %                C3Cylinder =
                [~,C3Cylinder]= ismember(  comefrom(1:2), obj.RelateTableOrr(:,1:2),'rows') ;
                if xor(strcmp(Info{k,5},'3'''),  mod(obj.RelateTableOrr(C3Cylinder,4) ,2) == 0)   %
                    root= [obj.RelateTableOrr(C3Cylinder,5) ,comefrom(3) ] ; Dir = - Info{k,7}+1 ;
                else
                    root= [obj.RelateTableOrr(C3Cylinder,5) ,comefrom(4) ] ; Dir = Info{k,7}-1 ;
                end
                DecidePrims(2*k,:) =[root , Dir ] ;
            end
            count=0;
            ExtendScaf = zeros(2*size(DecidePrims,1) ,2 ) ; nn=1;
            TargetRCylandDecidePrims=[-1*ones(size(DecidePrims,1),1) , DecidePrims ] ;
            % later used for linking pseudo scaffold section
            for stpi= 1:length( obj.StapList3)
                stpR= obj.StapList3{stpi} ;
                OverHangCyl = ismember(stpR(:,1),AppRT(:,5) ) ;
                stpR=stpR(~OverHangCyl , :)  ;  % ignore overhange section
                
                if ismember(   stpR(1,:), DecidePrims(:,1:2) , 'rows')
                    [~, whichOverhangs] =   ismember(   stpR(1,:), DecidePrims(:,1:2) , 'rows') ;
                    %                     RootXY=obj.RelateTable( DecidePrims(whichOverhangs,1),6:7) ;% will have error for discontinuous cylinders
                    indsC55 = find(obj.RelateTable(:,5)==DecidePrims(whichOverhangs,1));indsC55=indsC55(1) ;
                    RootXY=obj.RelateTable( indsC55,6:7) ;
                    d= OverhangToNewCylinder(:,5:6)- ones(size(OverhangToNewCylinder,1),1)*RootXY;
                    d=d(:,1).^2 +d(:,2).^2 ;    d=d==1;
                    check2= DecidePrims(whichOverhangs,2)==OverhangToNewCylinder(:,3) ;
                    check3= DecidePrims(whichOverhangs,2)==OverhangToNewCylinder(:,4) ;
                    targetCylC5ind=find(and(or(check2,check3),d))   ;
                    TargetXY= OverhangToNewCylinder(targetCylC5ind,5:6) ;
                    [~,TargetC5ind ]=   ismember(   TargetXY, AppRT(:,6:7) , 'rows');
                    TargetC5=AppRT(TargetC5ind,5);
                    [~,qqind] =ismember(   stpR(1,:), DecidePrims(:,1:2) , 'rows') ;
                    TargetRCylandDecidePrims(qqind,1)=TargetC5 ;
                    NewR = [[TargetC5, stpR(1,2)+DecidePrims(whichOverhangs,3)];[TargetC5, stpR(1,2)];stpR    ];
                    %                     obj.StapList3{stpi}= NewR;
                    ExtendScaf(2*nn-1:2*nn, : ) =[[TargetC5, stpR(1,2)];[TargetC5, stpR(1,2)+DecidePrims(whichOverhangs,3)]   ] ;
                    nn=nn+1 ;
                end
            end
            for stpi= 1:length( obj.StapList3)
                stpR= obj.StapList3{stpi} ;
                OverHangCyl = ismember(stpR(:,1),AppRT(:,5) ) ;
                stpR=stpR(~OverHangCyl , :)  ;  % ignore overhange section
                if  ismember(   stpR(end,:), DecidePrims(:,1:2) , 'rows')
                    [~, whichOverhangs] =   ismember(   stpR(end,:), DecidePrims(:,1:2) , 'rows') ;
                    %                     RootXY=obj.RelateTable( DecidePrims(whichOverhangs,1),6:7) ;  % will have error for discontinuous cylinders
                    indsC55 = find(obj.RelateTable(:,5)==DecidePrims(whichOverhangs,1));indsC55=indsC55(1) ;
                    RootXY=obj.RelateTable( indsC55,6:7) ;
                    d= OverhangToNewCylinder(:,5:6)- ones(size(OverhangToNewCylinder,1),1)*RootXY;
                    d=d(:,1).^2 +d(:,2).^2 ;    d=d==1;
                    check2= DecidePrims(whichOverhangs,2)==OverhangToNewCylinder(:,3) ;
                    check3= DecidePrims(whichOverhangs,2)==OverhangToNewCylinder(:,4) ;
                    targetCylC5ind=find(and(or(check2,check3),d))   ;
                    TargetXY= OverhangToNewCylinder(targetCylC5ind,5:6) ;
                    [~,TargetC5ind ]=   ismember(   TargetXY, AppRT(:,6:7) , 'rows');
                    TargetC5=AppRT(TargetC5ind,5);
                    [~,qqind] =ismember(    stpR(end,:), DecidePrims(:,1:2) , 'rows');
                    TargetRCylandDecidePrims(qqind,1)=TargetC5 ;
                    NewR = [stpR ; [TargetC5, stpR(end,2)] ; [TargetC5, stpR(end,2)+DecidePrims(whichOverhangs,3)]];
                    %                     obj.StapList3{stpi}= NewR;
                    ExtendScaf(2*nn-1:2*nn, : ) =[[TargetC5, stpR(end,2)+DecidePrims(whichOverhangs,3)];[TargetC5, stpR(end,2)]   ] ;
                    nn=nn+1 ;
                end
            end
            %---------------
            %             obj.RelateTable=[obj.RelateTableOrr ; AppRT] ;  %update RTable
            [~,index]=sortrows(obj.RelateTable,5);    %update RVec
            Vec=obj.RelateTable(:,4);
            TTRelateVec=Vec(index);
            WQER =union(TTRelateVec,[],'stable');
            %--------------
            obj.ConvertScafG;
            obj.ConvertScafSQ;      % Get properties:ScafdigitSQ
            obj.ConvertScafHC;      % Get properties:ScafdigitHC
            %---------------update  %obj.ScafdigitHC    obj.ScafdigitSQ
            for k = 1 :  nn-1
                FromTo = ExtendScaf(2*k-1 :2*k,:) ;
                DDor= FromTo(2,2)-FromTo(1,2) ;  DDor=DDor/norm(DDor) ;
                Arr= (FromTo(1,2):DDor:FromTo(2,2))' ; Arr=[ FromTo(1,1)*ones(length(Arr),1)  ,Arr ];
                Arr2 = [circshift(Arr,1) , circshift(Arr,-1) ] ;
                Arr2(1,1:2)=-1 ; Arr2(end, 3:4)=-1 ;
                if  DDor==1
                    obj.ScafdigitSQ{FromTo(1,1)}(FromTo(1,2):FromTo(2,2),: ) =Arr2 ;
                    obj.ScafdigitHC{FromTo(1,1)}(FromTo(1,2):FromTo(2,2),: ) =Arr2 ;
                else
                    obj.ScafdigitSQ{FromTo(1,1)}(FromTo(2,2):1:FromTo(1,2),: ) =flip(Arr2 );
                    obj.ScafdigitHC{FromTo(1,1)}(FromTo(2,2):1:FromTo(1,2),: ) =flip(Arr2 );
                end
            end
            %----------------------||||| discuss case
            ForEasy =TargetRCylandDecidePrims ;
            ConnScaf2= cell(size(Info ,1),1) ;
            for k= 1:size(Info ,1)
                
                if strcmp(Info{k,8},'Connected')  &&  ~strcmp(Info{k,2},Info{k,5})
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , ForEasy(2*k-1:2*k ,3 )];
                elseif  strcmp(Info{k,8},'Free end')  &&  ~strcmp(Info{k,2},Info{k,5})
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , ForEasy(2*k-1:2*k ,3 )+ForEasy(2*k-1:2*k ,4 )];
                elseif strcmp(Info{k,8},'Connected')  &&  strcmp(Info{k,2},Info{k,5})
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , [ForEasy(2*k-1 ,3 );  ForEasy(2*k ,3 )+ForEasy(2*k ,4 )    ]   ];
                elseif strcmp(Info{k,8},'Free end')  &&  strcmp(Info{k,2},Info{k,5})
                    MM =  [ForEasy(2*k-1:2*k ,1 ) , [ForEasy(2*k-1 ,3 )+ForEasy(2*k-1 ,4 );  ForEasy(2*k ,3 )    ]   ];
                end
                MM=flip(MM) ;   % first row 5' end , second row 3' end ,otherwise cadnano not show
                %                 k
                FirstRow= obj.ScafdigitSQ{MM(1,1)}(MM(1,2),:);
                Inds=find(FirstRow==-1) ;
                FirstRow(Inds(1)) = MM(2,1) ;   FirstRow(Inds(2)) = MM(2,2) ;
                obj.ScafdigitSQ{MM(1,1)}(MM(1,2),:)=FirstRow ;
                obj.ScafdigitHC{MM(1,1)}(MM(1,2),:)=FirstRow ;
                
                SecondRow= obj.ScafdigitSQ{MM(2,1)}(MM(2,2),:);
                Inds2=find(SecondRow==-1) ;
                SecondRow(Inds2(1)) = MM(1,1) ;   SecondRow(Inds2(2)) = MM(1,2) ;
                obj.ScafdigitSQ{MM(2,1)}(MM(2,2),:)= SecondRow ;
                obj.ScafdigitHC{MM(2,1)}(MM(2,2),:)= SecondRow ;
                
                ConnScaf2{k} =MM ;
                
            end
            obj.ClosingStrand.ExtendScaf=ExtendScaf;
            obj.ClosingStrand.ConnScaf2=ConnScaf2;
        end  % end of  LoadJsonRecoverClosingStrand
        
        
        function Lengths=CalculateSegLengths(obj,Ex,Gcut)
            Lengths=zeros(length(find(Gcut>=1))+1,1);
            %                           Lengths=zeros(length(find(Gcut==1))+1,1);
            
            CC=0;k=1;
            for i=1:length(Ex)
                if Gcut(i)<1
                    CC=CC+Ex(i);
                elseif Gcut(i)>=1
                    A=round(Ex(i)/2);
                    B=Ex(i)-A;
                    CC=CC+A;
                    Lengths(k)=CC;
                    k=k+1;
                    CC=B;
                end
            end
            Lengths(k)=CC;
        end
        function  [RegularInterval,LongWrongInterval]=SPbreakfindleng2Interval(obj,Lengths,Interval,Ex,relax,LongInterval,LongGInfo,Bound)
            %              Bound=[30,60];
            OldInt=zeros(Interval.N+length(LongInterval),2);
            OldInt(1:Interval.N,:)=[Interval.L' Interval.SP'];
            OldInt(Interval.N+1:end,:)=[LongGInfo(LongInterval,2:3)];
            [NInt,index]=sortrows(OldInt,2);
            NIntercal.L=NInt(:,1); NIntercal.SP=NInt(:,2); NIntercal.N=length(NInt(:,1));
            
            QQ=cumsum(Lengths); %in terms of s, size=> total number of edge
            WW=cumsum(Ex);  %in terms of s,
            tableA=zeros(length(Lengths)-1,2);
            %             Lengths
            %             QQ
            Mix=sort([QQ(1),WW]);
            kindex=find(Mix==QQ(1));
            kindex=kindex(1) ;  %
            
            if kindex==1
                kindex=2;
            end
            tA1=find(WW==Mix(kindex-1));
            [size(Mix), kindex] ;
            %             if kindex==length(Mix)
            %             sdf=2;
            %             end
            tA2= find(WW==Mix(kindex+1));
            tableA(1,:)=[tA1, tA2];
            
            for iShortStp=2:length(Lengths)-1
                Mix=sort([QQ(iShortStp),WW]);kindex=find(Mix==QQ(iShortStp));
                tableA(iShortStp,:)=[find(WW==Mix(kindex-1)),  find(WW==Mix(kindex+1))] ;
            end
            TAA(:,1)=tableA(:,2);
            for kthBK=1:size(TAA,1)
                XX=TAA(kthBK,1);
                NextSP=XX<NIntercal.SP;          %subs Intervals with NInterval
                YY=find(NextSP==1);
                if isempty(YY)
                    YY=NIntercal.N+1;
                end
                TAA(kthBK,2)=YY(1)-1;
            end
            %TAA :   first columns cutting points belongs to which edge,
            %2th colums refers to which Intervals
            failLs=or( Lengths>Bound(2), Lengths<Bound(1));
            wrongL=find(failLs);
            if relax>0.7
                IncorrectKnives=union(union(union(-1+wrongL,wrongL),-2+wrongL),1+wrongL)   ;  %in terms of which intervals
            else
                IncorrectKnives=union(-1+wrongL,wrongL) ;  %in terms of which intervals
            end
            IncorrectKnives=setdiff(IncorrectKnives,[-1 0 size(failLs,1) size(failLs,1)+1]);
            
            IncorrectInterval=TAA(IncorrectKnives,2);
            %             IncorrectInterval=union(IncorrectInterval,[]);
            IncorrectInterval=setdiff(IncorrectInterval,0);
            %             IncorrectInterval
            
            ReverseAndDivide=index(IncorrectInterval);
            RegularInterval= ReverseAndDivide(ReverseAndDivide<=Interval.N);
            LongWrongInterval= ReverseAndDivide(ReverseAndDivide>Interval.N)-Interval.N;
            LongWrongInterval=LongInterval(LongWrongInterval);   %in global edge unit
        end   %end of SPbreakfindleng2Interval
        
        function plotDomain(obj,JsonOBJ)
              % use UseCadnano to get JsonOBJ 
             StackScafAllBase = [-1 -1] ;
             for k =1:length(obj.ScafAllBase)
              StackScafAllBase=[StackScafAllBase ;  obj.ScafAllBase{k}];  
             end
            StackScafAllBase=StackScafAllBase(2:end,:)  ;
            
            figure(2301); clf ; hold on ;
            
            
            for k = 1: length( JsonOBJ.DomainStapOnScaf)
              for j = 1: length(JsonOBJ.DomainStapOnScaf{k})
                IndOnScaf = JsonOBJ.DomainStapOnScaf{k}{j} ;
                ScafIndC5=StackScafAllBase(IndOnScaf(:,1) ,:) ;
                
                BaseScafR3DXYZ =obj.plotScafR_BaseByBase(ScafIndC5) ;
                plot3(BaseScafR3DXYZ(:,1),BaseScafR3DXYZ(:,2),BaseScafR3DXYZ(:,3)) ;
              end
            end
            
        end
        
    end
    
    
end

