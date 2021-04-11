%-----script for connecting staples in the ends and having half-crossover
%input: StapleCell    output:StapleCell

if scriptcase==1
    %     Nicks=[obj.StapNick(:,1:2) ;obj.StapNick(:,3:4) ];
    %              if ~isempty(obj.StapNick)
    %               Nicks=[obj.StapNick(:,1:2) ;obj.StapNick(:,3:4) ];
    %              else
    Nicks=[-1,-1];
    %              end
    
    StapleCell  ;
    OStapleCell=StapleCell;
    
    SaveGHelixStap=cell(1,length(obj.containBundle));
    for Bundlei=1:length(obj.containBundle)
        QQWithPM10Bases =obj.containBundle{Bundlei}.HelixXYZGStap;
        SaveGHelixStap{Bundlei}= QQWithPM10Bases;
    end
    HeadToTail=[-1,-1,0] ;epson=1.4;   %1.8   index for OStapleCell,  not StapleCell,
    HalfXover=[-1,-1,-1];
    for stpi=1:length(StapleCell)
        for stpj=stpi+1:length(StapleCell)
            StpI=StapleCell{stpi};
            StpJ=StapleCell{stpj};
            
            %             stpi
            %
            %             if stpj==11
            %             stpsdfj=0;
            %             end
            StpIBCB=zeros(size(StpI,1),3);  StpIBCB(:,3)=StpI(:,2);
            for k=1:size(StpI,1)
                Inds=obj.RelateTable(:,5)==StpI(k,1) ;
                if sum(Inds)>1
                    Bundle =unique(obj.RelateTable(Inds,1) );
                    Clyinders =obj.RelateTable(Inds,2) ;
                    for cylk=1: length(Clyinders)
                        Val =(obj.containBundle{Bundle}.Zbase1(Clyinders(cylk))- StpI(k,2))*(obj.containBundle{Bundle}.Zbase2(Clyinders(cylk))- StpI(k,2)) ;
                        if  Val<=0
                            StpIBCB(k,1:2)=[Bundle ,Clyinders(cylk)];   % for changing crosssection
                        end
                    end
                    
                else
                    StpIBCB(k,1:2)= obj.RelateTable(Inds,1:2);
                end
                %             StpIBCB(k,1:2)= obj.RelateTable(obj.RelateTable(:,5)==StpI(k,1),1:2);
            end
            
            StpJBCB=zeros(size(StpJ,1),3); StpJBCB(:,3)=StpJ(:,2);
            for k=1:size(StpJ,1)
                Inds=obj.RelateTable(:,5)==StpJ(k,1) ;
                if sum(Inds)>1
                    Bundle =unique(obj.RelateTable(Inds,1) );
                    Clyinders =obj.RelateTable(Inds,2) ;
                    for cylk=1: length(Clyinders)
                        Val =(obj.containBundle{Bundle}.Zbase1(Clyinders(cylk))- StpJ(k,2))*(obj.containBundle{Bundle}.Zbase2(Clyinders(cylk))- StpJ(k,2)) ;
                        if  Val<=0
                            StpJBCB(k,1:2)=[Bundle ,Clyinders(cylk)];   % for changing crosssection
                        end
                    end
                else
                    StpJBCB(k,1:2)= obj.RelateTable(Inds,1:2);
                end
            end
            
            %----check I-head with J-tail
            RelativeIH=StpIBCB(1,3)-obj.containBundle{StpIBCB(1,1)}.Zbase1(StpIBCB(1,2))+11 ;
            P_Ihead=SaveGHelixStap{StpIBCB(1,1)}{StpIBCB(1,2)}(RelativeIH,:) ;
            
            RelativeJT=StpJBCB(end,3)-obj.containBundle{StpJBCB(end,1)}.Zbase1(StpJBCB(end,2))+11 ;
            P_JTail=SaveGHelixStap{StpJBCB(end,1)}{StpJBCB(end,2)}(RelativeJT,:) ;
            if norm(P_Ihead-P_JTail) < epson   &&  StpJBCB(end,1)== StpIBCB(1,1)  &&  StpJBCB(end,2)~= StpIBCB(1,2) && abs( StpJBCB(end,3)- StpIBCB(1,3))<2     % probably need change
                if  ~ismember(StpI(1,:),Nicks,'rows') && ~ismember(StpJ(end,:),Nicks,'rows') && ~ismember(stpi,HeadToTail(:,1)) && ~ismember(stpj,HeadToTail(:,2))  %&& StpJBCB(1,1)~=6  % hard
                    StpI_CylinderEnds= [ obj.containBundle{StpIBCB(1,1)}.Zbase1(StpIBCB(1,2)) ,  obj.containBundle{StpIBCB(1,1)}.Zbase2(StpIBCB(1,2)) ] ;
                    if sum(abs( StpIBCB(1,3) -StpI_CylinderEnds)<5 )>=1 % 05222019, prevent internal halXover
                        HeadToTail=union(HeadToTail,[stpi,stpj,0 ],'rows') ;
                        HalfXover=union(HalfXover,[StpJ(end,1), StpI(1,:)],'rows' ) ;
                    end
                end
            end
            %%%  add && abs( StpJBCB(end,3)- StpIBCB(1,3))<2 this criterion on
            %%%  Oct 30 2018, filter out half xover where bases are within
            %%%  certain range.
            %----check I-tail with J-head
            RelativeIT=StpIBCB(end,3)-obj.containBundle{StpIBCB(end,1)}.Zbase1(StpIBCB(end,2))+11 ;
            P_ITail=SaveGHelixStap{StpIBCB(end,1)}{StpIBCB(end,2)}(RelativeIT,:) ;
            
            RelativeJH=StpJBCB(1,3)-obj.containBundle{StpJBCB(1,1)}.Zbase1(StpJBCB(1,2))+11 ;
            P_JHead=SaveGHelixStap{StpJBCB(1,1)}{StpJBCB(1,2)}(RelativeJH,:) ;
            if norm(P_ITail-P_JHead) < epson  &&  StpJBCB(1,1)== StpIBCB(end,1)  &&  StpJBCB(1,2)~= StpIBCB(end,2) && abs( StpJBCB(1,3)- StpIBCB(end,3))<2    % probably need change
                if  ~ismember(StpI(end,:),Nicks,'rows') && ~ismember(StpJ(1,:),Nicks,'rows') && ~ismember(stpi,HeadToTail(:,2)) && ~ismember(stpj,HeadToTail(:,1))  % && StpJBCB(1,1)~=6
                    StpJ_CylinderEnds= [ obj.containBundle{StpJBCB(1,1)}.Zbase1(StpJBCB(1,2)) ,  obj.containBundle{StpJBCB(1,1)}.Zbase2(StpJBCB(1,2)) ] ;
                    if sum(abs( StpJBCB(1,3) -StpJ_CylinderEnds)<5 )>=1   % 05222019, prevent internal halXover
                        HeadToTail=union(HeadToTail,[stpj,stpi,1 ],'rows') ;
                        HalfXover=union(HalfXover,[StpI(end,1), StpJ(1,:)],'rows'   ) ;
                    end
                end
            end
        end
    end
    HeadToTail=setdiff(HeadToTail,[-1,-1,0], 'rows');
    HalfXover=setdiff(HalfXover,[-1,-1,-1], 'rows');
    
    %     Gr = digraph(HeadToTail(:,1),HeadToTail(:,2));
    %--------June 14
    CellGroup=cell(length(StapleCell),1);
    for puti=1:length(StapleCell); CellGroup{puti}=puti;end
    
    for k=1 :size(HeadToTail,1)
        Element1=HeadToTail(k,1);
        Element2=HeadToTail(k,2);
        CellIndex=zeros(1,2);
        for sCell=1:length(CellGroup)
            if ismember(Element1, CellGroup{sCell}); CellIndex(1)=sCell ;end
            if ismember(Element2, CellGroup{sCell}); CellIndex(2)=sCell ;end
        end
        
        if CellIndex(1)~=CellIndex(2)
            
            CellGroup{CellIndex(1)}= [CellGroup{CellIndex(1)} ,CellGroup{CellIndex(2)}];
            CellGroup{CellIndex(2)}=[];
        end
    end
    CellGroup=CellGroup(~cellfun('isempty',CellGroup)) ;
    
    cc=[];
    for qq=1:length(CellGroup)
        cc=union(cc,CellGroup{qq});
    end
    
    for celli=1:length(CellGroup)
        Vec=CellGroup{celli};
        for connj=1:length(Vec)-1
            if StapleCell{Vec(connj) }(end,2)~= StapleCell{ Vec(connj+1)}(1,2)
                sdfs=34;
            end
            
            
            StapleCell{ Vec(connj+1)}=[ StapleCell{Vec(connj+1) } ;   StapleCell{ Vec(connj)}];
            StapleCell{ Vec(connj) };
            StapleCell{ Vec(connj) }=0;
        end
        
    end
    shouldberemove=[];
    for iC=1:length(StapleCell)
        if length(StapleCell{iC}(:))==1
            shouldberemove=union(shouldberemove,iC);
        end
    end
    
    StapleCell=StapleCell(setdiff(1:length(StapleCell),shouldberemove)  );
    %     %------
    %
    %
    %     for Conni=1:size(HeadToTail,1)
    %       StapleCell{ HeadToTail(Conni,2)}=[ OStapleCell{ HeadToTail(Conni,2)} ;  OStapleCell{ HeadToTail(Conni,1)}];
    % %       if isempty(StapleCell{ HeadToTail(Conni,2)} )
    % %           sdsf=3
    % %       elseif isempty(StapleCell{ HeadToTail(Conni,1)} )
    % %           sdfsf=34
    % %       else
    % %       StapleCell{ HeadToTail(Conni,2)}=[ StapleCell{ HeadToTail(Conni,2)} ;  StapleCell{ HeadToTail(Conni,1)}];
    % %       end
    %
    %       StapleCell{ HeadToTail(Conni,1)}=[];
    %
    %     end
    %
    
    if ~exist('StpJ')
        StpJ=[];
        %     StpJ
    else
        StpJBCB=zeros(size(StpJ,1),3);
        StpJBCB(:,3)=StpJ(:,2);
    end
    
    %     StapleCell=StapleCell(~cellfun('isempty',StapleCell)) ;
    %-----staple self-looping
    stapSelfLoop= -1;
    for stpk=1:  length(StapleCell)
        stpK=StapleCell{stpk};
        StpKBCB=zeros(size(stpK,1),3);  StpKBCB(:,3)=stpK(:,2);
        for k=1:size(stpK,1)
            Inds=obj.RelateTable(:,5)==stpK(k,1) ;
            if sum(Inds)>1
                Bundle =unique(obj.RelateTable(Inds,1) );
                Clyinders =obj.RelateTable(Inds,2) ;
                for cylk=1: length(Clyinders)
                    Val =(obj.containBundle{Bundle}.Zbase1(Clyinders(cylk))- stpK(k,2))*(obj.containBundle{Bundle}.Zbase2(Clyinders(cylk))- stpK(k,2)) ;
                    if  Val<=0
                        StpKBCB(k,1:2)=[Bundle ,Clyinders(cylk)];   % for changing crosssection
                    end
                end
                
            else
                StpKBCB(k,1:2)= obj.RelateTable(Inds,1:2);
            end
            
            
            %             StpKBCB(k,1:2)= obj.RelateTable(obj.RelateTable(:,5)==stpK(k,1),1:2);
        end
        
        
        for k=1:size(StpJ,1)
            Inds=obj.RelateTable(:,5)==StpJ(k,1) ;
            if sum(Inds)>1
                Bundle =unique(obj.RelateTable(Inds,1) );
                Clyinders =obj.RelateTable(Inds,2) ;
                for cylk=1: length(Clyinders)
                    Val =(obj.containBundle{Bundle}.Zbase1(Clyinders(cylk))- StpJ(k,2))*(obj.containBundle{Bundle}.Zbase2(Clyinders(cylk))- StpJ(k,2)) ;
                    if  Val<=0
                        StpJBCB(k,1:2)=[Bundle ,Clyinders(cylk)];   % for changing crosssection
                    end
                end
            else
                StpJBCB(k,1:2)= obj.RelateTable(Inds,1:2);
            end
            
            
            %             StpJBCB(k,1:2)= obj.RelateTable(obj.RelateTable(:,5)==StpJ(k,1),1:2);
        end
        
        %----check I-head with J-tail
        RelativeIH=StpKBCB(1,3)-obj.containBundle{StpKBCB(1,1)}.Zbase1(StpKBCB(1,2))+11 ;
        P_Ihead=SaveGHelixStap{StpKBCB(1,1)}{StpKBCB(1,2)}(RelativeIH,:) ;
        %----check I-tail with J-head
        RelativeIT=StpKBCB(end,3)-obj.containBundle{StpKBCB(end,1)}.Zbase1(StpKBCB(end,2))+11 ;
        P_ITail=SaveGHelixStap{StpKBCB(end,1)}{StpKBCB(end,2)}(RelativeIT,:) ;
        
        HT=StpKBCB([1,end],:) ;
        if norm(P_ITail-P_Ihead) < epson && size( unique(HT(:,1:2),'rows'),1 )==2
            stapSelfLoop=union(stapSelfLoop,stpk ) ;
        end
    end
    
    
    stapSelfLoop=setdiff(stapSelfLoop, -1)  ;
    for k=1:length(stapSelfLoop)
        stpK=StapleCell{stapSelfLoop(k)};
        Uvec= stpK(1,:)-stpK(2,:) ; LU=norm(Uvec)  ;Uvec=Uvec/LU ;
        NewS= [ stpK(1,:) ;  stpK(1,:)-ceil(LU/2)*Uvec ;  stpK(1,:)-(1+ceil(LU/2))*Uvec  ;stpK(2:end,:) ] ;
        
        if abs(diff(stpK(1:2,2)))<4      % b.c. new clr can be 2nt, avoid bugs from face-face Xovers       
          StapleCell{stapSelfLoop(k)} =    stpK ;
        else
         StapleCell{stapSelfLoop(k)}= circshift(NewS,-2) ;
           
        end
        
    end
    
%     for ddbug = 1:length(StapleCell)
%         if StapleCell{ddbug}(1,1)==StapleCell{ddbug}(2,1) && StapleCell{ddbug}(1,2)==StapleCell{ddbug}(2,2)
%         ddbug
%         end
%     end
        
    AddHalfXoverPart1=1;
    
elseif scriptcase==2 %----------------------------------------------
    OOstapBP=obj.stapBP;
    
    %-------Z1 sdie
    AdjM=obj.CylinderC5AdjM;    %use CylAdjM in hyperbundle
    [U,V]=find(AdjM~=0);   %means cylinder in C5
    U2=U;V2=V;  %means cylinder in C5
    Edgeof5Index=  unique([U2 V2],'rows');  %C5
    k=1;
    EdgeList=cell(size(Edgeof5Index,1)/2,2);
    
    for i=1:size(Edgeof5Index,1)
        Ind=find(Edgeof5Index(i,1)==obj.RelateTable(:,5));Ind=Ind(1);
        Bundle=obj.containBundle{obj.RelateTable(Ind,1)};
        if ~xor(ismember( obj.RelateTable(Ind,2)   ,Bundle.AGroup) ,Bundle.AGroupGoUp==1)
            EdgeList{k,1}=[ Edgeof5Index(i,1) , Edgeof5Index(i,2)];   %in 5rd index;
            k=k+1;
        end
    end
    ExtraStapBP=[];
    CCZ=[];
    for edgei=1:size(EdgeList,1)
        CylMaster=  EdgeList{edgei,1}(2) ;         CylInSlave= EdgeList{edgei,1}(1);  %C5 Express
        Bundle= obj.containBundle{obj.RelateTable(obj.RelateTable(:,5)==CylMaster,1)} ;
        CylM= obj.RelateTable(obj.RelateTable(:,5)==CylMaster,2) ;  % bundle -"Cylinder"
        CylS= obj.RelateTable(obj.RelateTable(:,5)==CylInSlave,2) ;
        
        CinMoveUp=  ~xor( ismember(CylM,Bundle.AGroup),Bundle.AGroupGoUp);
        ZLAll=[];
        mmZ= min(Bundle.Zbase1([CylM,CylS])) ;
        for Zscan=mmZ: mmZ +20
            %       zleft= FindZ2( CylM,CylS,Bundle.CylInplanePosition, Zscan ,~CinMoveUp,[],[] ) ;
            if strcmp(Bundle.Lattice, 'Square')            %change for hybrid structure
                zleft= FindZ2SQ( CylM,CylS,Bundle.CylInplanePosition,Zscan,~CinMoveUp);
            else
                zleft= FindZ2HC( CylM,CylS,Bundle.CylInplanePosition,Zscan,~CinMoveUp);
            end
            
            
            ZLAll=union(ZLAll,zleft) ;
        end
        ZLAll(ZLAll-mmZ>15)=[]  ;    %----------++++++++++++++++++++++++
        ZLAll(ZLAll-mmZ<4)=[]  ;     %----------++++++++++++++++++++++++
        for i=1:length( ZLAll)
            ExtraStapBP=[ExtraStapBP; CylInSlave,ZLAll(i); CylMaster,ZLAll(i) ; CylInSlave,ZLAll(i)+1; CylMaster,ZLAll(i)+1 ];
            CCZ=[CCZ ; CylMaster,CylInSlave,ZLAll(i)];
        end
        EdgeList{edgei,2}= ZLAll;
    end
    %-----------Z2 side
    CCZ2=[];
    for edgei=1:size(EdgeList,1)
        CylMaster=  EdgeList{edgei,1}(2) ;         CylInSlave= EdgeList{edgei,1}(1);  %C5 Express
        Bundle= obj.containBundle{obj.RelateTable(obj.RelateTable(:,5)==CylMaster,1)} ;
        CylM= obj.RelateTable(obj.RelateTable(:,5)==CylMaster,2) ;  % bundle -"Cylinder"
        CylS= obj.RelateTable(obj.RelateTable(:,5)==CylInSlave,2) ;
        
        CinMoveUp=  ~xor( ismember(CylM,Bundle.AGroup),Bundle.AGroupGoUp);
        ZLAll2=[];
        MMZ= max(Bundle.Zbase2([CylM,CylS])) ;
        for Zscan=MMZ-20: MMZ
            %       zleft= FindZ2( CylM,CylS,Bundle.CylInplanePosition, Zscan ,~CinMoveUp,[],[] ) ;
            if strcmp(Bundle.Lattice, 'Square')            %change for hybrid structure
                zleft= FindZ2SQ( CylM,CylS,Bundle.CylInplanePosition,Zscan,~CinMoveUp);
            else
                zleft= FindZ2HC( CylM,CylS,Bundle.CylInplanePosition,Zscan,~CinMoveUp);
            end
            
            
            ZLAll2=union(ZLAll2,zleft) ;
        end
        ZLAll2(-ZLAll2+MMZ>15)=[]  ;    %----------++++++++++++++++++++++++
        ZLAll2(-ZLAll2+MMZ<4)=[]  ;     %----------++++++++++++++++++++++++
        for i=1:length( ZLAll2)
            ExtraStapBP=[ExtraStapBP; CylInSlave,ZLAll2(i); CylMaster,ZLAll2(i) ; CylInSlave,ZLAll2(i)+1; CylMaster,ZLAll2(i)+1 ];
            CCZ2=[CCZ2 ; CylMaster,CylInSlave,ZLAll2(i)];
        end
        EdgeList{edgei,2}= ZLAll2;
    end
    %--------------------
    
    
    
    
    
    HalfXover ;
    %     HeadToTail;
    MoreStapBP=[];
    for Bi=1 : length(obj.containBundle)
        %-----------------Z1 side
        AllCylInC5= obj.RelateTable(obj.RelateTable(:,1)== Bi,5) ;
        VecToStartOne=AllCylInC5;VecToStartOne(:,2)=VecToStartOne(:,1);
        VecToStartOne(:,2)=VecToStartOne(:,2)-min(VecToStartOne(:,2))+1 ;
        HalfXoverInThis=HalfXover(ismember(HalfXover(:,1),AllCylInC5),:) ;
        
        HalfXoverInThisZ1=       HalfXoverInThis ;
        HalfXoverInThisZ1(HalfXoverInThis(:,3)>100,:)=[] ;%------------------
        
        InplaneOri=unique(HalfXoverInThisZ1(:,1:2),'rows' ) ;
        IndFlip=InplaneOri(:,1)>InplaneOri(:,2);InplaneOri(IndFlip,:)= flip(InplaneOri(IndFlip,:),2) ;
        InplaneOri=unique(InplaneOri(:,1:2),'rows') ;
        
        %        InPool=intersect(unique(HalfXover(:,1:2)), AllCylInC5);
        %        InPool=unique(InplaneOri);  %----------
        
        PXoverList= or( ismember(CCZ(:,1),AllCylInC5),ismember(CCZ(:,2),AllCylInC5)) ;   %CM
        XoverList=CCZ(PXoverList,:) ;       IndFlip3=XoverList(:,1)>XoverList(:,2);
        XoverList(IndFlip3,1:2)= flip(XoverList(IndFlip3,1:2),2) ;
        
        InplaneList=unique(XoverList(:,1:2),'rows') ;
        IndFlip2=InplaneList(:,1)>InplaneList(:,2);InplaneList(IndFlip2,:)= flip(InplaneList(IndFlip2,:),2) ;
        InplaneList=setdiff(InplaneList,InplaneOri,'rows');
        
        nEdgeRequired= length(AllCylInC5) -  size(InplaneOri,1)-1 ;
        nwh=0; NCase=1:size(InplaneList,1);
        while 1
            RInd= randsample(NCase,nEdgeRequired) ; %index for InplaneList
            CylindersInvolved= union(unique(InplaneList(RInd,:)) , unique(InplaneOri));
            if length(CylindersInvolved)==length(AllCylInC5) % all cylinder involed
                AdjList=[InplaneOri;InplaneList(RInd,:)]; AdjListStartOne=zeros(size(AdjList));
                for k=1:size(VecToStartOne,1)
                    IndSS=AdjList==VecToStartOne(k,1);
                    AdjListStartOne(IndSS)=VecToStartOne(k,2);
                end
                Gr=graph(AdjListStartOne(:,1),AdjListStartOne(:,2));
                bins = conncomp(Gr);
                if sum(bins==1)==length(bins)  %no separate group
                    break
                end
            end
            nwh=nwh+1;
            if nwh>=50
                sdfsf=0;
            end
            
        end
        AdjList=InplaneList(RInd,:);
        
        checkM=[AdjList;InplaneOri];
        IInd=checkM(:,2)> checkM(:,1) ;checkM(IInd,:)= flip(checkM(IInd,:),2) ;
        checkM;
        if  size(unique(checkM, 'rows'),1)~= length(AllCylInC5)-1
            sdfsf=34;
        end
        
        sdfsff=34;
        %---------AdjList
        NeedExtra=setdiff(AdjList,InplaneOri,'rows');
        [~,Locb] =ismember(XoverList(:,1:2),NeedExtra,'rows') ; nwhile=1;
        while 1
            ExtraXover=zeros(max(Locb),3);
            for ck=1:max(Locb)
                Ind= find(Locb==ck) ;
                Ind=Ind(randi(length(Ind)));
                ExtraXover(ck,:)=XoverList(Ind,:);
            end
            All=[HalfXoverInThis; ExtraXover ] ;
            Good=1;
            for Cylk=1:length(AllCylInC5)
                InvolvedInd=or( ismember(All(:,1),Cylk),ismember(All(:,2),Cylk)) ;
                ZSep=All(InvolvedInd,3);
                if length(ZSep)>2
                    ZSep=sort(ZSep);
                    if min(ZSep(2:end)-ZSep(1:end-1))<5;  Good=0;   end
                end
            end
            if  Good==1 || nwhile>=50 ;  break;   end
            nwhile=nwhile+1;
        end
        QQ=zeros(size(ExtraXover,1)*4,2);
        for qq=1:size(ExtraXover,1)
            QQ(4*qq-3:4*qq,:)=[ExtraXover(qq,1),ExtraXover(qq,3);ExtraXover(qq,2),ExtraXover(qq,3);ExtraXover(qq,1),1+ExtraXover(qq,3);ExtraXover(qq,2),1+ExtraXover(qq,3)];
        end
        stapBPExtra=[stapBPExtra ; QQ];
        
        %        %----------------Z2 side
        % %
        AllCylInC5= obj.RelateTable(obj.RelateTable(:,1)== Bi,5) ;
        VecToStartOne=AllCylInC5;VecToStartOne(:,2)=VecToStartOne(:,1);
        VecToStartOne(:,2)=VecToStartOne(:,2)-min(VecToStartOne(:,2))+1 ;
        HalfXoverInThis=HalfXover(ismember(HalfXover(:,1),AllCylInC5),:) ;
        
        HalfXoverInThisZ2=       HalfXoverInThis ;
        HalfXoverInThisZ2(HalfXoverInThis(:,3)<100,:)=[] ;%------------------
        
        
        %        HalfXoverInThis(HalfXoverInThis(:,3)>200,:)=[] ;
        InplaneOri=unique(HalfXoverInThisZ2(:,1:2),'rows' ) ;
        IndFlip=InplaneOri(:,1)>InplaneOri(:,2);InplaneOri(IndFlip,:)= flip(InplaneOri(IndFlip,:),2) ;
        InplaneOri=unique(InplaneOri(:,1:2),'rows') ;
        
        %        InPool=unique(InplaneOri);  %----------
        PXoverList= or( ismember(CCZ2(:,1),AllCylInC5),ismember(CCZ2(:,2),AllCylInC5)) ;
        XoverList=CCZ2(PXoverList,:) ;
        IndFlip3=XoverList(:,1)>XoverList(:,2);XoverList(IndFlip3,1:2)= flip(XoverList(IndFlip3,1:2),2) ;
        
        InplaneList=unique(XoverList(:,1:2),'rows') ;
        IndFlip2=InplaneList(:,1)>InplaneList(:,2);InplaneList(IndFlip2,:)= flip(InplaneList(IndFlip2,:),2) ;
        InplaneList=setdiff(InplaneList,InplaneOri,'rows');
        
        nEdgeRequired= length(AllCylInC5) -  size(InplaneOri,1)-1 ;
        nwh=0; NCase=1:size(InplaneList,1);
        while 1
            RInd= randsample(NCase,nEdgeRequired) ; %index for InplaneList
            CylindersInvolved= union(unique(InplaneList(RInd,:)) , unique(InplaneOri));
            if length(CylindersInvolved)==length(AllCylInC5) % all cylinder involed
                AdjList=[InplaneOri;InplaneList(RInd,:)]; AdjListStartOne=zeros(size(AdjList));
                for k=1:size(VecToStartOne,1)
                    IndSS=AdjList==VecToStartOne(k,1);
                    AdjListStartOne(IndSS)=VecToStartOne(k,2);
                end
                Gr=graph(AdjListStartOne(:,1),AdjListStartOne(:,2));
                bins2 = conncomp(Gr);
                if sum(bins2==1)==length(bins2)  %no separate group
                    break
                end
            end
            
            nwh=nwh+1;
            if nwh>=50
                sdfsf=0;
            end
            
        end
        AdjList=InplaneList(RInd,:);
        
        
        
        %        RemainCyl=setdiff(AllCylInC5,InPool) ;
        %        AddEdge=RemainCyl; AddEdge(:,end+1)=0;
        %        for k=1:size(RemainCyl)
        %            Cyl=AddEdge(k,1);
        %            PossibleParents=find(obj.CylinderC5AdjM(Cyl,:)==1) ;
        %            PossibleParents=setdiff(PossibleParents,RemainCyl) ;
        %            AddEdge(k,2)=  PossibleParents(randi(length(PossibleParents))) ;
        %        end
        %        AdjList=AddEdge;
        %        IndFlip=AdjList(:,1)>AdjList(:,2);AdjList(IndFlip,:)= flip(AdjList(IndFlip,:),2) ;
        %
        %
        
        %        AllPossible = nchoosek(1:size(InplaneList,1),length(AllCylInC5)-1-size(InplaneOri,1));
        %        for rrp=1:100
        %            RandExtra=AllPossible(randi(size(AllPossible,1)),: ) ;
        %            AdjList=[InplaneOri ; InplaneList(RandExtra,:) ] ;
        %
        %           GraphEdge = graph(AdjList(:,1)',AdjList(:,2)');
        %            bins = conncomp(GraphEdge);
        %            if max(bins)==1
        %                break
        %            end
        %        end
        
        
        NeedExtra=setdiff(AdjList,InplaneOri,'rows');
        [~,Locb] =ismember(XoverList(:,1:2),NeedExtra,'rows') ; nwhile=1;
        while 1
            ExtraXover=zeros(max(Locb),3);
            for ck=1:max(Locb)
                Ind= find(Locb==ck) ;   Ind=Ind(randi(length(Ind)));
                ExtraXover(ck,:)=XoverList(Ind,:);
            end
            All=[HalfXoverInThis; ExtraXover ] ;
            Good=1;
            for Cylk=1:length(AllCylInC5)
                InvolvedInd=or( ismember(All(:,1),Cylk),ismember(All(:,2),Cylk)) ;
                ZSep=All(InvolvedInd,3);
                if length(ZSep)>2
                    ZSep=sort(ZSep);
                    if min(ZSep(2:end)-ZSep(1:end-1))<5;  Good=0;   end
                end
            end
            if  Good==1 || nwhile>=50 ;  break;   end
            nwhile=nwhile+1;
        end
        QQ=zeros(size(ExtraXover,1)*4,2);
        for qq=1:size(ExtraXover,1)
            QQ(4*qq-3:4*qq,:)=[ExtraXover(qq,1),ExtraXover(qq,3);ExtraXover(qq,2),ExtraXover(qq,3);ExtraXover(qq,1),1+ExtraXover(qq,3);ExtraXover(qq,2),1+ExtraXover(qq,3)];
        end
        stapBPExtra=[stapBPExtra ; QQ];
        
        
        %
        
        
        
        
        
        
        
        
    end
    
    
    
end




