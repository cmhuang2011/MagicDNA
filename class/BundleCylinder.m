classdef BundleCylinder < handle
    %description:
    %     put all common function and properties of SQ and HC here.
    %  as one of the components in hyperbundle
    %   Detailed explanation goes here
    
    properties
        Zbase1=[];     %unit base
        Zbase2=[];
        CylInplanePosition=[];
        %--------------------below are new from MagicDNA2
        skipPosition =[]; %from caDNAno ,or use OOP and DNACurve to improve % [Bundle Cylinder C5 base]
        Default_skipPattern1=[] ;   % skip mod2 =0
        Default_skipPattern2=[] ;
        
        InsertPosition =[]; %from caDNAno ,or use OOP and DNACurve to improve % [Bundle Cylinder C5 base]
        Default_Insert=[] ;   % skip mod2 =0
        
        LocalCoordinatFromLineModel =[];
        %----
        
    end
    
    properties (Hidden = true)  %starting from MagicDNA2
        SR_2p5 = 2.5/2.0 ;
        CylAdjMat=[];
        NarrowAdj=[];
        BottomAdj=[];
        TopAdj=[];
        
%         Tol=15;  %pair tolerence
                Tol=40;  %pair tolerence
        maxTol=40;
        CylRadius=1;
        ssScafOver=5;
        
        AGroup=[];
        BGroup=[];
        CylCategory=[];
        CylABbyCategory=[];
        %         AGroupGoUp=1;   %default value
        %----
        ExternalXoverAsFB=[];         %external xover position as forced-connection of each bundle
        TransformMatrix2=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
        SimulateTransMFromTM2=[];
    end
    
    properties (Dependent, Hidden)
        Z1;     %in nm
        Z2;
        CylinderXYZGlobal;
        TransformMatrix;  %4 by 4
        HelixXYZG;
        HelixXYZGBVec;
        HelixXYZGStap;
        ScaledCylInplanePosition ;
    end
    
    
    methods
        function ScaledCylInplanePosition = get.ScaledCylInplanePosition(obj)
            meanCenter = mean(obj.CylInplanePosition) ;
            ScaledCylInplanePosition = (obj.CylInplanePosition - meanCenter)*obj.SR_2p5 + meanCenter;
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
        function HCoor=findHelixQ(obj,QCylinder,QBasesinCell)
            if length(QCylinder)~=length(QBasesinCell)
                HCoor=[];
                return
            end
            rr=0.85;
            nQ=length(QBasesinCell);      %
            HCoor=cell(nQ,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=obj.CorrectComst;   % check with Cadnano 12/6/2016
            SInplain = obj.ScaledCylInplanePosition ;
            for i=1:nQ
                Cly=QCylinder(i);
                ZZrange=QBasesinCell(i);
                if ~xor(ismember(Cly,obj.AGroup) ,obj.AGroupGoUp)
                    CC=mod(ZZrange-correct(1),obj.period(1))*obj.period(2)*2*pi/obj.period(1);  %theta
                    x0=SInplain(Cly,1) ; y0=SInplain(Cly,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                else
                    CC=mod(ZZrange-correct(2),obj.period(1))*obj.period(2)*2*pi/obj.period(1)  ; %theta
                    x0=SInplain(Cly,1) ; y0=SInplain(Cly,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                end
                HCoor{i}=GlobalXYZ(:,1:3);
            end
        end  % end of findHelixQ
        function HCoor=get.HelixXYZGStap(obj)
            n=length(obj.Zbase1);      % orientation haven't check
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=flip(obj.CorrectComst);
            bpnmConst=0.34;
            rr=0.85;
            %             rr=1;
            EnCylSpace=1;   % 2.3/2
            SInplain = obj.ScaledCylInplanePosition ;
            for i=1:n
                if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                    %                  ZZrange=(obj.Zbase1(i):1:obj.Zbase2(i))';
                    ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                    CC=mod(ZZrange-correct(1),obj.period(1))*obj.period(2)*2*pi/obj.period(1);  %theta
                    %                  x0=EnCylSpace*obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=EnCylSpace*obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1) ; y0=EnCylSpace*SInplain(i,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*bpnmConst ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                    %                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;
                else
                    %                  ZZrange=(obj.Zbase1(i):1:obj.Zbase2(i))';
                    ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                    CC=mod(ZZrange-correct(2),obj.period(1))*obj.period(2)*2*pi/obj.period(1)  ; %theta
                    %                  x0=EnCylSpace*obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=EnCylSpace*obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1) ; y0=EnCylSpace*SInplain(i,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*bpnmConst ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                    %                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;
                end
                HCoor{i}=GlobalXYZ(:,1:3);
            end
        end  % end of HelixXYZGStap
        function HCoor=HelixXYZGStapNoPM(obj,rr,Isstap)
            n=length(obj.Zbase1);      % orientation haven't check
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            if Isstap==1     %if is staple domain
                correct=flip(obj.CorrectComst);
            else       %if is scaffold domain
                correct=obj.CorrectComst;
            end
            bpnmConst=0.34;
            EnCylSpace=1;   % 2.3/2
            SInplain = obj.ScaledCylInplanePosition ;
            for i=1:n
                if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                    ZZrange=(obj.Zbase1(i):1:obj.Zbase2(i))';
                    CC=mod(ZZrange-correct(1),obj.period(1))*obj.period(2)*2*pi/obj.period(1);  %theta
                    %                  x0=EnCylSpace*obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=EnCylSpace*obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1) ; y0=EnCylSpace*SInplain(i,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*bpnmConst ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                else
                    ZZrange=(obj.Zbase1(i):1:obj.Zbase2(i))';
                    CC=mod(ZZrange-correct(2),obj.period(1))*obj.period(2)*2*pi/obj.period(1)  ; %theta
                    %                  x0=EnCylSpace*obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=EnCylSpace*obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1) ; y0=EnCylSpace*SInplain(i,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*bpnmConst ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                end
                HCoor{i}=GlobalXYZ(:,1:3);
            end
        end  % end of HelixXYZGStap
        function XYZ=HelixRegardlessCylinder(obj,rr,Isstap,InplaneXY,queryCylBase,ThefSimilarDirCylder)
            %   queryCylBase:, Nx1 , [ BaseArr] , a series of bases on one
            %   single cylinder(specify by InplaneXY in RTable), i.e. 'a section'
            %             n=size(queryCylBase);

            if Isstap==1     %if is staple domain
                correct=flip(obj.CorrectComst);
            else       %if is scaffold domain  , should not happen so far
                correct=obj.CorrectComst;
            end
            %             rr=0.6;  %-----------debug
            PP = obj.period  ;
%             PP =[72 7 ] ;  % hard
            InplaneXYOri =InplaneXY ;
            meanCenter = mean(obj.CylInplanePosition) ;
            InplaneXY = (InplaneXY - meanCenter)*obj.SR_2p5 + meanCenter;
            
            bpnmConst=0.34;     EnCylSpace=1;   % 2.3/2
            x0=EnCylSpace*InplaneXY(:,1); y0=EnCylSpace*InplaneXY(:,2);
            ZZrange=reshape(queryCylBase, length(queryCylBase) ,1);
            if ~xor(ismember(ThefSimilarDirCylder,obj.AGroup) ,obj.AGroupGoUp)
                CC=mod(ZZrange-correct(1),PP(1))*PP(2)*2*pi/PP(1);  %theta
            else
                CC=mod(ZZrange-correct(2),PP(1))*PP(2)*2*pi/PP(1);  %theta , fixed 08/24
            end
            [x,y] = pol2cart(CC,rr*0.95);
            LocalXYZ=[x+x0, y+y0, ZZrange*bpnmConst ,ones(size(x))];
            XYZ=   transpose(obj.TransformMatrix2*LocalXYZ' );
        end
        
        
        function HCoor=HelixXYZGStapNoPMLocal(obj,rr,Isstap)
            n=length(obj.Zbase1);
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            if Isstap==1     %if is staple domain
                correct=flip(obj.CorrectComst);
            else       %if is scaffold domain
                correct=obj.CorrectComst;
            end
            bpnmConst=0.34;
            EnCylSpace=1;   % 2.3/2
            SInplain = obj.ScaledCylInplanePosition ;
            for i=1:n
                if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                    ZZrange=(obj.Zbase1(i):1:obj.Zbase2(i))';
                    CC=mod(ZZrange-correct(1),obj.period(1))*obj.period(2)*2*pi/obj.period(1);  %theta
                    %                  x0=EnCylSpace*obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=EnCylSpace*obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1) ; y0=EnCylSpace*SInplain(i,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*bpnmConst ,ones(size(x))];
                    %                 GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                else
                    ZZrange=(obj.Zbase1(i):1:obj.Zbase2(i))';
                    CC=mod(ZZrange-correct(2),obj.period(1))*obj.period(2)*2*pi/obj.period(1)  ; %theta
                    %                  x0=EnCylSpace*obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=EnCylSpace*obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1) ; y0=EnCylSpace*SInplain(i,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*bpnmConst ,ones(size(x))];
                    %                  GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                end
                HCoor{i}=LocalXYZ(:,1:3);
            end
        end  % end of HelixXYZGStap
        
        function HCoor=get.HelixXYZGBVec(obj)
            n=length(obj.Zbase1);        % Bvec, pointing to the axis
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=obj.CorrectComst;
            rr=[1,2] ;
            EnCylSpace=1;   % 2.3/2
            SInplain = obj.ScaledCylInplanePosition ;
            for i=1:n
                if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                    ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                    CC=mod(ZZrange-correct(1),obj.period(1))*obj.period(2)*2*pi/obj.period(1);  %theta
                    %                  x0=obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1) ; y0=EnCylSpace*SInplain(i,2);
                    [xr2,yr2] = pol2cart(CC,rr(2));
                    [xr1,yr1]= pol2cart(CC,rr(1));
                    x=xr2-xr1; y=yr2-yr1;
                    LocalXYZ=[x, y, ZZrange*0 ,zeros(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                    %                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;
                else
                    ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                    CC=mod(ZZrange-correct(2),obj.period(1))*obj.period(2)*2*pi/obj.period(1)  ; %theta
                    %                  x0=obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1) ; y0=EnCylSpace*SInplain(i,2);
                    [xr2,yr2] = pol2cart(CC,rr(2));
                    [xr1,yr1]= pol2cart(CC,rr(1));
                    x=xr2-xr1; y=yr2-yr1;
                    
                    LocalXYZ=[x, y, ZZrange*0 ,zeros(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                    %                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;
                end
                HCoor{i}=GlobalXYZ(:,1:3);
            end
        end     % end of    get.HelixXYZGBVec
        
        function obj=BundleCylinder(type,varargin)
            if type==1
                AZ=varargin{2};
                BZ=varargin{3};
                mmAZ=min(AZ);
                if mmAZ-21<0
                    AZ=AZ+21;
                    BZ=BZ+21;
                end
                obj.Zbase1=AZ;
                obj.Zbase2=BZ;
                obj.CylInplanePosition=varargin{4};
                if nargin==6
                    if ~isempty(varargin{5})
                        obj.Tol=  varargin{5};
                    end
                end
                
                Adj=zeros(length(obj.Zbase1),length(obj.Zbase1));
                for i=1:length(obj.Zbase1)
                    for j=1:length(obj.Zbase1)
                        if abs(norm(obj.CylInplanePosition(i,:)-obj.CylInplanePosition(j,:))-2*obj.CylRadius)<1e-3
                            H=[obj.Zbase1(i)  obj.Zbase2(i) obj.Zbase1(j)  obj.Zbase2(j)];
                            if max(H)-min(H)<= abs(obj.Zbase1(i)-obj.Zbase2(i))+ abs(obj.Zbase1(j)-obj.Zbase2(j))
                                Adj(i,j)=1; Adj(j,i)=1;
                            end
                        end
                    end
                end
                obj.CylAdjMat=Adj;     obj.AGroup=union(obj.AGroup,1);AB=1;    obj.BGroup=[];
                test=0;
                while 1
                    test=test+1;
                    switch AB
                        case 1
                            currentCyl=obj.AGroup;
                            NextCyl=[];
                            for i=1:length(currentCyl)
                                NextCyl=union(NextCyl,find(obj.CylAdjMat(:,currentCyl(i))==1) );
                            end
                            obj.BGroup=union(obj.BGroup,NextCyl);
                            AB=2;
                            
                        case 2
                            currentCyl=obj.BGroup;
                            NextCyl=[];
                            for i=1:length(currentCyl)
                                NextCyl=union(NextCyl,find(obj.CylAdjMat(:,currentCyl(i))==1) );
                            end
                            obj.AGroup=union(obj.AGroup,NextCyl);
                            AB=1;
                    end
                    if length(obj.AGroup)+length(obj.BGroup)==size(obj.CylAdjMat,1)
                        obj.AGroup=obj.AGroup';
                        obj.BGroup=obj.BGroup';
                        break
                    end
                    if test>500   %isolate cylinder
                        %                        IsolatedCyl= setdiff(setdiff(1:length(varargin{2}),obj.BGroup),obj.AGroup);
                        RemainCyl=union(obj.BGroup,obj.AGroup);
                        obj.Zbase1=obj.Zbase1(RemainCyl);obj.Zbase2=obj.Zbase2(RemainCyl);obj.CylInplanePosition=obj.CylInplanePosition(RemainCyl,:);
                        obj.CylAdjMat= obj.CylAdjMat(RemainCyl,RemainCyl);
                        NAGp=zeros(size(obj.AGroup));NBGp=zeros(size(obj.BGroup));
                        for i=1:length(NAGp)
                            NAGp(i)=find(obj.AGroup(i)==RemainCyl);
                        end
                        for i=1:length(NBGp)
                            NBGp(i)=find(obj.BGroup(i)==RemainCyl);
                        end
                        obj.AGroup=NAGp';obj.BGroup=NBGp';
                        break;
                    end
                end   %end of wihle
                %                  test
                %                 Tol=5;
                tolerance=obj.Tol;
                obj=obj.FindNarrowAdj;
                
                %------if Matlab Graph and Network Toolbox works;
                AdjM=obj.NarrowAdj;
                gg=graph(AdjM);
                obj.CylCategory=conncomp(gg);
                CC=cell(max(obj.CylCategory),2);
                for i=1:length(obj.CylCategory)
                    if ismember(i,obj.AGroup)
                        CC{  obj.CylCategory(i),1}=union(CC{  obj.CylCategory(i),1},i);
                    elseif ismember(i,obj.BGroup)
                        CC{  obj.CylCategory(i),2}=union(CC{  obj.CylCategory(i),2},i);
                    end
                end
                obj.CylABbyCategory=CC;
                
                
            elseif type==2   %union two parts into one
                MasterPart=varargin{1};
                SlavePart=varargin{2};
                AdjMNew=MasterPart.CylAdjMat;
                NewZ1=MasterPart.Zbase1;NewZ2=MasterPart.Zbase2;
                NewCylInplanePosition=MasterPart.CylInplanePosition;
                for i=1:length(SlavePart.Zbase1)
                    XY=SlavePart.CylInplanePosition(i,:);
                    Compare=zeros(size(NewCylInplanePosition,1),1);
                    for j=1:size(NewCylInplanePosition,1)
                        Compare(j)=sum(XY==NewCylInplanePosition(j,:))==2;
                    end
                    if sum(Compare)>=1  % has same inplane XY
                        indexofMaster=find(Compare==1);
                        iM=indexofMaster;
                        WW=[];
                        for k=1:length(indexofMaster)
                            WW=[WW NewZ1(indexofMaster(k))   NewZ2(indexofMaster(k))];
                        end
                        [ result ] = FindBelongRange( WW,[SlavePart.Zbase1(i)  SlavePart.Zbase2(i)]);
                        %                            height=sort([SlavePart.Zbase1(i)  SlavePart.Zbase2(i) NewZ1(indexofMaster)  NewZ2(indexofMaster)]);
                        %                            if max(height)-min(height)<= abs(SlavePart.Zbase1(i)-SlavePart.Zbase2(i))+abs(NewZ1(indexofMaster)-NewZ2(indexofMaster))
                        if ~isempty(result)
                            % Do union and check Adjacent
                            if length(result)==1    %just intersect with 1 cyl of Master
                                H=[SlavePart.Zbase1(i)  SlavePart.Zbase2(i) NewZ1(indexofMaster(result))  NewZ2(indexofMaster(result))];
                                NewZ1(indexofMaster(result))=min(H);NewZ2(indexofMaster(result))=max(H);
                                A=NewCylInplanePosition-ones(size( NewCylInplanePosition,1),1)*XY;
                                B=sqrt(A(:,1).*A(:,1)+A(:,2).*A(:,2));
                                NeiborCylInXY=find(B==4);
                                if ~isempty(NeiborCylInXY)
                                    for k=1:length(NeiborCylInXY)
                                        H2=sort([SlavePart.Zbase1(i)  SlavePart.Zbase2(i) NewZ1(NeiborCylInXY(k))  NewZ2(NeiborCylInXY(k))]);
                                        if max(H2)-min(H2)<= abs(SlavePart.Zbase1(i)-SlavePart.Zbase2(i))+abs(NewZ1(NeiborCylInXY(k))-NewZ2(NeiborCylInXY(k))) %contact
                                            AdjMNew(iM,NeiborCylInXY(k))=1; AdjMNew(NeiborCylInXY(k),iM)=1;
                                        end
                                    end
                                end
                                
                            else   %multiple intersect
                                H=[SlavePart.Zbase1(i)  SlavePart.Zbase2(i) NewZ1(indexofMaster(result))  NewZ2(indexofMaster(result))];
                                NewZ1(indexofMaster(result(1)))=min(H);NewZ2(indexofMaster(result(1)))=max(H);
                                DeleteCyl=iM(result(2:end));
                                AdjMNew(DeleteCyl,:)=[]; AdjMNew(:,DeleteCyl)=[];
                                NewZ1(DeleteCyl)=[];NewZ2(DeleteCyl)=[]; NewCylInplanePosition(DeleteCyl,:)=[];
                            end
                            %----
                        else  %create new cyl   and check Adjacent
                            NewZ1(end+1)=SlavePart.Zbase1(i);  NewZ2(end+1)=SlavePart.Zbase2(i);NewCylInplanePosition(end+1,:)=XY;
                            A=NewCylInplanePosition-ones(size( NewCylInplanePosition,1),1)*XY;
                            B=sqrt(A(:,1).*A(:,1)+A(:,2).*A(:,2));
                            NeiborCylInXY=find(B==4);
                            if ~isempty(NeiborCylInXY)
                                for k=1:length(NeiborCylInXY)
                                    H2=sort([SlavePart.Zbase1(i)  SlavePart.Zbase2(i) NewZ1(NeiborCylInXY(k))  NewZ2(NeiborCylInXY(k))]);
                                    if max(H2)-min(H2)<= abs(SlavePart.Zbase1(i)-SlavePart.Zbase2(i))+abs(NewZ1(NeiborCylInXY(k))-NewZ2(NeiborCylInXY(k))) %contact
                                        AdjMNew(length(NewZ1),NeiborCylInXY(k))=1; AdjMNew(NeiborCylInXY(k),length(NewZ1))=1;
                                    end
                                end
                            end
                        end
                        
                    else   %create new cyl in Master
                        NewZ1(end+1)=SlavePart.Zbase1(i);  NewZ2(end+1)=SlavePart.Zbase2(i);NewCylInplanePosition(end+1,:)=XY;
                        A=NewCylInplanePosition-ones(size( NewCylInplanePosition,1),1)*XY;
                        B=sqrt(A(:,1).*A(:,1)+A(:,2).*A(:,2));
                        NeiborCylInXY=find(B==4);
                        if ~isempty(NeiborCylInXY)
                            for k=1:length(NeiborCylInXY)
                                H2=sort([SlavePart.Zbase1(i)  SlavePart.Zbase2(i) NewZ1(NeiborCylInXY(k))  NewZ2(NeiborCylInXY(k))]);
                                if max(H2)-min(H2)<= abs(SlavePart.Zbase1(i)-SlavePart.Zbase2(i))+abs(NewZ1(NeiborCylInXY(k))-NewZ2(NeiborCylInXY(k))) %contact
                                    AdjMNew(length(NewZ1),NeiborCylInXY(k))=1; AdjMNew(NeiborCylInXY(k),length(NewZ1))=1;
                                end
                            end
                        end
                    end
                end
                
                if isa(obj,'BundleCylinderHC')
                    obj=BundleCylinderHC(1,[],NewZ1,NewZ2,NewCylInplanePosition);
                elseif  isa(obj,'BundleCylinderSQ')
                    obj=BundleCylinderSQ(1,[],NewZ1,NewZ2,NewCylInplanePosition);
                else
                    iGotError=11
                end
            end
        end  % end of BundleCylinder, constructor
        
        function  obj=FindNarrowAdj(obj)
            NAdj=obj.CylAdjMat;
            BAdj=obj.CylAdjMat;
            TAdj=obj.CylAdjMat;
            [Cu,Cv]=find(NAdj==1);
            for i=1:nnz(NAdj)
                uu=Cu(i);
                vv=Cv(i);
                if abs(obj.Zbase1(uu)-obj.Zbase1(vv))>obj.Tol  || abs(obj.Zbase2(uu)-obj.Zbase2(vv))>obj.Tol
                    NAdj(uu,vv)=0;
                end
                if abs(obj.Zbase1(uu)-obj.Zbase1(vv))>obj.Tol
                    BAdj(uu,vv)=0;
                end
                if abs(obj.Zbase2(uu)-obj.Zbase2(vv))>obj.Tol
                    TAdj(uu,vv)=0;
                end
            end
            obj.NarrowAdj=NAdj;
            obj.BottomAdj=BAdj;
            obj.TopAdj=TAdj;
        end % end of FindNarrowAdj
        
        function BaseIndex=findNeiborBaseIndexofCylinder(obj,QueryCylinders)
            % will be called by hyperB
            if length(QueryCylinders)>length(obj.Zbase1)
                BaseIndex=[];   return;
            end
            if isempty(QueryCylinders)
                QueryCylinders=1:length(obj.Zbase1) ;
            end
            
            
            
%             EliminateFromSides=6;   %unit base , available bridge points from two sides
            EliminateFromSides=16;   %unit base , available bridge points from two sides
              % use 6 for rotor
%             fprintf('In bundle class, EliminateFromSides = %i. Default =16 \n', EliminateFromSides);
            
            
            nQ=length(QueryCylinders);
            templateHC=obj.template;
            TT=[];
            BaseIndex=cell(nQ,1);
            for Clyi=1:nQ
                Cylinder=QueryCylinders(Clyi);
                ZZrange=[obj.Zbase1(Cylinder), obj.Zbase2(Cylinder)];
                
                repNum=floor(ZZrange(2)/obj.period(1));
                TT=[];
                
                for k=1:repNum
                    TT=union(TT,templateHC+obj.period(1)*k*ones(1,length(templateHC)));
                end
                TT(TT<ZZrange(1)+EliminateFromSides)=[];
                TT(TT>ZZrange(2)-EliminateFromSides)=[];
                BaseIndex{Clyi}=TT;
            end
        end  % end of findNeiborBaseIndexofCylinder
        
        function ZZ1=get.Z1(obj)
            ZZ1=0.34*obj.Zbase1;
        end
        function ZZ2=get.Z2(obj)
            ZZ2=0.34*obj.Zbase2;
        end
        function TMatrix=get.TransformMatrix(obj)
            TMatrix=obj.TransformMatrix2;
        end
        
        function  plotcylinder(obj, RTable )
            figure(555) ; clf ; hold on;
            SInplain = obj.ScaledCylInplanePosition ;
            for k= 1 : size(obj.CylInplanePosition ,1)
                x= SInplain(k,1)  ;
                y = SInplain(k,2)  ;
                plot3([x;x], [y;y] , [obj.Z1(k) ; obj.Z2(k)  ], 'k') ;
            end
            
            HCoor= HelixXYZGStapNoPMLocal(obj,1,1) ;
            for k= 1 :size(HCoor , 1)
                plot3(HCoor{k}(:,1) , HCoor{k}(:,2) ,HCoor{k}(:,3)    , 'b') ;
            end
            
            [ HClattice ] = findHClattice( 1 ,[20 20]) ;
            scatter( HClattice.HCcenter(:,1) , HClattice.HCcenter(:,2)  ,'r' )  ;
            
            ranRefCyl =randi(9)+14   ;
            %             RefXY=SInplain(RTable(ranRefCyl,2),:);
            RefXY =  obj.CylInplanePosition(RTable(ranRefCyl,2),:)  ;
            RefColRow =  RTable(ranRefCyl, 6:7)  ;
            HCMapping= findHClatticeMapping( 1 ,[50 50]) ;
            [~,ind] =ismember(RefXY ,  HCMapping(:,1:2) ,'rows') ;
            
            RefCyldColRow_shift =  RefColRow- HCMapping(ind,3:4) ; %
            %                 scatter( HCMapping(:,1),HCMapping(:,2)  ,64,'g.' );
            %                 text( HCMapping(:,1),HCMapping(:,2) ,num2str(HCMapping(:,3:4) ));
            %             return
            
            scatter(RefXY(1), RefXY(2) , 84,'s') ;
            text(RefXY(1), RefXY(2) ,num2str(RefColRow)) ;
            
            %             for cc= 1 :10
            %                 RandColRow = [5+randi(6), 3+randi(6) ]  ;   % assume as imput
            % %                 RandColRow - RefCyldColRow_shift   ;
            %                 [~,ind2] =ismember( RandColRow - RefCyldColRow_shift ,  HCMapping(:,3:4) ,'rows')   ;
            %
            %                 RandCyldXY =  HCMapping(ind2,1:2)  ; %
            %
            %                 scatter(RandCyldXY(1),RandCyldXY(2) ,'m' );
            %                 text( RandCyldXY(1),RandCyldXY(2) ,num2str(RandColRow ));
            %             end
        end
        
        
        
        
        function HCoor=findHelixWithR(obj,rr,QCylinder,QBasesinCell)
            if length(QCylinder)~=length(QBasesinCell)
                HCoor=[];
                return
            end
            nQ=length(QBasesinCell);      %
            HCoor=cell(nQ,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=obj.CorrectComst;   % check with Cadnano 12/6/2016
            SInplain = obj.ScaledCylInplanePosition ;
            for i=1:nQ
                Cly=QCylinder(i);
                ZZrange=QBasesinCell{i};
                if size(ZZrange,2)~=1
                    ZZrange=ZZrange';
                end
                if ~xor(ismember(Cly,obj.AGroup) ,obj.AGroupGoUp)
                    CC=mod(ZZrange-correct(1),obj.period(1))*obj.period(2)*2*pi/obj.period(1);  %theta
                    %                     x0=obj.CylInplanePosition(Cly,1)*obj.SR_2p5; y0=obj.CylInplanePosition(Cly,2)*obj.SR_2p5;
                    x0=SInplain(Cly,1); y0=SInplain(Cly,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                else
                    CC=mod(ZZrange-correct(2),obj.period(1))*obj.period(2)*2*pi/obj.period(1)  ; %theta
                    x0=obj.CylInplanePosition(Cly,1)*obj.SR_2p5; y0=obj.CylInplanePosition(Cly,2)*obj.SR_2p5;
                    x0=SInplain(Cly,1); y0=SInplain(Cly,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                end
                HCoor{i}=GlobalXYZ(:,1:3);
            end
        end % end of findHelixWithR
        
        function HCoor=get.HelixXYZG(obj)
            n=length(obj.Zbase1);      % orientation haven't check
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=obj.CorrectComst;
            rr=0.85;    %Backbone cylindric location
            EnCylSpace=1;   % 2.3/2
            SInplain = obj.ScaledCylInplanePosition ;
            %             [ismember(9,obj.AGroup) ,ismember(8,obj.AGroup) ]
            for i=1:n
                if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                    ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                    CC=mod(ZZrange-correct(1),obj.period(1))*obj.period(2)*2*pi/obj.period(1);  %theta
                    %                  x0=EnCylSpace*obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=EnCylSpace*obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1); y0=EnCylSpace*SInplain(i,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                    %                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;
                else
                    ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                    CC=mod(ZZrange-correct(2),obj.period(1))*obj.period(2)*2*pi/obj.period(1)  ; %theta
                    %                  x0=EnCylSpace*obj.CylInplanePosition(i,1)*obj.SR_2p5; y0=EnCylSpace*obj.CylInplanePosition(i,2)*obj.SR_2p5;
                    x0=EnCylSpace*SInplain(i,1); y0=EnCylSpace*SInplain(i,2);
                    [x,y] = pol2cart(CC,rr);
                    LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))];
                    GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                    %                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;
                end
                HCoor{i}=GlobalXYZ(:,1:3);
            end
        end
        
        function GCoor=get.CylinderXYZGlobal(obj)
            n=length(obj.Zbase1);
            GCoor=zeros(n,6) ;   %[x1 y1 z1   x2 y2 z2]  in unit of nm
            SInplain = obj.ScaledCylInplanePosition ;
            for i=1:n
                LocalXYZ=[SInplain(i,1:2)  obj.Z1(i)  SInplain(i,1:2) obj.Z2(i)]; %2.5nm
                %                                 LocalXYZ=[obj.CylInplanePosition(i,1:2)*obj.SR_2p5  obj.Z1(i)  obj.CylInplanePosition(i,1:2)*obj.SR_2p5 obj.Z2(i)]; %2.5nm
                
                %                 LocalXYZ=[obj.CylInplanePosition(i,1:2)  obj.Z1(i)
                %                 obj.CylInplanePosition(i,1:2) obj.Z2(i)]; % 2nm spacing.
                %                 MagicDNA2
                
                R=obj.TransformMatrix2(1:3,1:3);
                T=obj.TransformMatrix2(1:3,4);
                
                S1=R*[LocalXYZ(1);LocalXYZ(2);LocalXYZ(3)] +T  ;
                S2=R*[LocalXYZ(4);LocalXYZ(5);LocalXYZ(6)] +T ;
                X1=S1(1);  Y1=S1(2);  ZZ1=S1(3);
                X2=S2(1);  Y2=S2(2);  ZZ2=S2(3);
                %                 [X1,Y1,ZZ1] = transformPointsForward(obj.TransformMatrix2,LocalXYZ(1),LocalXYZ(2),LocalXYZ(3));
                %                 [X2,Y2,ZZ2] = transformPointsForward(obj.TransformMatrix2,LocalXYZ(4),LocalXYZ(5),LocalXYZ(6));
                GCoor(i,:)=[X1,Y1,ZZ1,X2,Y2,ZZ2];
            end
        end
        
        
        
    end
    
end

