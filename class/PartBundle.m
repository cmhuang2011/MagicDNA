classdef PartBundle %< handle
   % single Part contains cylinders information
    %------update 08222016
    
    %last edit: 11202016
    properties
        
        Zbase1=[];     %unit base
        Zbase2=[];
        CylInplanePosition=[];
        CylAdjMat=[];
        NarrowAdj=[];
        BottomAdj=[];
        TopAdj=[];
        
        Tol=5;
        maxTol=40;
        CylRadius=1;
        ssScafOver=5;
        
        AGroup=[];
        BGroup=[];
        CylCategory=[];
        CylABbyCategory=[];
        AGroupGoUp=1;   %default value
        
        %------
        Rx=0;  %rad
        Ry=0;
        Rz=0;
        LocalOinG=[0 0 0];
        PhaseAngleByBase=3;
        %----
        ExternalXoverAsFB=[];         %external xover position as forced-connection of each bundle
        TransformMatrix2=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
        SimulateTransMFromTM2=[];
    end
    
    properties (Dependent)
       Z1;     %in nm
       Z2;
       CylinderXYZGlobal;
       TransformMatrix;  %4 by 4
       HelixXYZG;
       HelixXYZGBound30 ;
       HelixXYZGBVec;
       HelixXYZGStap;
    end
    
    methods
        function ZZ1=get.Z1(obj)
            ZZ1=0.34*obj.Zbase1;            
        end
        function ZZ2=get.Z2(obj)
            ZZ2=0.34*obj.Zbase2;            
        end
        function TMatrix=get.TransformMatrix(obj)
            TMatrix=obj.TransformMatrix2;
%             TRx=obj.rotx(obj.Rx);
%             TRy=obj.roty(obj.Ry);
%             TRz=obj.rotz(obj.Rz);
%             TR=TRx*TRy*TRz;
%             TMatrix2=[TR ,[0; 0; 0];obj.LocalOinG, 1]; 
%             TMatrix=affine3d(TMatrix2);
        end
        
        function GCoor=get.CylinderXYZGlobal(obj)
            n=length(obj.Zbase1);
            GCoor=zeros(n,6) ;   %[x1 y1 z1   x2 y2 z2]  in unit of nm
            for i=1:n
                LocalXYZ=[obj.CylInplanePosition(i,1:2)  obj.Z1(i)  obj.CylInplanePosition(i,1:2) obj.Z2(i)]; 
                
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
        
        function HCoor=get.HelixXYZG(obj)
            n=length(obj.Zbase1);     
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=[5.5 ,0]-0.6;
            rr=0.85;    %Backbone cylindric location
            EnCylSpace=1;   % 2.3/2
%             [ismember(9,obj.AGroup) ,ismember(8,obj.AGroup) ]
             for i=1:n
                 if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                 ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                 CC=mod(ZZrange-correct(1),32)*3*2*pi/32;  %theta
                 x0=EnCylSpace*obj.CylInplanePosition(i,1); y0=EnCylSpace*obj.CylInplanePosition(i,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))]; 
                GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;     
                 else
                 ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                 CC=mod(ZZrange-correct(2),32)*3*2*pi/32  ; %theta
                 x0=EnCylSpace*obj.CylInplanePosition(i,1); y0=EnCylSpace*obj.CylInplanePosition(i,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))]; 
                 GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;                         
                 end              
                HCoor{i}=GlobalXYZ(:,1:3);
             end
        end
        function HCoor=get.HelixXYZGBound30(obj)
            extr=30;
            n=length(obj.Zbase1);    
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=[5.5 ,0]-0.6;
            rr=0.85;    %Backbone cylindric location
            EnCylSpace=1;   % 2.3/2
%             [ismember(9,obj.AGroup) ,ismember(8,obj.AGroup) ]
             for i=1:n
                 if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                 ZZrange=(obj.Zbase1(i)-extr:1:obj.Zbase2(i)+extr)';
                 CC=mod(ZZrange-correct(1),32)*3*2*pi/32;  %theta
                 x0=EnCylSpace*obj.CylInplanePosition(i,1); y0=EnCylSpace*obj.CylInplanePosition(i,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))]; 
                GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;     
                 else
                 ZZrange=(obj.Zbase1(i)-extr:1:obj.Zbase2(i)+extr)';
                 CC=mod(ZZrange-correct(2),32)*3*2*pi/32  ; %theta
                 x0=EnCylSpace*obj.CylInplanePosition(i,1); y0=EnCylSpace*obj.CylInplanePosition(i,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))]; 
                 GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;                         
                 end              
                HCoor{i}=GlobalXYZ(:,1:3);
             end
        end        
        
        
        
        
       function HCoor=get.HelixXYZGStap(obj)
            n=length(obj.Zbase1);      % orientation haven't check
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=[0 ,5.5]-0.6;
            bpnmConst=0.34;
            rr=0.85;
%             rr=0;  %used to draw staple cylinder in 3D
              EnCylSpace=1;   % 2.3/2
%             [ismember(9,obj.AGroup) ,ismember(8,obj.AGroup) ]
             for i=1:n
                 if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                 ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                 CC=mod(ZZrange-correct(1),32)*3*2*pi/32;  %theta
                 x0=EnCylSpace*obj.CylInplanePosition(i,1); y0=EnCylSpace*obj.CylInplanePosition(i,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*bpnmConst ,ones(size(x))]; 
                GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;     
                 else
                 ZZrange=(obj.Zbase1(i)-10:1:obj.Zbase2(i)+10)';
                 CC=mod(ZZrange-correct(2),32)*3*2*pi/32  ; %theta
                 x0=EnCylSpace*obj.CylInplanePosition(i,1); y0=EnCylSpace*obj.CylInplanePosition(i,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*bpnmConst ,ones(size(x))]; 
                 GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;                         
                 end              
                HCoor{i}=GlobalXYZ(:,1:3);
             end
       end
        
       
       
        function HCoor=get.HelixXYZGBVec(obj)
            extra=30;
            n=length(obj.Zbase1);        % Bvec, pointing to the axis 
            HCoor=cell(n,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=[5.5 ,0]-0.6;
            rr=[1,2] ;
%             [ismember(9,obj.AGroup) ,ismember(8,obj.AGroup) ]
             for i=1:n
                 if ~xor(ismember(i,obj.AGroup) ,obj.AGroupGoUp)
                 ZZrange=(obj.Zbase1(i)-extra:1:obj.Zbase2(i)+extra)';
                 CC=mod(ZZrange-correct(1),32)*3*2*pi/32;  %theta
                 x0=obj.CylInplanePosition(i,1); y0=obj.CylInplanePosition(i,2);
                 [xr2,yr2] = pol2cart(CC,rr(2));
                  [xr1,yr1]= pol2cart(CC,rr(1));
                 x=xr2-xr1; y=yr2-yr1;
                LocalXYZ=[x, y, ZZrange*0 ,zeros(size(x))]; 
                GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix2.T;     
                 else
                 ZZrange=(obj.Zbase1(i)-extra:1:obj.Zbase2(i)+extra)';
                 CC=mod(ZZrange-correct(2),32)*3*2*pi/32  ; %theta
                 x0=obj.CylInplanePosition(i,1); y0=obj.CylInplanePosition(i,2);
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
        
        
 
                
        
        function HCoor=findHelixQ(obj,QCylinder,QBasesinCell)
            if length(QCylinder)~=length(QBasesinCell)
               HCoor=[]; 
               return 
            end
             rr=0.85;

            
            nQ=length(QBasesinCell);      % 
            HCoor=cell(nQ,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=[5.5 ,0]-0.6;   % check with Cadnano 12/6/2016
%             Over=0;
             for i=1:nQ
%                  i;
                 Cly=QCylinder(i);
                 ZZrange=QBasesinCell(i);
                 if length(ZZrange)==106
                     sdfsf=34234;
                 end
                 if ~xor(ismember(Cly,obj.AGroup) ,obj.AGroupGoUp)
                 CC=mod(ZZrange-correct(1),32)*3*2*pi/32;  %theta
                 x0=obj.CylInplanePosition(Cly,1); y0=obj.CylInplanePosition(Cly,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))]; 
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix;    
                GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                  GlobalXYZ=LocalXYZ*obj.TransformMatrix.T;
%                  %before 2/1/2017
                 else
                 CC=mod(ZZrange-correct(2),32)*3*2*pi/32  ; %theta
                 x0=obj.CylInplanePosition(Cly,1); y0=obj.CylInplanePosition(Cly,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))]; 
%                 size(LocalXYZ)
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix.T;                 
%                  GlobalXYZ=LocalXYZ*obj.TransformMatrix.T;      
                   GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                 end              
                HCoor{i}=GlobalXYZ(:,1:3);
             end                     
        end
        
                function HCoor=findHelixWithR(obj,rr,QCylinder,QBasesinCell)
            if length(QCylinder)~=length(QBasesinCell)
               HCoor=[]; 
               return 
            end
         
            
            nQ=length(QBasesinCell);      % 
            HCoor=cell(nQ,1) ;   %[x1 y1 z1   ]  in unit of nm
            correct=[5.5 ,0]-0.6;   % check with Cadnano 12/6/2016
%             Over=0;
             for i=1:nQ
%                  i;
                 Cly=QCylinder(i);
                 ZZrange=QBasesinCell{i};
                 if size(ZZrange,2)~=1
                     ZZrange=ZZrange';
                 end
                 
                 
%                  if length(ZZrange)==106
%                      sdfsf=34234
%                  end
                 if ~xor(ismember(Cly,obj.AGroup) ,obj.AGroupGoUp)
                 CC=mod(ZZrange-correct(1),32)*3*2*pi/32;  %theta
                 x0=obj.CylInplanePosition(Cly,1); y0=obj.CylInplanePosition(Cly,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))]; 
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix;    
% % obj.TransformMatrix2
% % LocalXYZ'
% % if isempty(LocalXYZ)
% %     sdfsf=324
% % end
% obj.TransformMatrix2
% LocalXYZ
% if size(LocalXYZ,2)>300
%     dfgdg=345
% end

                GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
%                  GlobalXYZ=LocalXYZ*obj.TransformMatrix.T;
%                  %before 2/1/2017
                 else
                 CC=mod(ZZrange-correct(2),32)*3*2*pi/32  ; %theta
                 x0=obj.CylInplanePosition(Cly,1); y0=obj.CylInplanePosition(Cly,2);
                 [x,y] = pol2cart(CC,rr);                
                LocalXYZ=[x+x0*ones(size(x)), y+y0*ones(size(y)), ZZrange*0.34 ,ones(size(x))]; 
%                 size(LocalXYZ)
%                 GlobalXYZ=LocalXYZ*obj.TransformMatrix.T;                 
%                  GlobalXYZ=LocalXYZ*obj.TransformMatrix.T;      
                   GlobalXYZ=transpose(obj.TransformMatrix2*LocalXYZ' );
                 end              
                HCoor{i}=GlobalXYZ(:,1:3);
             end                     
        end
        
        function BaseIndex=findNeiborBaseIndexofCylinder(obj,QueryCylinders)
              if length(QueryCylinders)>length(obj.Zbase1)
                 BaseIndex=[];   return;
              end
            if isempty(QueryCylinders)
                QueryCylinders=1:length(obj.Zbase1) ;
            end
          
            EliminateFromSides=25;   %unit base , available bridge points from two sides
            nQ=length(QueryCylinders);  
            template=[2 ,3, 4, 5,7 8,10 ,11, 12,13,15 16,18 ,19,20 21,23, 24,26,27,28,29,31 ,32];
            TT=[];
            BaseIndex=cell(nQ,1);
             for Clyi=1:nQ
                 Cylinder=QueryCylinders(Clyi);
                 ZZrange=[obj.Zbase1(Cylinder), obj.Zbase2(Cylinder)];
                 
                 repNum=floor(ZZrange(2)/32);
                        TT=[];
                 
                 for k=1:repNum
                 TT=union(TT,template+32*k*ones(1,length(template)));
                 end
                 TT(TT<ZZrange(1)+EliminateFromSides)=[];
                 TT(TT>ZZrange(2)-EliminateFromSides)=[];
                 BaseIndex{Clyi}=TT;
             end
        end
        
        
        function obj= changeTol(obj, Cval)
            OriTol=obj.Tol;
            OriZ1=obj.Zbase1;
            OriZ2=obj.Zbase2;
            OXY=obj.CylInplanePosition;
            if OriTol+Cval>=obj.maxTol
               return 
            end                     
            if Cval~=0
               NTol=OriTol+Cval;
               obj=PArtBundle(1,[],OriZ1,OriZ2,OXY,NTol);   
               obj.Tol=NTol;
            end
            
        end   %end of change tol
            
        function obj=PartBundle(type,varargin)
            if type==1
%                 obj.CylAdjMat=sparse(varargin{1});
%                 obj.CylAdjMat=[];
                AZ=varargin{2};
                BZ=varargin{3};
                mmAZ=min(AZ);
%                 if mmAZ-32<0
%                 AZ=AZ+32;
%                 BZ=BZ+32;    
%                 end
                    
                    
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
                obj.CylAdjMat=Adj;                
                 obj.AGroup=union(obj.AGroup,1);AB=1;
                 obj.BGroup=[];
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
                   if length(obj.AGroup)+length(obj.BGroup)==size(obj.CylAdjMat,1);
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
                obj=PArtBundle(1,[],NewZ1,NewZ2,NewCylInplanePosition);
            end
        end
        


        function path=findsecondPath(obj,stem1,side,SP2,EP2)
            repeated=0;
            n=size(obj.CylAdjMat,1);
            parent=zeros(n,1);
            visited=zeros(n,1);
            TAdj=obj.TopAdj;
            BAdj=obj.BottomAdj;
            S1=stem1.getPath;
            
           if side==1;AB=true; else AB=false;end
           for i=1:length(S1)-1     %eliminate used edges from Stem 1
              if  xor(AB,mod(i,2)==0)
               TAdj(S1(i),S1(i+1))=0;
               TAdj(S1(i+1),S1(i))=0;   
              else
               BAdj(S1(i),S1(i+1))=0;
               BAdj(S1(i+1),S1(i))=0;  
              end
           end
            currentNode=SP2;
            randomPath=currentNode;
            visited(currentNode)=1;
            AB2=xor(ismember(S1(1),obj.AGroup),xor(AB,ismember(SP2(1),obj.AGroup)));
            while 1               
                          switch AB2
                             case 1
                                neighbor=find(TAdj(currentNode,:) == 1); 
                                if ~isempty(neighbor)
                                    child=neighbor(randi(length(neighbor)));
                                    TAdj(currentNode,child)=0;
                                    TAdj(child,currentNode)=0;
                                    if visited(child) > repeated
                                        continue;
                                    else
                                        parent(child)=currentNode;
                                        visited(child)=visited(child)+1;
                                        randomPath(end+1)=child;
                                        currentNode=child;                   
                                    end                    
                                else
                                    randomPath(end)=[];   %back to previous
                                    currentNode=parent(currentNode);
                                end
                               AB2=~AB2 ;
                             case 0
                                neighbor=find(BAdj(currentNode,:) == 1); 
                                if ~isempty(neighbor)
                                    child=neighbor(randi(length(neighbor)));
                                    BAdj(currentNode,child)=0;
                                    BAdj(child,currentNode)=0;
                                    if visited(child) > repeated
                                        continue;
                                    else
                                        parent(child)=currentNode;
                                        visited(child)=visited(child)+1;
                                        randomPath(end+1)=child;
                                        currentNode=child;                   
                                    end                    
                                else
                                    randomPath(end)=[];   %back to previous
                                    currentNode=parent(currentNode);
                                end
                                    AB2=~AB2;
                          end  %end of switch
                     if currentNode == EP2  %|| currentNode==0
                                                                           
               
                         break;  
                     end             
            end 
                path=Path([]);
                for i=1:length(randomPath)-1
                    path=path.addEdge([randomPath(i) randomPath(i+1)]); 
                end
           
   
        end
        
        function [path]=findRandomPathSS2(obj,startNode,endNode,repeated,side)
            %find a random path between two paths
            n=size(obj.CylAdjMat,1);
            parent=zeros(n,1);
            visited=zeros(n,1);
%             adj=obj.CylAdjMat;
            TAdj=obj.TopAdj;
            BAdj=obj.BottomAdj;
            currentNode=startNode;
            randomPath=currentNode;
            visited(currentNode)=1;
                    if side==1;
                    AB=true;
                    else 
                    AB=false;
                    end
            while 1               
                          switch AB
                             case 1
                                neighbor=find(TAdj(currentNode,:) == 1); 
                                if ~isempty(neighbor)
                                    child=neighbor(randi(length(neighbor)));
                                    TAdj(currentNode,child)=0;
                                    TAdj(child,currentNode)=0;
                                    if visited(child) > repeated
                                        continue;
                                    else
                                        parent(child)=currentNode;
                                        visited(child)=visited(child)+1;
                                        randomPath(end+1)=child;
                                        currentNode=child;                   
                                    end                    
                                else
                                    randomPath(end)=[];   %back to previous
                                    currentNode=parent(currentNode);
                                end
                               AB=~AB; 
                             case 0
                                neighbor=find(BAdj(currentNode,:) == 1); 
                                if ~isempty(neighbor)
                                    child=neighbor(randi(length(neighbor)));
                                    BAdj(currentNode,child)=0;
                                    BAdj(child,currentNode)=0;
                                    if visited(child) > repeated
                                        continue;
                                    else
                                        parent(child)=currentNode;
                                        visited(child)=visited(child)+1;
                                        randomPath(end+1)=child;
                                        currentNode=child;                   
                                    end                    
                                else
                                    randomPath(end)=[];   %back to previous
                                    currentNode=parent(currentNode);
                                end
                                    AB=~AB;
                          end  %end of switch
                     if currentNode == endNode
                         break;
                     end               
            end            
                path=Path([]);
                for i=1:length(randomPath)-1
                    path=path.addEdge([randomPath(i) randomPath(i+1)]); 
                end
            
        end
        
        
        
        

        function  TF=checkstem(obj,stem,whichside)          
                    path=stem.getPath;
                    nStem=length( stem.getPath);
                    checkresult=zeros(1,nStem-1);
                    for i=1:nStem-1
                       Pcyl=path(i); Ncyl=path(i+1);
                       switch mod(i+whichside,2)
                           case 0
                               if abs(obj.Zbase2(Pcyl)-obj.Zbase2(Ncyl))<=obj.Tol
                                   checkresult(i)=1;
                               end
                           case 1
                               if abs(obj.Zbase1(Pcyl)-obj.Zbase1(Ncyl))<=obj.Tol
                                   checkresult(i)=1;
                               end
                       end                       
                    end
                    
                    if nnz(checkresult)==length(checkresult);TF=1;else TF=0;end     
              
        end
        
        function  [Doable,PairList,Unpair]=checkpairable(obj,varargin)
            adjM=obj.NarrowAdj;
            path=[];
            for i=1:length(varargin)           
                stem=varargin{i};
                path=union(stem.getPath,path);
            end
            adjM(path,:)=0;adjM(:,path)=0;
            PairList=zeros(ceil((size(obj.CylAdjMat,1)-length(path))/2) ,2);
            nP=0;
            N=0;
            Doable=false;
            Count=0;           
         while N<5000
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
                potent=find(adjM(tt(1),:)==1);                
                wa=potent(randi([1 length(potent)]));
                if abs(obj.Zbase1(wa)-obj.Zbase1(tt))<=obj.Tol && abs(obj.Zbase2(wa)-obj.Zbase2(tt))<=obj.Tol
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
                 xx=xx(randi([1 length(xx)]));
                  if abs(obj.Zbase1(yy)-obj.Zbase1(xx))<=obj.Tol && abs(obj.Zbase2(yy)-obj.Zbase2(xx))<=obj.Tol
                 nP=nP+1;
                 PairList(nP,1)=yy;
                 PairList(nP,2)=xx;
                 adjM(yy,:)=0;adjM(:,yy)=0;
                 adjM(xx,:)=0;adjM(:,xx)=0;    
                 Anything=1;
                  end
                 end
             end
             if Anything==0;
                 Count=Count+1;
             end

             if nnz(PairList)==size(PairList,1)*size(PairList,2) ;
                 Doable=true;    break; end
         end   %end of while
            used=PairList(:);
            Unpair=setdiff(setdiff(union(obj.AGroup,obj.BGroup),used),path);
            PairList(sum(PairList,2)==0,: )=[];
        end
        
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
            
           
        end
        
        function  same=isequaltopart(obj, targetpart)           
         same=true  ; 
         if length(obj.Zbase1)~=length(targetpart.Zbase1)  || length(obj.Zbase2)~=length(targetpart.Zbase2) ||size(obj.CylInplanePosition,1)~=size(targetpart.CylInplanePosition,1)
             same=false;return;
         end
         
         
%          c1=nnz(obj.CylAdjMat==targetpart.CylAdjMat)~=size(targetpart.CylAdjMat,1)*size(targetpart.CylAdjMat,2);
         c2=nnz(obj.Zbase1==targetpart.Zbase1)~=size(targetpart.Zbase1,1)*size(targetpart.Zbase1,2);
         c3=nnz(obj.Zbase2==targetpart.Zbase2)~=size(targetpart.Zbase2,1)*size(targetpart.Zbase2,2);
         c4=nnz( obj.CylInplanePosition==targetpart.CylInplanePosition)~=size(targetpart.CylInplanePosition,1)*size(targetpart.CylInplanePosition,2);
           if  c2 || c3 || c4
           same=false;
           end
           
        end
        function   [Result,affectedcyl]=isintersect(obj,testpart)
            Result=false;       
            xy=obj.CylInplanePosition;
            cyl=[];
            affectedcyl=[];
            for k=1:length(testpart.Zbase1)  %if testpart have mulit cylinder
            
                for i=1:size(xy,1)
                    if sum(xy(i,:)==testpart.CylInplanePosition(k,:))==2
                      cyl=union(cyl,i);
                    end
                end

                for j=1:length(cyl)
                   h= [obj.Zbase1(cyl(j)) obj.Zbase2(cyl(j)) testpart.Zbase1(k) testpart.Zbase2(k)];
                    if max(h)-min(h)<= abs(h(1)-h(2))+abs(h(3)-h(4)) 
                        Result=true;
                        affectedcyl=union(affectedcyl,cyl(j));
%                         return
                    end
                end
                
            end
        end
        
        function Cyl=ReturnCylInfo(obj,cylindex)
            Cyl.NZ1=obj.Zbase1(cylindex);Cyl.NZ2=obj.Zbase2(cylindex);
            Cyl.NCylInplanePosition=obj.CylInplanePosition(cylindex,:);                                
        end
        
        function  rr=ReturnCylIndex(obj,Zbase1,Zbase2,XY)
            if size(XY,1)~=length(Zbase1) || size(XY,1)~=length(Zbase2)
                rr=0;
                return
            end
            rr=[];
           for i=1:length(Zbase1)
            iZ1=find(obj.Zbase1==Zbase1(i));
            iZ2=find(obj.Zbase2==Zbase2(i)) ;
            iXY=intersect(find(obj.CylInplanePosition(:,1)==XY(i,1)),find(obj.CylInplanePosition(:,2)==XY(i,2)));
            CC=intersect(intersect(iZ1,iZ2),iXY);
            if isempty(CC);CC=0;end
           rr= [rr CC];
           end
        end
        
        function    [outputgood,choices,Zbase1,Zbase2]=AddExtraCyl(obj,cylindex,samedir)
            NPCyl.Zbase1=obj.Zbase1(cylindex);  NPCyl.Zbase2=obj.Zbase2(cylindex);
            NPCyl.CylInplanePosition=obj.CylInplanePosition(cylindex,:);
           Categ=obj.CylCategory(cylindex);
            Allgroup=obj.CylCategory==Categ;
            XYexistdata=obj.CylInplanePosition(Allgroup,:);
           NewXY=[999 999];
           for i=1:size(XYexistdata,1)
              XYi= XYexistdata(i,:);  
              EWSN=[XYi+[4 0];XYi+[-4 0];XYi+[0 4];XYi+[0 -4]];
              NewXY=union(NewXY,EWSN,'rows');
           end
           NewXY=setdiff(NewXY,[ 999 999],'rows');
            NewXY=setdiff(NewXY,XYexistdata,'rows');
            u=abs(NewXY(:,1)-NPCyl.CylInplanePosition(1))/4;
            v=abs(NewXY(:,2)-NPCyl.CylInplanePosition(2))/4;w=u+v;
           NewXYS= NewXY(mod(w+samedir,2)==1,:);
           order=randperm(size(NewXYS,1));
           Yesinter=zeros(size(NewXYS,1),1);
           for i=1:size(NewXYS,1)
              test=PArtBundle(1,[],NPCyl.Zbase1 ,NPCyl.Zbase2,NewXYS(order(i),:));
              [Yesinter(order(i)),~]=obj.isintersect(test);              
           end
           if sum(Yesinter)<length(Yesinter)
               outputgood=1;   %find good adding cylinders 
               choices=NewXYS(Yesinter==0,:);
           else
               outputgood=0;    %adding these will casuse interfere, requiring further dealing.
               choices=NewXYS;
           end
           Zbase1=NPCyl.Zbase1;Zbase2=NPCyl.Zbase2;
            %
        end
        
        function  Result=FindRoute(obj,StartP,EndP,side)            
         Doable=false;
         trial=0;
         SaveNoOfUnpair=union(obj.AGroup,obj.BGroup);
         while trial<=1000;
             Teststem=obj.findRandomPathSS2(StartP,EndP,0,side);
              [Doable,PairList,Unpair]=obj.checkpairable(Teststem);
              if length(Unpair)<length(SaveNoOfUnpair) %|| RENEW==1
              SaveNoOfUnpair=Unpair;
              SaveStem=Teststem;
              SavePair=PairList;
              end

             if  (Doable==1   ||  trial>=800)
                 break;
             end
             trial=trial+1;

         end        
         Result.stem=SaveStem;Result.Pair=SavePair;Result.SaveNoOfUnpair=SaveNoOfUnpair;
         Result.Doable=Doable;
        end   %End of FindRoute
        
        function [CycleConnectTo,SavePair]= FindBridgeBtwPair(obj,TempResult)
            SavePair=TempResult.Pair;
            SaveStem=TempResult.stem;
            stemcyl= SaveStem.getPath;
            AdjBind=zeros(size(SavePair,1)+1,size(SavePair,1)+1);
            [u,v]=find(obj.CylAdjMat==1);
            for i=1:length(u)
              if ismember(u(i), SavePair)==1 && ismember(v(i), SavePair)==1 
                  [x1,~]=find(u(i)==SavePair); [y1,~]=find(v(i)==SavePair);
                  AdjBind(1+x1,1+y1)=1;AdjBind(1+y1,1+x1)=1;
              end
               if ismember(u(i), stemcyl)==1 && ismember(v(i), SavePair)==1 
                   x1=1;[y1,~]=find(v(i)==SavePair);
                   AdjBind(x1,1+y1)=1;AdjBind(1+y1,x1)=1;
               end 
            end

            for i=1:size(AdjBind,1);AdjBind(i,i)=0;end
            gAdjB=graph(AdjBind);
            [~,pred] = minspantree(gAdjB);
            if  sum(isnan(pred))>0
               fsdfsf=23313; 
            end
            
            CycleConnectTo=zeros(size(SavePair,1),1);
            for i=1:size(SavePair,1)
                if pred(i+1)==1
                precyls=stemcyl;
                else
%                     i
%                     pred
                 precyls=  SavePair(pred(i+1)-1,:);
                end    
                C1=SavePair(i,1);C2=SavePair(i,2);
                if sum(obj.CylAdjMat(C1,precyls))>0
                    ind=find(obj.CylAdjMat(C1,precyls)==1);ind= ind(randi(length(ind)));
                    CycleConnectTo(i)=precyls(ind);
                else
                    SavePair(i,2)=C1; SavePair(i,1)=C2;
                    ind=find(obj.CylAdjMat(C2,precyls)==1);
                    ind= ind(randi(length(ind)));
                    CycleConnectTo(i)=precyls(ind);
                end
            end     
        end %end of findbridgebtwpair
        
        function  [obj,AffiVar]=FindRTable(obj,side,StartP)
                AffiVar.InplaneXYRef=setdiff(union(obj.CylInplanePosition,[10000 10000 ],'rows'),[10000 10000 ],'rows');
                AffiVar.InplaneXYRef=sortrows(AffiVar.InplaneXYRef,2);   %unnecs
                Pcoor=[-1 -1];
                for i=1:length(obj.AGroup); Pcoor=union(Pcoor, obj.CylInplanePosition(obj.AGroup(i),:),'rows');end
                Pcoor=setdiff(Pcoor,[-1 -1],'rows');
                Ncoor=[-1 -1];
                for i=1:length(obj.BGroup);Ncoor=union(Ncoor, obj.CylInplanePosition(obj.BGroup(i),:),'rows');end
                Ncoor=setdiff(Ncoor,[-1 -1],'rows');
                RelateTable=zeros(2,length(obj.Zbase1));
                        AffiVar.DeterCylDirect=~xor(side~=1,ismember(StartP,obj.AGroup));      
                        switch AffiVar.DeterCylDirect
                            case 1    %PList =even|move up     ,   NList=odd|move down   
                        RRTABLE=[ [Pcoor ,linspace(0,2*size(Pcoor,1)-2,size(Pcoor,1))'] ; [Ncoor ,linspace(1,2*size(Ncoor,1)-1,size(Ncoor,1))']];        
                        RelateTable(1,:)=[obj.AGroup obj.BGroup];
                        RT21=zeros(size(obj.AGroup));
                        RT22=zeros(size(obj.BGroup));
                        for k=1:length(RT21) ;
                            [~,~,CylinderIndex] =intersect( obj.CylInplanePosition(RelateTable(1,k),:) ,RRTABLE(:,1:2),'rows');
                            RT21(k)=RRTABLE(CylinderIndex,3)  ;
                        end
                        for k=1:length(RT22) ; 
                            [~,~,CylinderIndex] =intersect( obj.CylInplanePosition(  RelateTable(1,length(RT21)+k),:) ,RRTABLE(:,1:2),'rows');
                            RT22(k)=RRTABLE(CylinderIndex,3)  ;
                        end
                        RelateTable(2,:)=[RT21 RT22];
                        RelateTable=sortrows(RelateTable',1)';                
                         AffiVar.ShiftCol=10;  
                         AffiVar.Shiftrow=10;

                            case 0    %NList =even|move up      PList=odd|move down   
                        RRTABLE=[ [Pcoor ,linspace(1,2*size(Pcoor,1)-1,size(Pcoor,1))'] ; [Ncoor ,linspace(0,2*size(Ncoor,1)-2,size(Ncoor,1))']];
                         RelateTable(1,:)=[obj.BGroup obj.AGroup ]  ;
                        RT21=zeros(size(obj.BGroup));
                        RT22=zeros(size(obj.AGroup));
                        for k=1:length(RT21) ; 
                             [~,~,CylinderIndex] =intersect( obj.CylInplanePosition(RelateTable(1,k),:) ,RRTABLE(:,1:2),'rows');
                            RT21(k)=RRTABLE(CylinderIndex,3)  ;
                        end
                        for k=1:length(RT22) ;
                            [~,~,CylinderIndex] =intersect( obj.CylInplanePosition(RelateTable(1,length(RT21)+k),:) ,RRTABLE(:,1:2),'rows');
                            RT22(k)=RRTABLE(CylinderIndex,3)  ;       
                        end       
                           RelateTable(2,:)=[RT21 RT22];      
                           RelateTable=sortrows(RelateTable',1)';
                         AffiVar.ShiftCol=9;  
                         AffiVar.Shiftrow=10;  
                        end
                obj.XYCoorOrder=setdiff(union(obj.CylInplanePosition,[10000 10000 ],'rows'),[10000 10000 ],'rows');
                for j=1:size(RelateTable,2)
                   CheckXY=obj.CylInplanePosition( RelateTable(1,j),:);
                   [~,~,CylIndex]= intersect( CheckXY ,obj.XYCoorOrder,'rows');
                   RelateTable(3,j)=CylIndex;
                end
                obj.RTable=  RelateTable;        % 1st rows -> index in PArtBundle ; 2th rows-> cadnano index(start from 0, even means go up)  ;   3rd-> Matlab use(start from 1) 
                [~,index]=sortrows(obj.RTable',3);
                Vec=obj.RTable(2,:);
                TTRelateVec=Vec(index);
                obj.RelateVec=union(TTRelateVec,[],'stable');
        end  %%end of FindRTable and RVec.
        
        function obj=FindScafRouting(obj,stem,side,SavePair,CycleConnectTo,DeterCylDirect )  %scaf
            path=stem.getPath;
            CZdata=zeros(2*length(path),2);
            %side=0  -> start from bottom,      side=1  -> start from top
            ExtPath=[path(1) path path(end)];
            for i=2:length(ExtPath)-1   %jigzaw
                Movedir=mod(i+side,2);
                switch Movedir
                    case 0  %move up
                        CZdata(2*i-3,1)=ExtPath(i);
                        RefH=min([obj.Zbase1(ExtPath(i-1))  obj.Zbase1(ExtPath(i))]);
%                         CZdata(2*i-3,2)=FindZ2V2( ExtPath(i-1),ExtPath(i),obj.CylInplanePosition,RefH,0,obj.AGroup,DeterCylDirect,0 );
                         CZdata(2*i-3,2)=RefH-obj.ssScafOver;
                        CZdata(2*i-2,1)=ExtPath(i);
                        RefH=max([obj.Zbase2(ExtPath(i))  obj.Zbase2(ExtPath(i+1))]);
                        CZdata(2*i-2,2)=FindZ2V2( ExtPath(i),ExtPath(i+1),obj.CylInplanePosition,RefH,1,obj.AGroup,DeterCylDirect,1 ); 
                        CZdata(2*i-2,2)=RefH+obj.ssScafOver;
                    case 1
                        CZdata(2*i-3,1)=ExtPath(i);
                        RefH=max([obj.Zbase2(ExtPath(i-1))  obj.Zbase2(ExtPath(i))]);
%                         CZdata(2*i-3,2)=FindZ2V2( ExtPath(i-1),ExtPath(i),obj.CylInplanePosition,RefH,1,obj.AGroup,DeterCylDirect,1 );
                        CZdata(2*i-3,2)=RefH+obj.ssScafOver;
                        CZdata(2*i-2,1)=ExtPath(i);
                        RefH=min([obj.Zbase1(ExtPath(i))  obj.Zbase1(ExtPath(i+1))]);
                        CZdata(2*i-2,2)=FindZ2V2( ExtPath(i),ExtPath(i+1),obj.CylInplanePosition,RefH,0,obj.AGroup,DeterCylDirect,0 );
                        CZdata(2*i-2,2)=RefH-obj.ssScafOver;
                end  
            end
            cyclefinish=zeros(size(SavePair,1),1);

            while 1
               UsedCyl=  union( CZdata(:,1),[]);
               CC= ismember(CycleConnectTo,UsedCyl);
               CC=and(CC,~cyclefinish);
                ind=find(CC==1);
                for i=1:length(ind)
                    mergeto=CycleConnectTo(ind(i));       
                    cycle=SavePair(ind(i),:);
                   MH=sort([obj.Zbase1(mergeto) obj.Zbase2(mergeto) obj.Zbase1(cycle(1)) obj.Zbase2(cycle(1)) ]);
                    MH=floor(0.5*(MH(2)+MH(3)));
                    MH = FindZ2( cycle(1),mergeto, obj.CylInplanePosition,MH,0,obj.AGroup,DeterCylDirect );        
                    BreakIN=find(CZdata(:,1)==mergeto);
                    if length(BreakIN)>2   %find belonging segment
                        for j=1:length(BreakIN)/2
                           BINs=BreakIN(2*j-1:2*j);
                           CBin=CZdata(BINs,2);
                           if find(sort([CBin;MH])==MH)==2
                               BreakIN=BINs;
                               break
                           end
                        end
                    end
                    Ref=min([obj.Zbase1(cycle(1)) obj.Zbase1(cycle(2))]);
%                     CylcelZ1= FindZ2V2( cycle(1),cycle(2),obj.CylInplanePosition,Ref,0,obj.AGroup,DeterCylDirect,0 );
                    CylcelZ1=Ref-obj.ssScafOver;
                    Ref=max([obj.Zbase2(cycle(1)) obj.Zbase2(cycle(2))]);
%                     CylcelZ2= FindZ2V2( cycle(1),cycle(2),obj.CylInplanePosition,Ref,1,obj.AGroup,DeterCylDirect,1 );
                    CylcelZ2=Ref+obj.ssScafOver;
                    if CZdata(BreakIN(1),2)>CZdata(BreakIN(2),2)
                    CZdata=[CZdata(1:BreakIN(1),:)  ;mergeto MH ;cycle(1) MH ; cycle(1) CylcelZ2; cycle(2) CylcelZ2;cycle(2) CylcelZ1;cycle(1) CylcelZ1;cycle(1) MH-1;mergeto MH-1 ;  CZdata(BreakIN(2):end,:)];  % ff
                    else
                    CZdata=[CZdata(1:BreakIN(1),:)  ;mergeto MH-1 ;cycle(1) MH-1 ; cycle(1) CylcelZ1; cycle(2) CylcelZ1;cycle(2) CylcelZ2;cycle(1) CylcelZ2;cycle(1) MH;mergeto MH ;  CZdata(BreakIN(2):end,:)] ;           
                    end
                 cyclefinish(ind(i))=1;
                end

                if nnz(cyclefinish)==length(cyclefinish)
                    break;
                end    
            end

            XYZdata=zeros(size(CZdata,1),3);
            for i=1:size(CZdata,1)
                XYZdata(i,1:2)=obj.CylInplanePosition(CZdata(i,1),:);
                XYZdata(i,3)=CZdata(i,2);
            end
            obj.ScafRoutingXYZ=XYZdata;
            obj.ScafRoutingCZ=CZdata;
        end  %end of FindScafRouting
        
        function [obj,MM,MCylinderIndex]=ConvertScaf(obj)
            RRTable=obj.RTable
            SeqI=obj.ScafRoutingCZ;
            SeqI(:,2)=SeqI(:,2);    %+32*ones(size(SeqI(:,2)));  %in case of reach bottom bound
            Seq=[SeqI(1:end-1,1:2) SeqI(2:end,1:2)];
           Seq(:,2)= round(Seq(:,2));Seq(:,4)= round(Seq(:,4));
            
            QQ=obj.CylInplanePosition;
            MCylinderIndex=size(union(QQ,[10000 10000 ],'rows'),1)-1;
            DCellList=cell(MCylinderIndex,1);
            MM=max(max(Seq)); 
            MM=32*ceil(MM/32);
            for i=1:MCylinderIndex   %create shell to store information
              DCellList{i}=zeros(MM,4);
            end
            C0=-1;
            Z0=-1;
            for i1=1:2:size(Seq,1)     % i1=1 3 5 7 ......
               CylinderIndex=obj.RTable(3,Seq(i1,1));   
                if Seq(i1,2)>Seq(i1,4)  %move down
                       for j1=Seq(i1,2):-1:Seq(i1,4)
                           if j1==0  
                               break
                           end
                           if (i1==size(Seq,1))  && (j1==Seq(i1,4))
                            DCellList{CylinderIndex}(j1,:)=[C0,Z0,-1,-1];    
                              break 
                           end
                            if    j1==Seq(i1,4)                    
                            C1=   obj.RTable(3,Seq(i1+1,3));    %  C1=  Seq(i1+1,3);
                            ZNext=  j1;
                            DCellList{CylinderIndex}(j1,:)=[C0,Z0,C1,ZNext]    ;
                            C0=CylinderIndex;
                            Z0=j1;             
                            else   
                            C1= obj.RTable(3,(Seq(i1,1)));
                            ZNext=j1-1;
                            DCellList{CylinderIndex}(j1,:)=[C0,Z0,C1,ZNext];
                             C0=CylinderIndex;
                             Z0=j1;
                            end
                       end   
                else     %move up
                      for j1=Seq(i1,2):1:Seq(i1,4)
                         if (i1==size(Seq,1))  && (j1==Seq(i1,4))
                         DCellList{CylinderIndex}(j1,:)=[C0,Z0,-1,-1];    
                              break 
                         end
                        if    j1==Seq(i1,4)
                        C1=  obj.RTable(3,Seq(i1+1,3));
                        ZNext=  j1;
                        DCellList{CylinderIndex}(j1,:)=[C0,Z0,C1,ZNext] ; 
                        C0=CylinderIndex;
                        Z0=j1;
                        else   
                        C1= obj.RTable(3,Seq(i1,1));
                        ZNext=j1+1;
                        DCellList{CylinderIndex}(j1,:)=[C0,Z0,C1,ZNext];
                       C0=CylinderIndex;
                       Z0=j1;
                        end
                      end
                end
               if j1==0
                  break
               end
            end
            for x=1:MCylinderIndex
                for y=1:length(DCellList{x})
                   if sum(DCellList{x}(y,:))==0 
                    DCellList{x}(y,:)=[-1 -1 -1 -1];
                   end
                end  
            end     
            obj.DigitScaf=DCellList;
        end  %end of converscaf
        
        function obj=FindStap(obj,type,side)   % Get properties: StapList , stapBP ,stapBPinCell
                 InitialStappEnds=zeros(2*length(obj.Zbase1),2);
                 for i=1:length(obj.Zbase1)
                 InitialStappEnds(2*i-1:2*i,1)= obj.RTable(3,  i);
                 InitialStappEnds(2*i-1,2)=obj.Zbase1(i);
                 InitialStappEnds(2*i,2)=obj.Zbase2(i);
                 end 
                 n=size(InitialStappEnds,1)/2;
                 StapleList=cell(n,1);
                 for u=1:size(InitialStappEnds,1)/2   
                    StapleList{u}=InitialStappEnds(2*u-1:2*u,:) ; 
                 end             
                obj.StapList=StapleList;  %cylinder index use 3rd;
                
            if type==1  %no xover between cyl (initial stap)
                return;
            else
                DCellList=obj.DigitScaf;
                G2=[15 16]  ; %West
                G4=[31 32]   ;%East
                G3=[7 8]   ;%South
                G1=[23 24]   ; %North

                if side==1
                 G4=[15 16] ;  
                 G2=[31 32] ;
                 G1=[7 8];
                 G3=[23 24];  
                end
                period=32;
                AdjM=obj.CylAdjMat;    %use CylAdjM in PArtBundle
                AdjM2=AdjM;
                AdjM2= AdjM2.*triu(ones(size(AdjM2)),1);
                
                
                clearance=8;
                [U,V]=find(AdjM~=0);   %means cylinder in 1st
                
                U2=obj.RTable(3,U)';
                V2=obj.RTable(3,V)';
                Edgeof3Index=setdiff(union([U2 V2],[10000 10000 ],'rows'),[10000 10000 ],'rows');                
                k=1;
                EdgeList=cell(size(Edgeof3Index,1)/2,2);

                for i=1:size(Edgeof3Index,1)
                 if ismember(find(Edgeof3Index(i,1)==obj.RTable(3,:))   ,obj.AGroup) 
                    EdgeList{k,1}=[ Edgeof3Index(i,1) , Edgeof3Index(i,2)];   %in 3rd index;
                    k=k+1;
                 end
                end
                TwoDExpre=[];
                for edge=1:size(EdgeList,1)
                    PosCyl=EdgeList{edge,1}(1);
                    NegCyl=EdgeList{edge,1}(2);
                    EndsOfPosCyl=InitialStappEnds(find(InitialStappEnds(:,1)==PosCyl),2);
                    EndsOfNegCyl=InitialStappEnds(find(InitialStappEnds(:,1)==NegCyl),2);
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


                    NWSE=FindNWSE(PosCyl,NegCyl,obj.XYCoorOrder);
                    switch NWSE
                        case 1
                            for zz=G1(1):period:maxAll
                                zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                                zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                                if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4)  && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied 
                                       EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                                       TwoDExpre(end+1,:)=[PosCyl zz];
                                       TwoDExpre(end+1,:)=[NegCyl zz];
                                       TwoDExpre(end+1,:)=[PosCyl zz+1];                      
                                       TwoDExpre(end+1,:)=[NegCyl zz+1];
                                end
                            end
                         case 2
                            for zz=G2(1):period:maxAll
                                 zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                                zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                                if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4) && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied 
                                       EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                                       TwoDExpre(end+1,:)=[PosCyl zz];
                                       TwoDExpre(end+1,:)=[NegCyl zz];
                                       TwoDExpre(end+1,:)=[PosCyl zz+1];                      
                                       TwoDExpre(end+1,:)=[NegCyl zz+1];
                                end
                            end               
                         case 3
                            for zz=G3(1):period:maxAll
                                zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                                zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                                if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4) && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied 
                                       EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                                       TwoDExpre(end+1,:)=[PosCyl zz];
                                       TwoDExpre(end+1,:)=[NegCyl zz];
                                       TwoDExpre(end+1,:)=[PosCyl zz+1];                      
                                       TwoDExpre(end+1,:)=[NegCyl zz+1];
                                end
                            end
                          case 4
                            for zz=G4(1):period:maxAll
                                zzk=zz+0.1;[~,Izzk]=sort([zzk;EndsOfNegCyl2]);Izzk=find(Izzk==1);
                                zzl=zz+0.1;[~,Izzl]=sort([zzl;EndsOfPosCyl2]);Izzl=find(Izzl==1);
                                if (sum(DCellList{PosCyl}(zz,:))~=-4) &&(sum(DCellList{NegCyl}(zz,:))~=-4) && mod(Izzk(1),2)==0  && mod(Izzl(1),2)==0 %Both two sides are occupied 
                                       EdgeList{edge,2}=[EdgeList{edge,2} zz zz+1];
                                       TwoDExpre(end+1,:)=[PosCyl zz];
                                       TwoDExpre(end+1,:)=[NegCyl zz];
                                       TwoDExpre(end+1,:)=[PosCyl zz+1];                      
                                       TwoDExpre(end+1,:)=[NegCyl zz+1];
                                end
                            end
                    end
                end
                obj.stapBP=TwoDExpre;
                obj.OstapBP=TwoDExpre;
                obj.stapBPinCell=EdgeList;
            end
        end %end of FindStap
        
        function obj=GetScafXOver(obj)    %Get property: ScafXover
            Seq=[obj.ScafRoutingCZ(1:end-1,:) obj.ScafRoutingCZ(2:end,:)];  
            for i=1:size(Seq,1)
                Seq(i,1)= obj.RTable(3, Seq(i,1));
                Seq(i,3)= obj.RTable(3, Seq(i,3));
            end
            
            
            BPL=obj.stapBPinCell;
            MCyl=max(max(Seq(:,[1 3])));
            Corner=cell(MCyl,1);
            MidCorner=cell(MCyl,1);
            SeqP=[Seq(:,[1 2]) ;Seq(end,3) Seq(end,4)];
            for i=1:size(Seq,1)
               cly= Seq(i,1);
               Position=Seq(i,2);   
               Corner{cly}=union(Corner{cly}, Position);     
               if i==size(Seq,1)
                cly= Seq(i,3);
               Position=Seq(i,4);   
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
            OSBP=obj.stapBP;
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
        

        
        function obj=FindStapStep2(obj)  %---------
            StapleCell=obj.StapList; obj.stapBP
            AppOrder=randperm(size(obj.stapBP,1)/4);
            for i=1:size(obj.stapBP,1)/4    
                InputM= obj.stapBP(4*AppOrder(i)-3:4*AppOrder(i),:);   
                [ ~,StapleCell] = AddXoverInStapple(StapleCell,InputM );     
            end
            tt=StapleCell;
            [ StapleCell ] = CalibStapDir( StapleCell,obj.RelateVec);
            [ NewStapList ] = OrganizeStapList( StapleCell,obj ,tt);
%             [ NewStapList] = CheckCycleStapAndBreak( NewStapList,obj.stapBP );
            obj.StapList2=NewStapList;
%             obj.StapList2=StapleCell;

        end %end of FindStapStep2
        
        
        function obj=ConvertStap(obj,MM,NoOfCyls)
%              RelateTable=obj.RTable;
             StappCellList=cell(NoOfCyls,1);
             for i=1:NoOfCyls   %create shell to store information
              StappCellList{i}=-1*ones(MM,4);
             end
             QQ=obj.XYCoorOrder;
             ZZZ=setdiff(union(QQ,[10000 10000 ],'rows'),[10000 10000 ],'rows');
             heads=[];
             for j=1:size(obj.StapList2,1)    
                 C0=-1  ;Z0=-1;
                 StappStrand2=  obj.StapList2{j};
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
                     [~,~,CylinderIndex] =intersect( QQ(  StappStrand(k,1),:) ,ZZZ,'rows');
                    if   StappStrand(k,1)~=CylinderIndex
                        sdfs=123;
                    end
                        
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
                                    for z=StappStrand(k,2):1:StappStrand(k+1,2);
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
            obj.DigitStap=StappCellList;
            obj.HeadOfStap=heads;
        end %end of ConvertStap
        
        function obj=ExportJSON(obj,AffiVar)
            dat=loadjson('CadNanoScaffOnly.json');
            NNdat=dat;
            jj=randi([10,99],1,1);
            JJ=strcat('Test',int2str(jj));
            JJ=strcat(JJ,'.json');
            file_name=JJ;
            NNdat.name=file_name;
            ColorChoice=[29184      243362     1507550     3355443     5749504     7536862     8947848    11184640    12060012    13369344    16225054];
%             CCSTR={'1.206e+07' '5.7495e+06' '5.7495e+06' '1.1185e+07' '8.9478e+06' '1.1185e+07' '7.5369e+06'};
            for cylindex=2:size(obj.DigitScaf,1)
               NNdat.vstrands{cylindex} =NNdat.vstrands{1};       %initialized
            end
            DCellList2=obj.DigitScaf;
            StappCellList2=obj.DigitStap;
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
            for cylindex=1:size(obj.DigitScaf,1)
                Colorstap= find(obj.HeadOfStap(:,1)==cylindex);
                ColorM=zeros(length(Colorstap),2);
                ColorM(:,1)=obj.HeadOfStap(Colorstap,2) + 0*ones(length(Colorstap),1)   ;
                ColorM(:,2)= ColorChoice(randi([1 length(ColorChoice)],size(Colorstap)))   ;  
                ColorM2=sortrows(ColorM,1);
                NNdat.vstrands{cylindex}.stap_colors ={ColorM2};
            %        NNdat.vstrands{cylindex}.stap_colors =cell(0,1);
               NNdat.vstrands{cylindex}.num =obj.RelateVec(cylindex);
               NNdat.vstrands{cylindex}.scafLoop =cell(0,1);
                 NNdat.vstrands{cylindex}.stap=StappCellList2{cylindex};
%                 NNdat.vstrands{cylindex}.stap=-1*ones(size(DCellList2{cylindex}));
               NNdat.vstrands{cylindex}.skip=zeros(1,size(obj.DigitScaf{1},1));
               NNdat.vstrands{cylindex}.scaf=DCellList2{cylindex};     %   
                NNdat.vstrands{cylindex}.stapLoop=cell(0,1);  
                NNdat.vstrands{cylindex}.loop=zeros(1,size(obj.DigitScaf{1},1));
                NNdat.vstrands{cylindex}.col=(obj.XYCoorOrder(cylindex,1)-obj.XYCoorOrder(1,1) )/2+AffiVar.ShiftCol;
                 
                NNdat.vstrands{cylindex}.row=(obj.XYCoorOrder(cylindex,2)-obj.XYCoorOrder(1,2) )/2+AffiVar.Shiftrow;
                 
                 [NNdat.vstrands{cylindex}.col  NNdat.vstrands{cylindex}.row]
                 sdf=1323424;
            end
            TTtext=savejson('Title',NNdat,'ArrayIndent',0,'Compact',1 );
            TTtext(1:10)=[];
            TTtext(end-1:end)=[];
            TTtext(TTtext==' ')=[];
            TTtext(TTtext=='	')=[];
            IOfSC=strfind(TTtext, '"stap_colors"');
            for i=1:length(IOfSC)
                P=IOfSC(i);
                if  ~strcmp(TTtext(P+15),'[')      
                   FirstBracket=P+14;
                   Segment=TTtext(P:P+30);
                   NextBracket=strfind(Segment,']')-1+P;
                   BeSubed=TTtext(FirstBracket:NextBracket);
                   Subing=strcat('[',BeSubed,']');
                   TTtext=strrep(TTtext, BeSubed, Subing);
                   IOfSC=strfind(TTtext, '"stap_colors"'); 
                end
            end
            fileID = fopen(file_name,'w');
            fprintf(fileID,TTtext);
            fclose(fileID);
            ReadTest=loadjson(file_name)      
            obj.DigitScafCAD=DCellList2;
            obj.DigitStapCAD=StappCellList2;
        end %end of ExportJSON
        
        function T=rotx(obj,rad)
            T=[1 0 0;0 cos(rad) -sin(rad);0 sin(rad) cos(rad)];
        end
        function T=roty(obj,rad)
            T=[cos(rad) 0 sin(rad);0 1 0; -sin(rad) 0 cos(rad)];
        end
        function T=rotz(obj,rad)
            T=[cos(rad) -sin(rad) 0;sin(rad) cos(rad)  0 ;0  0 1];
        end
        
    end
    
end


