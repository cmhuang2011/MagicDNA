classdef part < handle
    %description:
    % this class mainly handles the graphic in main only , 
    % later convert data into class BundleCylinder as part designs
    %   Detailed explanation goes here
    
    properties
        allGraphics
        Zbase1, Zbase2 ,CylInplanePosition, hScatter,DetectArray,DetectArrayIn1,DetectArrayInFinal
        OVList , InnerLoop, InnerDetectResult,unit, Nxy,cylh,HaveCylZ1,HaveCylZ2,CrossSection,
        suggesth,PartSec,       
        OuterVertexs=[];
        OuterLoop=[];
        Cylradius=1;
        
        GUISavePart ;
    end
    
    methods
        function obj=part(allGraphics)
            obj.allGraphics=allGraphics;
            
%             return
            %-----setup graphic functions
             if isa(obj,'partHC')                
                 subcls= 'HC' ;
             elseif isa(obj,'partSQ')
                 subcls= 'SQ' ;             
             end
                        
             currentfield = strcat(subcls,'_','Clear');
             allGraphics.allbtn.(currentfield).Callback=@(src,evn)clearAxe(obj,src,evn)         ;    
             
             field2 = strcat(subcls,'_','Cross');
             allGraphics.allaxes.(field2).ButtonDownFcn=@(src,evn)axesCrossCallBack(obj,src,evn)         ;    %                  

             field3 = strcat(subcls,'_','MakeColseLoop');
             allGraphics.allbtn.(field3).Callback=@(src,evn)MakeCloseLoop(obj,src,evn)         ;   
             
             field = strcat(subcls,'_','Show3D');
             allGraphics.allbtn.(field).Callback=@(src,evn)Show3D(obj,src,evn)         ;        
             
             field = strcat(subcls,'_','Adjust');
             allGraphics.allbtn.(field).Callback=@(src,evn)AdjustHeight(obj,src,evn)         ;                   

             %-----------
             field = strcat(subcls,'_','SaveSection');
             allGraphics.allbtn.(field).Callback=@(src,evn)SaveSection(obj,src,evn)         ;                   
             
             fd = strcat(subcls,'_','SaveSection');
             allGraphics.allpop.(fd).Callback=@(src,evn)popCallback(obj,src,evn)        ;   
             %----------------------
             fd = strcat(subcls,'_','SavePart');
             allGraphics.allbtn.(fd).Callback=@(src,evn)SavePart(obj,src,evn)        ;   

             fd = strcat(subcls,'_','LoadPart');
             allGraphics.allbtn.(fd).Callback=@(src,evn)LoadPart(obj,src,evn)        ;   

             fd = strcat(subcls,'_','ClearPart');
             allGraphics.allbtn.(fd).Callback=@(src,evn)ClearPart(obj,src,evn)        ;   

             
        end
        
        function  SavePart(obj,src,evn)
            Psave=obj;
            uisave({'Psave'},'Partx');            
        end  % end of SavePart
        function  LoadPart(obj,src,evn)
            uiopen('load')
           if ~strcmp(class(obj),class(Psave))
            error('import object is not correct !!') ;    return ;
           end            
            ListProperties=properties(obj) ;
            Exclude={'type','allGraphics','Cylradius'};
            for k=1:length(ListProperties)   
                if ~ismember( ListProperties{k} , Exclude)
                obj.( ListProperties{k})=Psave.( ListProperties{k});    
                end
            end
            fd2 = strcat(obj.type,'_','CylModel');           
            axes(obj.allGraphics.allaxes.(fd2));   cla;
            hold on;view(3); xlabel('nm'); ylabel('nm'); zlabel('base');
            DrawPart( obj.GUISavePart,[1 0 1],1 );       
            fd = strcat(obj.type,'_','SaveSection');
             obj.allGraphics.allpop.(fd).Value=1;
             obj.allGraphics.allpop.(fd).String={'New'};
             iPartSec=length(obj.PartSec);
             for i=1:iPartSec
                obj.allGraphics.allpop.(fd).String{i}=strcat('Sec',num2str(i));
             end
             obj.allGraphics.allpop.(fd).String{end+1}=strcat('New');
            
        end  % end of SavePart        
        
        function ClearPart(obj,src,evn) 
            ListProperties=properties(obj) ;
            Exclude={'type','allGraphics','Cylradius'};
            for k=1:length(ListProperties)
                if ~ismember( ListProperties{k} , Exclude)
                obj.( ListProperties{k})=[];
                end
            end
            subcls=obj.type;
            fd = strcat(subcls,'_','Cross');           
            cla(obj.allGraphics.allaxes.(fd));
            fd = strcat(subcls,'_','Extrude');           
            cla(obj.allGraphics.allaxes.(fd));
            fd = strcat(subcls,'_','CylModel');           
            cla(obj.allGraphics.allaxes.(fd));
            fd = strcat(subcls,'_','Section');           
            cla(obj.allGraphics.allaxes.(fd));
        end
        
        function popCallback(obj,src,evn)
            fd = strcat(obj.type,'_','Section');
            axes(obj.allGraphics.allaxes.(fd)); cla; hold on;
            try
            DrawPart( obj.PartSec{src.Value},[1 0 1],1 );    
            catch
            end            
        end
        
        function AdjustHeight(obj,src,evn)

            
            field = strcat(obj.type,'_','Extrude');
            fd_nm = strcat(obj.type,'_','Unit_nm');
            fd_cylAxes = strcat(obj.type,'_','Extrude');
            ed3=  strcat(obj.type,'_','Adjustlower');  ed4=  strcat(obj.type,'_','Adjustupper');
            if obj.allGraphics.allrb.(fd_nm).Value==1  %convert to base
                coeff=0.34;
            else
                coeff=1;
            end
            IZ1_e=str2double(get(obj.allGraphics.alledit.(ed3) ,'String'));
            IZ2_e=str2double(get(obj.allGraphics.alledit.(ed4)  ,'String'));
            IZ1_e=round(IZ1_e/coeff);   %interger unit bp
            IZ2_e=round(IZ2_e/coeff);
%             axes(handles.axes2);
            axes(obj.allGraphics.allaxes.(fd_cylAxes)); cla;
            %-------------
 HaveCyl=find(obj.DetectArrayInFinal==1);
 iselect=find(obj.hScatter.SizeData==80);
 HaveCylZ1=obj.HaveCylZ1;
 HaveCylZ2=obj.HaveCylZ2;
 for i=1:length(HaveCyl)
     if ismember(HaveCyl(i),iselect)
   HaveCylZ1(i)=HaveCylZ1(i)-  IZ1_e   ;
   HaveCylZ2(i)=HaveCylZ2(i)+  IZ2_e  ;
     end
     RowPlusColumn= mod(HaveCyl(i),obj.Nxy(2))+floor(HaveCyl(i)/obj.Nxy(2));
     if mod(RowPlusColumn,2)==1
   cylh{i}= cylinder2P(obj.allGraphics.allaxes.(fd_cylAxes),1, 20,[obj.DetectArray(1,HaveCyl(i)) obj.DetectArray(2,HaveCyl(i)) 0.34*HaveCylZ1(i) ],[obj.DetectArray(1,HaveCyl(i)) obj.DetectArray(2,HaveCyl(i)) 0.34*HaveCylZ2(i)],[1 0 1]) ;      
     else
   cylh{i}= cylinder2P(obj.allGraphics.allaxes.(fd_cylAxes),1, 20,[obj.DetectArray(1,HaveCyl(i)) obj.DetectArray(2,HaveCyl(i)) 0.34*HaveCylZ1(i)],[obj.DetectArray(1,HaveCyl(i)) obj.DetectArray(2,HaveCyl(i)) 0.34*HaveCylZ2(i)],[0 1 0]) ;
     end
 end
 obj.cylh=cylh;
obj.HaveCylZ1=HaveCylZ1;obj.HaveCylZ2=HaveCylZ2;
%------------
%             sdfsf=3
            axis square
        end
        
        function obj=SaveSection(obj,src,evn)       
             fd = strcat(obj.type,'_','SaveSection');
             fd2 = strcat(obj.type,'_','CylModel');
             HaveCyl=find(obj.DetectArrayInFinal==1);
             
             if  isa(obj,'partSQ')       %if SQ
               Part1=BundleCylinderSQ(1,[],obj.HaveCylZ1',obj.HaveCylZ2',obj.DetectArray(:,HaveCyl)')   ;
            elseif  isa(obj,'partHC')         %if HC
               Part1=BundleCylinderHC(1,[],obj.HaveCylZ1',obj.HaveCylZ2',obj.DetectArray(:,HaveCyl)')   ;
        %        PartEx.HCXY=findHClattice( handles.Cylradius ,[50 50]) ;
             end
             obj.PartSec{obj.allGraphics.allpop.(fd).Value}=Part1;
             obj.OVList{obj.allGraphics.allpop.(fd).Value}=obj.OuterVertexs;

             if length(obj.PartSec)==1
                obj.GUISavePart=obj.PartSec{1};
            elseif length(obj.PartSec)>1
                TemPart=obj.PartSec{1};
                for i=1:length(obj.PartSec)-1
                     if  isa(obj,'partSQ')       %if SQ
                    TemPart=BundleCylinderSQ(2,TemPart,obj.PartSec{i+1});
                     elseif isa(obj,'partHC')     %if HC
                    TemPart=BundleCylinderHC(2,TemPart,obj.PartSec{i+1});
                     end

                end
                obj.GUISavePart=TemPart;
             end  
             
             axes(obj.allGraphics.allaxes.(fd2));cla;
             hold on;view(3); xlabel('nm'); ylabel('nm'); zlabel('nm');
             DrawPart( obj.GUISavePart,[1 0 1],1 );
             
             obj.allGraphics.allpop.(fd).Value=1;
             obj.allGraphics.allpop.(fd).String={'New'};
             iPartSec=length(obj.PartSec);
             for i=1:iPartSec
                obj.allGraphics.allpop.(fd).String{i}=strcat('Sec',num2str(i));
             end
             obj.allGraphics.allpop.(fd).String{end+1}=strcat('New');
%              handles.OVList{handles.popupmenu1.Value}=handles.OuterVertexs;
             axis square ;
             sdfsf=3;
             
        end
        
        function obj=Show3D(obj,src,evn)
            field = strcat(obj.type,'_','Extrude');
            fd_nm = strcat(obj.type,'_','Unit_nm');
            AddCyl=zeros(size(obj.DetectArrayInFinal));
            DeleteCyl=zeros(size(obj.DetectArrayInFinal));
            for i=1:length(DeleteCyl)
               if  sum(obj.hScatter.CData(i,:)==[0 0 1])==3 || sum(obj.hScatter.CData(i,:)==[0 0 0])==3  ; AddCyl(i)=1;   %black or blue
               elseif sum(obj.hScatter.CData(i,:)==[1 0 0])==3;DeleteCyl(i)=1;
               end
            end    

            HaveCyl=find(AddCyl==1);      NewDAIF=zeros(size(obj.DetectArrayInFinal));
            NewDAIF(HaveCyl)=1;    obj.DetectArrayInFinal=NewDAIF;
            ax1= obj.allGraphics.allaxes.(strcat(obj.type,'_','Cross')) ;
            axes(obj.allGraphics.allaxes.(field));  ax=gca; cla; hold on;
            xlabel('nm');ylabel('nm');zlabel('nm');
            
            fd1 = strcat(obj.type,'_','show3Dlower');
            fd2 = strcat(obj.type,'_','show3Dupper');
%             hAdd=handles.hAdd;
            IZ1=str2double(get(obj.allGraphics.alledit.(fd1) ,'String'));
            IZ2=str2double(get(obj.allGraphics.alledit.(fd2)  ,'String'));
            if obj.allGraphics.allrb.(fd_nm).Value==1  %convert to base
                coeff=0.34;
            else
                coeff=1;
            end
            IZ1=round(IZ1/coeff);   %interger unit bp
            IZ2=round(IZ2/coeff);
            obj.cylh=cell(length(HaveCyl),1);
             for i=1:length(HaveCyl)   %alwasy show in nm in axes2
               RowPlusColumn= mod(HaveCyl(i),obj.Nxy(2))+floor(HaveCyl(i)/obj.Nxy(2));
               if mod(RowPlusColumn,2)==1
               obj.cylh{i}= cylinder2P(ax,obj.Cylradius, 12,[obj.DetectArray(1,HaveCyl(i)) obj.DetectArray(2,HaveCyl(i)) 0.34*IZ1],[obj.DetectArray(1,HaveCyl(i)) obj.DetectArray(2,HaveCyl(i)) 0.34*IZ2],[0.8 0.5 0.5]) ;  
               else
               obj.cylh{i}= cylinder2P(ax,obj.Cylradius, 12,[obj.DetectArray(1,HaveCyl(i)) obj.DetectArray(2,HaveCyl(i)) 0.34*IZ1],[obj.DetectArray(1,HaveCyl(i)) obj.DetectArray(2,HaveCyl(i)) 0.34*IZ2],[0.8 0.5 0.5]) ;      
               end
             end
             axis auto ; light('Position',[-1 0 0],'Style','local') ;
            obj.HaveCylZ1=ones(length(HaveCyl),1)*IZ1;
            obj.HaveCylZ2=ones(length(HaveCyl),1)*IZ2;
            obj.CrossSection=flip(reshape(obj.DetectArrayInFinal,flip(obj.Nxy)),1);
            
            [u,v]=find(obj.CrossSection==1);w=u+v;
            NofGreen=sum(mod(w,2)==0);   NofPurple=sum(mod(w,2)==1);
            
            PosAddCyl=Neibor3(obj.CrossSection);
            BinaryAddMatrix=zeros(size(PosAddCyl));
            [u1,v1]=find(PosAddCyl==1);
            for i=1:length(u1)
                if mod(u1(i)+v1(i),2)==1;        BinaryAddMatrix(u1(i),v1(i))=2;
                elseif  mod(u1(i)+v1(i),2)==0;    BinaryAddMatrix(u1(i),v1(i))=1;
                end        
            end
            BinaryAddMatrix=flip(BinaryAddMatrix,1);
             
            N=24;
            th=linspace(0,360,N);
            R=repmat([ 2 1.2],[1 N/2]);
            if ~isempty(obj.suggesth)
            delete(obj.suggesth(:));
            end             
            
            if  NofGreen>NofPurple   %Reduce green or increase purple
                AddPurpleCyli=find(BinaryAddMatrix==2);
                 plotxydata=zeros(N+1,length(AddPurpleCyli),2);
                for i=1:N
                    plotxydata(i,:,1)=obj.hScatter.XData(AddPurpleCyli)+ones(1,length(AddPurpleCyli))*R(i)*cosd(th(i));
                    plotxydata(i,:,2)=obj.hScatter.YData(AddPurpleCyli)+ones(1,length(AddPurpleCyli))*R(i)*sind(th(i));
                end
                plotxydata(N+1,:,:)=plotxydata(1,:,:);
                obj.suggesth=plot(ax1, plotxydata(:,:,1), plotxydata(:,:,2),'-','color',[1 0 1]);
            elseif NofGreen<NofPurple
                AddGreenCyli=find(BinaryAddMatrix==1);
                plotxydata=zeros(N+1,length(AddGreenCyli),2);
                for i=1:N
                    plotxydata(i,:,1)=obj.hScatter.XData(AddGreenCyli)+ones(1,length(AddGreenCyli))*R(i)*cosd(th(i));
                    plotxydata(i,:,2)=obj.hScatter.YData(AddGreenCyli)+ones(1,length(AddGreenCyli))*R(i)*sind(th(i));
                end
                plotxydata(N+1,:,:)=plotxydata(1,:,:);
                obj.suggesth=plot(ax1, plotxydata(:,:,1), plotxydata(:,:,2),'-','color',[0 1 0]);
            end            
            view(3)       ;   axis square;     
        end  % end of Show3D
        
        function obj=MakeCloseLoop(obj,src,evn )
%              field2 = strcat(obj.type,'_','Loop_Outer') ;
             field = strcat(obj.type,'_','Cross');
             field2 = strcat(obj.type,'_','Loop_Outer') ;
             field3 = strcat(obj.type,'_','Loop_Inner') ;

             axes(obj.allGraphics.allaxes.(field));
             
             if get(obj.allGraphics.allrb.(field2) ,'Value')==1  && isempty(obj.OuterLoop)
                new2Points=obj.OuterVertexs(1,:); oldvertexs=obj.OuterVertexs;
                NOuterVertexs=CheckIntersect(new2Points, oldvertexs,2);
                if sum(obj.OuterVertexs(1,:)==NOuterVertexs(end,:))==2
                    OuterVertexs=NOuterVertexs;
                     plot(OuterVertexs(end-1:end,1),OuterVertexs(end-1:end,2),'-*b','HitTest','off','tag','Remain' );
                     text(OuterVertexs(end-1,1),OuterVertexs(end-1,2),num2str(size(OuterVertexs,1)-1));
                if  strcmp(obj.type,'SQ')  % handles.SelectedTab2==1   %if SQ
                   [TestPointInG,in,~,Nxy]= DetectCylArray2( obj.OuterVertexs) ;
                elseif  strcmp(obj.type,'HC')     %if HC
                   [TestPointInG,in,~,Nxy]= DetectCylArray2HC( obj.OuterVertexs,obj.Cylradius) ;
                end

                   obj.Nxy=Nxy;
                   ScatterColor=zeros(length(in),3);ScatterColor(find(in==1),3)=1;ScatterColor(find(in==0),1)=1;
                   obj.hScatter=scatter(TestPointInG(1,:)',TestPointInG(2,:)',20*ones(length(TestPointInG(2,:)),1),ScatterColor);      
                end
                obj.OuterLoop=OuterVertexs;  
                obj.DetectArray=TestPointInG;
                obj.DetectArrayIn1=in;
                obj.DetectArrayInFinal= in; 
             end
            
             if get(obj.allGraphics.allrb.(field3) ,'Value')==1   && ~isempty(obj.InnerLoop{end})
              InnerVertexs=obj.InnerLoop{end};
              new2Points=InnerVertexs(1,:); oldvertexs=InnerVertexs;
              NInnerVertexs=CheckIntersect(new2Points, oldvertexs,2);
              if sum(InnerVertexs(1,:)==NInnerVertexs(end,:))==2
                      InnerVertexs=NInnerVertexs;
                      plot(InnerVertexs(end-1:end,1),InnerVertexs(end-1:end,2),'-.b','HitTest','off');
                      [in2,~] = inpolygon(obj.DetectArray(1,:),obj.DetectArray(2,:),InnerVertexs(:,1)',InnerVertexs(:,2)');
                      obj.InnerDetectResult(length(obj.InnerLoop),:)=in2;
                      InOut=zeros(size(obj.DetectArray,2),1);
                      for i=1:length(obj.InnerLoop)
                         InOut=or(InOut, obj.InnerDetectResult(i,:)');
                      end
                      InOut(obj.DetectArrayIn1==0)=0;         
            %           hScatter=cell(1,size(handles.DetectArray,2));  
                      delete(findall(gcf,'Type','Scatter'));
                      FinalInorOut=xor(InOut,obj.DetectArrayIn1');
                      ScatterColor=zeros(length(FinalInorOut),3);ScatterColor(find(FinalInorOut==1),3)=1;ScatterColor(find(FinalInorOut==0),1)=1;
                      obj.hScatter=scatter(obj.DetectArray(1,:)',obj.DetectArray(2,:)',20*ones(length(obj.DetectArray(1,:)),1),ScatterColor);
                      obj.DetectArrayInFinal= xor(InOut(:),obj.DetectArrayIn1(:));     
                      obj.InnerLoop{end+1}=[];
              end
             end            
             obj.hScatter.ButtonDownFcn=@(src,evnt)selectpoint(obj,src,evnt);

        end   % end of  MakeCloseLoop       
        
        function obj=selectpoint(obj,src,evnt)
        field = strcat(obj.type,'_','Cyl_Select');
        field2 = strcat(obj.type,'_','Cyl_AddOrDel') ;
        if get(obj.allGraphics.allrb.(field) ,'Value')==1
        ax = gca; SS=ax.CurrentPoint;SS=[SS(1,1);SS(1,2)];
        XYdata=[obj.hScatter.XData;obj.hScatter.YData];
        dis=zeros(size(XYdata,2),1);   %distances to every point-> find the closest one
            for i=1:size(XYdata,2)
            dis(i)=norm(XYdata(:,i)-SS);
            end
        idis=find(dis==min(dis));
            if sum(obj.hScatter.CData(idis,:)==[0 0 1])==3 && evnt.Button==1
            obj.hScatter.CData(idis,:)=[0 0 0];
            obj.hScatter.SizeData(idis)=80;   
            end
            if sum(obj.hScatter.CData(idis,:)==[0 0 0])==3 && evnt.Button==3
            obj.hScatter.CData(idis,:)=[0 0 1];
            obj.hScatter.SizeData(idis)=20;   
            end
        end

        if get(obj.allGraphics.allrb.(field2) ,'Value')==1
        ax = gca;
        SS=ax.CurrentPoint;SS=[SS(1,1);SS(1,2)];
        XYdata=[obj.hScatter.XData;obj.hScatter.YData];
        dis=zeros(size(XYdata,2),1);   %distances to every point-> find the closest one
            for i=1:size(XYdata,2)
            dis(i)=norm(XYdata(:,i)-SS);
            end
        idis=find(dis==min(dis));
            if sum(obj.hScatter.CData(idis,:)==[0 0 1])==3 && evnt.Button==3
           obj.hScatter.CData(idis,:)=[1 0 0];
            elseif sum(obj.hScatter.CData(idis,:)==[1 0 0])==3 && evnt.Button==1
            obj.hScatter.CData(idis,:)=[0 0 1];
            end
        end
            
        end % end of slectPoint
        
        function obj=axesCrossCallBack(obj,src,evn)
            currentfield = strcat(obj.type,'_','Cyl_Draw') ;
            field2 = strcat(obj.type,'_','Loop_Outer') ;
            field3 = strcat(obj.type,'_','Loop_Inner') ;
            ax=src;
%             ax = gca;
            if get(obj.allGraphics.allrb.(currentfield) ,'Value')==1 
                if isempty(obj.OuterLoop) && get(obj.allGraphics.allrb.(field2) ,'Value')==1
                    currentPoint=ax.CurrentPoint;
                    NewPointList=CheckIntersect(currentPoint(1,1:2), obj.OuterVertexs,1);
                    obj.OuterVertexs=NewPointList;
                    if size(obj.OuterVertexs,1) > 1
                        for i=1:size(obj.OuterVertexs,1)-1
                            plot(obj.OuterVertexs(i:i+1,1),obj.OuterVertexs(i:i+1,2),'-*b','HitTest','off','tag','Remain');
            %                 text(OuterVertexs(i,1),OuterVertexs(i,2),num2str(i));
                        end
                    end
                    NewPointListInn=[];
                elseif get(obj.allGraphics.allrb.(field3) ,'Value')==1
                    n=length(obj.InnerLoop);
                    currentPoint=ax.CurrentPoint;
                    if n==0
                       InnerVertexs=[];
                       NewPointListInn=CheckIntersect(currentPoint(1,1:2), InnerVertexs,1);
                        if size(NewPointListInn,1) > 1
                            for i=1:size(NewPointListInn,1)-1
                                plot(NewPointListInn(i:i+1,1),NewPointListInn(i:i+1,2),'-.b');
                            end
                        end
                        obj.InnerLoop{1}=NewPointListInn;
                    else
                        InnerVertexs=obj.InnerLoop{n};
                        NewPointListInn=CheckIntersect(currentPoint(1,1:2), InnerVertexs,1);
                        if size(NewPointListInn,1) > 1
                            for i=1:size(NewPointListInn,1)-1
                                plot(NewPointListInn(i:i+1,1),NewPointListInn(i:i+1,2),'-.b');
            %                     text(NewPointListInn(i,1),NewPointListInn(i,2),num2str(i));
                            end
                        end
                        obj.InnerLoop{n}=NewPointListInn;
                    end
                end
            end
           
        end  % end of axesCrossCallBack
        
        
        
        
        
        
        
    end
    
end

