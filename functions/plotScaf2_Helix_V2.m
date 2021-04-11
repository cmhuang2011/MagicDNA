function varargout=plotScaf2_Helix_V2( GetHyperB,CornerNotation,Isstap ,color , TM )

% old :function pScaf2H=plotScaf2_Helix( GetHyperB,CornerNotation  )
%Use in: cadnano and oxDNA panels
% Update for MagicDNA2, July 30 2020   
%
%


AllSkipBase =GetHyperB.skipBase ;
AllinsertBase =GetHyperB.insertBase ;

% Coeff= 0.85*0.4 ;

StapHelixCell=cell(1,length(CornerNotation) );
StapHelix2Cell=cell(1,length(CornerNotation) );
BasecCenterCell =  cell(1,length(CornerNotation) );

NVecCell =  cell(length(CornerNotation),1 );
BundleRoutCell=cell(length(CornerNotation),1 );
pScaf2H=cell(length(CornerNotation),1) ;
pScaf_center=cell(length(CornerNotation) ,1) ;
for stai=1:length(CornerNotation)   %actually mean scaffold, in case multi scaffolds in futures

    
    StapAll=CornerNotation{stai}  + [zeros(size(CornerNotation{stai},1) ,1 ) ,zeros(size(CornerNotation{stai},1),1 )]   ;
    %      StapAll=GetHyperB.StapList3{stai}    ;
    StapHelix=zeros(10000,3); kc=1;   % pre-allocate
    StapHelix2=zeros(10000,3);
    BaseCenterHelix =zeros(10000,3);    % oxdna center location
    NVec =zeros(10000,3);    % oxdna center location
    BundleRout = zeros(10000,2);    % BelongM
    for edgeinSCR=1:2:size(StapAll,1)
        C5Cyl=StapAll(edgeinSCR,1);
        bundle=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,1)   ;%multi-section needs be cautious
        bundle=unique(bundle);
        Bundle=GetHyperB.containBundle{bundle};
        Cyl=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,2) ;
        Cyl=Cyl(1);  %C2 express
        ColRow =  GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,6:7);
        IndInRT = find(GetHyperB.RelateTable(:,5)==C5Cyl) ; IndInRT=IndInRT(1); %   fixed Apr 25 2019,
        if mod(GetHyperB.RelateTable(IndInRT,4) ,2) ==0 %-------
            %                    bundle
            %                    GetHyperB.RelateTable(:,1)==bundle
            %                    mod(GetHyperB.RelateTable(:,4) ,2)==0
            SimilarDirCylders = and( mod(GetHyperB.RelateTable(:,4) ,2)==0 , GetHyperB.RelateTable(:,1)==bundle) ;
        else
            SimilarDirCylders = and( mod(GetHyperB.RelateTable(:,4) ,2)==1 , GetHyperB.RelateTable(:,1)==bundle) ;
        end
        fSimilarDirCylders = find(SimilarDirCylders) ; fSimilarDirCylders=fSimilarDirCylders(1) ;
        ThefSimilarDirCylder= GetHyperB.RelateTable(fSimilarDirCylders,2) ;
        RefCyl = and(GetHyperB.RelateTable(:,1)==bundle, GetHyperB.RelateTable(:,2)~=-1  ) ;
        fRefCyl= find(RefCyl) ;
        %                if bundle==3
        %                    sdfsf=3
        %                end
        
        fRefCyl=fRefCyl(1) ;
        % fRefCyl;
        InplaneXY=Bundle.findExtraCylInplanePosition(  GetHyperB.RelateTable, ColRow,fRefCyl)   ;%
        BaseStart=StapAll(edgeinSCR,2);   BaseEnd=StapAll(edgeinSCR+1,2);
        
%         if  strcmp(Bundle.Lattice, 'Square') && Cyl~=-1   % not on overhang
%             if xor(BaseStart<BaseEnd,Isstap~=1)  %go to left------------------stap, important, for scaffold it is opposite
%                 %                         skipP=skipPattern2;
%                 skipP=AllSkipBase(AllSkipBase(:,1)==C5Cyl,2) ;
%             else
%                 %                         skipP=skipPattern1;
%                 skipP=AllSkipBase(AllSkipBase(:,1)==C5Cyl,2) ;
%             end
%         else
%             skipP=[];        %only square lattice needs skip
%             
%         end

        skipP=AllSkipBase(AllSkipBase(:,1)==C5Cyl,2) ;
        if BaseStart<BaseEnd
            Vertical_NVec = [0 ,0 ,1] ;
        else
            Vertical_NVec = [0 ,0 ,-1] ;
        end
        
        
        
        BasesArr =  setdiff(linspace(BaseStart,BaseEnd , abs(BaseStart-BaseEnd) +1  ) , skipP,'stable') ;
                       %-----Introduce insert
               insertP=AllinsertBase(AllinsertBase(:,1)==C5Cyl,2) ;insertP=reshape(insertP, 1,[]);
               if ~isempty(intersect(insertP ,BasesArr ))
                   if BasesArr(1)>BasesArr(end)
                       BasesArr=sort([BasesArr,intersect(insertP ,BasesArr )] ,'descend')  ;
                   else
                       BasesArr=sort([BasesArr,intersect(insertP ,BasesArr )] ,'ascend')    ;
                   end
               end
               %----------
        
        Global_XYZ=Bundle.HelixRegardlessCylinder(1.06,Isstap,InplaneXY,BasesArr,ThefSimilarDirCylder) ;   %  already assign as scaffold domain (rr,isstap,......)
        % get 3d coordinate section by section
        PartHelix2= Global_XYZ(:,1:3) ;
        
        Global_XYZ_Center=Bundle.HelixRegardlessCylinder(0.59,Isstap,InplaneXY,BasesArr,ThefSimilarDirCylder) ;   %  already assign as scaffold domain (0.66(oxdna center),isstap,......)
        % has transformed to aseembly Tmatrix
        Global_XYZ_Center= Global_XYZ_Center(:,1:3) ;
        
        Vertical_NVec= ones( size(Global_XYZ_Center,1),1)*Vertical_NVec ;
        BundleR =ones( size(Global_XYZ_Center,1),1)*bundle;
        if ~isempty(GetHyperB.containBundle{bundle}.TransformMatrix2)  %if had done simulation ----------------
            if TM==2
                T_simulation=GetHyperB.containBundle{bundle}.SimulateTransMFromTM2;  % simulation matrix
            else
                %                  T_simulation=GetHyperB.containBundle{bundle}.TransformMatrix2;  %
                T_simulation=eye(4) ;
            end
            
            %                  T_assembly=GetHyperB.containBundle{bundle}.TransformMatrix2;  % transformation matrix  % already applied in "HelixRegardlessCylinder"
            %                  T_simulation=T_assembly;
            %                     QQQ=     transpose(  inv( T_assembly(1:3,1:3))*( PartHelix2' -T_assembly(1:3,4)*ones(1, size(PartHelix2,1)) )   );   % parallel position
            
            
            QQQ2= T_simulation(1:3,1:3)*PartHelix2'+ T_simulation(1:3,4)*ones(1, size(PartHelix2,1)  );
            %                  StapHelix2(kc:kc+size(PartHelix2,1)-1,: )=QQQ' ;
            StapHelix2(kc:kc+size(PartHelix2,1)-1,: )=transpose(QQQ2 ) ;
            
            
            WWW = T_simulation(1:3,1:3)*Global_XYZ_Center'+T_simulation(1:3,4)*ones(1, size(Global_XYZ_Center,1)  );
            BaseCenterHelix(kc:kc+size(PartHelix2,1)-1,: )=WWW' ;
            
            %                              XYZ=   transpose(obj.TransformMatrix2*LocalXYZ' );
            NVec(kc:kc+size(PartHelix2,1)-1,: )=transpose(Bundle.TransformMatrix2(1:3,1:3)*Vertical_NVec' )  ;
            BundleRout(kc:kc+size(PartHelix2,1)-1,: )=[BundleR ,BasesArr'] ;
            
%             if ~isempty(intersect(insertP ,BasesArr ))
%                 sdfsf=3
%             end
            
        end
        %          PartHelix2V2=PartHelix2;
        StapHelix(kc:kc+size(PartHelix2,1)-1,: )=PartHelix2 ;
        kc=kc+size(PartHelix2,1);    % accumulate
%         size(PartHelix2,1)
    end
    %
    
    StapHelix( sum(StapHelix,2)==0,:)=[];  %delete extra
    RemovedInd=sum(StapHelix2,2)==0;    % index ro remove over-allocate
    StapHelix2( RemovedInd,:)=[];
    BaseCenterHelix( RemovedInd,:)=[];
    NVec( RemovedInd,:)=[];
    BundleRout( RemovedInd,:)=[];
    
    
    
    StapHelixCell{stai}=StapHelix;
    StapHelix2Cell{stai}=StapHelix2;
    BasecCenterCell{stai}=BaseCenterHelix;
    NVecCell{stai}=NVec;
    BundleRoutCell{stai}=BundleRout;
    
    
    %     pScaf2H{stai}=plot3(  StapHelix2(:,1),StapHelix2(:,2),StapHelix2(:,3),'-k','Color',color);
    pScaf2H{stai}=plot3(  StapHelix2(:,1),StapHelix2(:,2),StapHelix2(:,3),'-k','LineWidth',1.5,'Color',color);
    
    
    pScaf_center{stai}=line(  BaseCenterHelix(:,1),BaseCenterHelix(:,2),BaseCenterHelix(:,3),'LineStyle','none', 'Marker','.');  % increase efficiency 08/26/2019
    %      pScaf_center{stai}=scatter3(  BaseCenterHelix(:,1),BaseCenterHelix(:,2),BaseCenterHelix(:,3),[], color,'.');
    
    %     pScaf2H{stai}=plot3(  StapHelix(:,1),StapHelix(:,2),StapHelix(:,3),'-k','LineWidth',1.5);
    %         text( StapHelix(1,1), StapHelix(1,2), StapHelix(1,3),num2str(stai) );
    %     dfsdg=4
end


varargout{1} = pScaf2H ;   % graphic object
varargout{2} = StapHelix2Cell ;  % 3D helix XYZ

varargout{3} = pScaf_center ;  % graphic object scatter
varargout{4} = BasecCenterCell ;  % 3D helix XYZ of oxdna center
varargout{5} = NVecCell ;  % N Vector of oxdna center

varargout{6} = BundleRoutCell ;  % N Vector of oxdna center


end

