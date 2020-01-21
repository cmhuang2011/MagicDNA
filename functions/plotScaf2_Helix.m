function varargout=plotScaf2_Helix( GetHyperB,CornerNotation  )

% old :function pScaf2H=plotScaf2_Helix( GetHyperB,CornerNotation  )
%Use in: cadnano and oxDNA panels
%   Detailed explanation goes here

% SacfR=GetHyperB.ScafRouting ;
% MaxBase=max(SacfR(:,3));

SacfR=GetHyperB.ScafRouting{1} ;
for k = 2 : length(GetHyperB.ScafRouting)
 SacfR= [SacfR ; GetHyperB.ScafRouting{k}    ] ;
end

MaxBase=max(SacfR(:,3));



skipPattern1=9:60:MaxBase;   % skip mod2 =0   %test 4/20
skipPattern2=39:60:MaxBase;


     SaveGHelixStap=cell(1,length(GetHyperB.containBundle));
    for Bundlei=1:length(GetHyperB.containBundle)
        QQWithPM10Bases =GetHyperB.containBundle{Bundlei}.HelixXYZG;       %scaf domain
%                 QQWithPM10Bases =GetHyperB.containBundle{Bundlei}.HelixXYZGStap;             
        SaveGHelixStap{Bundlei}= QQWithPM10Bases;
    end

    StapHelixCell=cell(1,length(GetHyperB.StapList3) );
    StapHelix2Cell=cell(1,length(GetHyperB.StapList3) );
    AppedixStapHelRCell=cell(1,length(GetHyperB.StapList3) );
    StapSeqATCGCell=cell(1,length(GetHyperB.StapList3) );
    AppedixBVecCell=cell(1,length(GetHyperB.StapList3) );
    AppedixNVecCell=cell(1,length(GetHyperB.StapList3) );
    skipinStap=0;

    pScaf2H=cell(length(GetHyperB.StapList3),1) ;
      for stai=1:length(CornerNotation)   %actually mean scaffold, in case multi scaffolds in futures
      StapAll=CornerNotation{stai}    ;
%      StapAll=GetHyperB.StapList3{stai}    ;
      
      
     StapHelix=zeros(10000,3); kc=1;   % pre-allocate
     StapHelix2=zeros(10000,3);
     AppedixStapHelR=zeros(10000,3);  % [Bundle , Cyl, Base];
        for edgeinSCR=1:2:size(StapAll,1)
             C5Cyl=StapAll(edgeinSCR,1);
             bundle=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,1)   ;%multi-section needs be cautious

             Bundle=GetHyperB.containBundle{bundle};
             Cyl=GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,2);   %C2 express
             
             if  Cyl~=-1   % overhang cylinders, Not used anymore
                 BaseStart=StapAll(edgeinSCR,2);   BaseEnd=StapAll(edgeinSCR+1,2);

                    if   length(Cyl)~=1   % not constant crosssection
                        bundle=unique(bundle);
                        for whicCyl=1:length(Cyl)
                        lB1=GetHyperB.containBundle{bundle}.Zbase1(Cyl(whicCyl));
                        lB2=GetHyperB.containBundle{bundle}.Zbase2(Cyl(whicCyl)) ;   
                        [lB1,lB2];
                        if BaseStart>=lB1 && BaseStart<=lB2
                            Cyl=Cyl(whicCyl);
                            break;
                        end
                        end
                    end
                 ExtendBases= 11;  %11  
                 RelativeBS=BaseStart-GetHyperB.containBundle{bundle}.Zbase1(Cyl)+ExtendBases ;
                 RelativeBE=BaseEnd-GetHyperB.containBundle{bundle}.Zbase1(Cyl)+ExtendBases;
                 QQ=linspace(RelativeBS,RelativeBE,abs(RelativeBE-RelativeBS)+1) ;
                 if  strcmp(Bundle.Lattice, 'Square') 
                     if BaseStart<BaseEnd  %go to left------------------stap
                         skipinStap=skipinStap+length(intersect(QQ, skipPattern2-GetHyperB.containBundle{bundle}.Zbase1(Cyl)+11,'stable'));
                        QQ=setdiff(QQ, skipPattern2-GetHyperB.containBundle{bundle}.Zbase1(Cyl)+11,'stable');
                        skipP=skipPattern2;
                     else
                          skipinStap=skipinStap+length(intersect(QQ, skipPattern1-GetHyperB.containBundle{bundle}.Zbase1(Cyl)+11,'stable')); 
                         QQ=setdiff(QQ, skipPattern1-GetHyperB.containBundle{bundle}.Zbase1(Cyl)+11,'stable');  
                           skipP=skipPattern1;
                     end
                 else
                     skipP=[];        %only square lattice needs skip
                 end     

                 size(SaveGHelixStap{bundle}{Cyl})   ;
                 PartHelix2=SaveGHelixStap{bundle}{Cyl}(QQ,:);
                 if ~isempty(GetHyperB.containBundle{bundle}.TransformMatrix2)  %if had done simulation ----------------
                     TAll=GetHyperB.containBundle{bundle}.TransformMatrix2;
            %          QQQ=TAll(1:3,1:3)'*PartHelix2'-TAll(1:3,4)*ones(1, size(PartHelix2,1)  );%---old
                     QQQ=TAll(1:3,1:3)*PartHelix2'+TAll(1:3,4)*ones(1, size(PartHelix2,1)  );
                     StapHelix2(kc:kc+size(PartHelix2,1)-1,: )=QQQ' ;
                     AppedixStapHelR(kc:kc+size(PartHelix2,1)-1,: )= [ones(size(PartHelix2,1),1)*bundle, ones(size(PartHelix2,1),1)*Cyl , setdiff(linspace(BaseStart,BaseEnd,1+abs(BaseStart- BaseEnd))', skipP ,'stable')  ]   ;
                 end

             else  % Use overhangs
              
                   ColRow =  GetHyperB.RelateTable(GetHyperB.RelateTable(:,5)==C5Cyl,6:7);
                   if mod(GetHyperB.RelateTable(C5Cyl,4) ,2) ==0
                       SimilarDirCylders = and( mod(GetHyperB.RelateTable(:,4) ,2)==0 , GetHyperB.RelateTable(:,1)==bundle) ;
                   else
                       SimilarDirCylders = and( mod(GetHyperB.RelateTable(:,4) ,2)==1 , GetHyperB.RelateTable(:,1)==bundle) ;
                   end
                   fSimilarDirCylders = find(SimilarDirCylders) ; fSimilarDirCylders=fSimilarDirCylders(1) ; 
                   ThefSimilarDirCylder= GetHyperB.RelateTable(fSimilarDirCylders,2) ;

                   RefCyl = and(GetHyperB.RelateTable(:,1)==bundle, GetHyperB.RelateTable(:,2)~=-1  ) ;
                   fRefCyl= find(RefCyl) ;  fRefCyl=fRefCyl(1) ;
                   InplaneXY=Bundle.findExtraCylInplanePosition(  GetHyperB.RelateTable, ColRow,fRefCyl)   ;%
                   BaseStart=StapAll(edgeinSCR,2);   BaseEnd=StapAll(edgeinSCR+1,2);
                   BasesArr =  linspace(BaseStart,BaseEnd , abs(BaseStart-BaseEnd) +1  ) ;
                   Global_XYZ=Bundle.HelixRegardlessCylinder(1,0,InplaneXY,BasesArr,ThefSimilarDirCylder) ;   %  already assign as scaffold domain
                   PartHelix2= Global_XYZ(:,1:3) ;
                 if ~isempty(GetHyperB.containBundle{bundle}.TransformMatrix2)  %if had done simulation ----------------
                     TAll=GetHyperB.containBundle{bundle}.TransformMatrix2;
            %          QQQ=TAll(1:3,1:3)'*PartHelix2'-TAll(1:3,4)*ones(1, size(PartHelix2,1)  );%---old
                     QQQ=TAll(1:3,1:3)*PartHelix2'+TAll(1:3,4)*ones(1, size(PartHelix2,1)  );
                     StapHelix2(kc:kc+size(PartHelix2,1)-1,: )=QQQ' ;
                     if Cyl==-1  % staple overhang
                         
                     AppedixStapHelR(kc:kc+size(PartHelix2,1)-1,: )= [ones(size(PartHelix2,1),1)*bundle, ones(size(PartHelix2,1),1)*Cyl , setdiff(linspace(BaseStart,BaseEnd,1+abs(BaseStart- BaseEnd))', [] ,'stable')  ]   ;
                     else
                     AppedixStapHelR(kc:kc+size(PartHelix2,1)-1,: )= [ones(size(PartHelix2,1),1)*bundle, ones(size(PartHelix2,1),1)*Cyl , setdiff(linspace(BaseStart,BaseEnd,1+abs(BaseStart- BaseEnd))', skipP ,'stable')  ]   ;
               
                     end      
                 end
                   
                   
             end
             
             
             
             
             PartHelix2V2=PartHelix2;
%              PartHelix2V2(:,1:2)=1.3*PartHelix2V2(:,1:2);       %---------debug
             StapHelix(kc:kc+size(PartHelix2,1)-1,: )=PartHelix2V2 ;     
             kc=kc+size(PartHelix2,1);    % accumulate
        end
        
        
        StapHelix( sum(StapHelix,2)==0,:)=[];  %delete extra
        RemovedInd=sum(StapHelix2,2)==0;
        StapHelix2( RemovedInd,:)=[];
        AppedixStapHelR(RemovedInd,:)=[];         
        StapHelixCell{stai}=StapHelix;
        StapHelix2Cell{stai}=StapHelix2;
        AppedixStapHelRCell{stai}=AppedixStapHelR;

%         [Incase,SeqInd]=ismember(AppedixStapHelR,AppedixAllHelR,'rows');
%             if nnz(Incase)~=length(Incase)   %has single-stranded staple
%                 sdfsf=234;
%             end
%         StapSeqATCGCell{stai}=  seqcomplement(ScafSeqATCG(SeqInd));   % complementary, A<->T, C<->G.------------
%         AppedixBVecCell{stai}=-AppedixBVec(SeqInd,:);   %find scaffold BVec, add negative
%         AppedixNVecCell{stai}=-AppedixNVec(SeqInd,:);   %find scaffold NVec, add negative
%         [stai, length(StapSeqATCGCell{stai})];
        pScaf2H{stai}=plot3(  StapHelix(:,1),StapHelix(:,2),StapHelix(:,3),'-k','LineWidth',1.5);
%         text( StapHelix(1,1), StapHelix(1,2), StapHelix(1,3),num2str(stai) );
      end


varargout{1} = pScaf2H ;


end

