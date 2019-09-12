classdef BundleCylinderSQ < BundleCylinder
    %description:
    %    
    %   Detailed explanation goes here
    
    properties
         type='SQ' ;
          AGroupGoUp=1;   %default value
    end
    properties (Constant)
        Lattice='Square';
        period=[32,3];
        template=[2 ,3, 4, 5,7 8,10 ,11, 12,13,15 16,18 ,19,20 21,23, 24,26,27,28,29,31 ,32];        
        CorrectComst= [5.5 ,0]-0.6;   % check with Cadnano 12/6/2016
    end
    
    methods
        function obj=BundleCylinderSQ(type,varargin)
            obj = obj@BundleCylinder(type, varargin{:} ) ;               
        end
       
        function XY=findExtraCylInplanePosition(obj, RTable, qColRow,fRefCyl)
%             qColRow
            qColRow=unique(qColRow,'rows','legacy') ;
            ranRefCyl =fRefCyl  ;
            RefXY =  obj.CylInplanePosition(RTable(ranRefCyl,2),:)  ;

            RefColRow =  RTable(ranRefCyl, 6:7)  ;
            SQMapping= findSQlatticeMapping( 1 ,[50 50]) ; 
            [~,ind] =ismember(RefXY ,  SQMapping(:,1:2) ,'rows') ;

%             if ind==0
%                 sdfsdg=3
%             end
            
            RefCyldColRow_shift =  RefColRow- SQMapping(ind,3:4) ; % 

            [~,ind2] =ismember( qColRow - RefCyldColRow_shift ,  SQMapping(:,3:4) ,'rows')   ;
            XY =  SQMapping(ind2,1:2)  ; %  
            
            
            
        end
        
    end
    
end

