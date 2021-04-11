classdef BundleCylinderHC < BundleCylinder
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type='HC' ;
        AGroupGoUp=0;   %default value
        %-----------MagicDNA2
        
        
        
    end
    properties (Constant,Hidden)
        Lattice='HoneyComb';
        period=[21, 2];
        template=[1, 2, 4, 5, 8, 9, 11, 12, 15, 16, 18, 19];
        %         AGroupGoUp=0 ;   %default value
        CorrectComst=[0,5.25]+0.4 ; %               [0,5.25]+0.88 ; % 06042019         correct=[5.5,0]+2
    end
    
    methods
        function obj=BundleCylinderHC(type,varargin)
            obj = obj@BundleCylinder(type, varargin{:} ) ;
            
        end
        
        function XY=findExtraCylInplanePosition(obj, RTable, qColRow,fRefCyl)
            %             qColRow=unique(qColRow,'rows','legacy') ;
            
            ranRefCyl =fRefCyl  ;
            RefXY =  obj.CylInplanePosition(RTable(ranRefCyl,2),:)  ;
            
            RefColRow =  RTable(ranRefCyl, 6:7)  ;
            HCMapping= findHClatticeMapping( 1 ,[50 50]) ;
            %             [~,ind] =ismember(RefXY ,  HCMapping(:,1:2) ,'rows') ;
            
            [~,ind] =ismembertol(RefXY ,  HCMapping(:,1:2) ,0.01 ,'Byrows',true) ;   % bug, 07/01/2020 overhang
            %             LIA = ismembertol(A,B,tol)
            
            RefCyldColRow_shift =  RefColRow- HCMapping(ind,3:4) ; %
            
            [~,ind2] =ismember( qColRow - RefCyldColRow_shift ,  HCMapping(:,3:4) ,'rows')   ;
            XY =  HCMapping(ind2,1:2)  ; %
            
            
            %             dsf=3
        end
        
    end
    
end

