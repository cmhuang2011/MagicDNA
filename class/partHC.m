classdef partHC < part
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type='HC' ;
    end
    
    methods
        function obj=partHC(allGraphics)            
        obj = obj@part( allGraphics ) ;
                
        
        
        end
        
        function clearAxe(obj,src,evn)         
%             LinesRemain=findobj(gca,'tag','Remain') ;
%             LinesRemain(:).XData
            currentfield = strcat(obj.type,'_','Cross');
            axes(obj.allGraphics.allaxes.(currentfield)); cla; hold on ;
            xlim([0 50]);ylim([0 50]);           
            Cylradius=1 ;  th=linspace(0,360,30);          
           [ HClattice ] = findHClattice( Cylradius ,[50 50]) ;
           for i=1:HClattice.nColumn
               for j=1:HClattice.nrow
                xx=HClattice.HCcenter( (i-1)*HClattice.nrow+j   ,1)+Cylradius*cosd(th);
                yy=HClattice.HCcenter( (i-1)*HClattice.nrow+j   ,2)+Cylradius*sind(th);
                plot(xx,yy,'Color',[0.8 0.8 0.8],'HitTest','off');
               end
           end     
%          
%            ;  % back to initial
            ListProperties=properties(obj) ;
            Exclude={'type','allGraphics','Cylradius','PartSec','OVList'};
            for k=1:length(ListProperties)
                if ~ismember( ListProperties{k} , Exclude)
                obj.( ListProperties{k})=[];
                end
            end            
%             sdfsf=3
            for k=1:length(obj.OVList)
               plot(obj.OVList{k}(:,1),obj.OVList{k}(:,2),'--k'); 
            end
          

            cla(obj.allGraphics.allaxes.HC_Extrude);


         
        end % end of clearAxe
        
        
        
        
        
        
    end
    
end

