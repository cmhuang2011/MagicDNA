classdef partSQ < part
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
          type='SQ' ;
    end
    
    methods
        function obj=partSQ(allGraphics)
            
        obj = obj@part( allGraphics ) ;
        
        end
        
        
        
        function clearAxe(obj,src,evn)   
%             LinesRemain=findobj(gca,'tag','Remain')
            currentfield = strcat(obj.type,'_','Cross');
            axes(obj.allGraphics.allaxes.(currentfield)); cla; hold on ;
            ax=gca;
            xlim([0 50]);ylim([0 50]);
            
            Cylradius=1 ;  th=linspace(0,360,30);          
            for i=0:2*Cylradius:50
                for j=0:2*Cylradius:50
                    x=i+Cylradius*cosd(th);
                    y=j+Cylradius*sind(th);
                    plot(x,y,'Color',[0.8 0.8 0.8],'HitTest','off');
                end
            end    
%         
%             back to initial
            %------------
            ListProperties=properties(obj) ;
            Exclude={'type','allGraphics','Cylradius','PartSec','OVList'};
            for k=1:length(ListProperties)
                if ~ismember( ListProperties{k} , Exclude)
                obj.( ListProperties{k})=[];
                end
            end
             for k=1:length(obj.OVList)
               plot(obj.OVList{k}(:,1),obj.OVList{k}(:,2),'--k'); 
            end
            
            %----------
%             sdfsfd=3
            cla(obj.allGraphics.allaxes.SQ_Extrude);
            
        end  % end of clearAxe
        
        
        
        
        
    end
    
end

