classdef UseCadnanoToDraw <jsonObject
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    properties (Dependent)
        Corner_k ;
    end
    methods
        function  ExampleUseInCommandLine(obj)
            JSON_routing = UseCadnanoToDraw ;
            XYZ_mapping = [zeros(2,1) , [cosd(0:90:270) ;sind(0:90:270)] ]' ;
            XYZ = JSON_routing.AssignXY(XYZ_mapping ) ;
            JSON_routing.ExportXYZ(XYZ) ;
        end
        
        
        function Out = get.Corner_k(obj)
            Out = obj.ScafCornerRep{1}  ;% assuming only one scaffold
            for k =1 : size(Out,1)
                Out(k,1) = find(obj.NumList(:,2) ==    Out(k,1) ) ;
            end
        end
        
        function XYZ = AssignXY(obj , XY) 
           Kroute =  obj.Corner_k ;
           XYZ  = zeros(size(Kroute,1) ,3) ;
           XYZ(:,3) = Kroute(:,2) ;
           XYZ(:,1:2) = XY( Kroute(:,1) ,: ) ;         
           
           MaxMinusMixInXYZDirs = max(XYZ)- min(XYZ)
           fprintf('Check if XYZ is desired. Modify if need. Then, export it by  JSON_routing.ExportXYZ(XYZ) ; \n')
        end
        
        function ExportXYZ(obj,XYZ) 
            
            WhereToSave = uigetdir ;
            
            prompt = {'Enter File name:'};  dlg_title = 'Input';
            num_lines = 1;   defaultans = {'XYZ_1'};
            file3_name = inputdlg(prompt,dlg_title,num_lines,defaultans);
            fileID2 = fopen(strcat([WhereToSave filesep file3_name{1}],'.xyz')  ,'w') ;
            fprintf(fileID2,'HeadNode 1\n'   );  % header

            
            for k = 1:size(XYZ ,1)
                fprintf(fileID2,'%6.2f %6.2f %6.2f \n'  , XYZ(k,1),XYZ(k,2),XYZ(k,3) );  % header
%                 fprintf(fileID2,'%6.2f %6.2f %6.2f \n'  , XYZ(k,3),XYZ(k,2),XYZ(k,1) );  % header
            end
            fclose(fileID2);
            
        end
        
        
    end
    
end

