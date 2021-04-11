classdef AssemblyCommand
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type=[];
        Trans_Rotate = []
        Global_Local = []
        Increment = []
        Direction = []
        Bundle=[];
    end
    
    methods
        
        function obj = AssemblyCommand(type, varargin )
            
            %TR,GL,Inc,Dir
            if strcmp(type,'XYZ_input')
                obj.type=type;
                obj.Trans_Rotate=varargin{1};
                obj.Global_Local=varargin{2};
                obj.Increment=varargin{3};
                
                switch varargin{4}
                    case 'q'
                        obj.Direction='a';
                    case 'a'
                        obj.Direction='q';
                    case 'w'
                        obj.Direction='s';
                    case 's'
                        obj.Direction='w';
                    case 'e'
                        obj.Direction='d';
                    case 'd'
                        obj.Direction='e';
                end
                obj.Bundle=varargin{5};
                
            elseif strcmp(type,'snapXYZ')
                obj.type=type;
                obj.Trans_Rotate=varargin{1};       % save rotation matrix which was calculated from orthoganal R
                obj.Bundle=varargin{2};
            end
            
        end
    end
    
end

