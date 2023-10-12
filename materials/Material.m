classdef Material < handle
    properties
        name;
    end
    
    methods
        function obj = Material(name)
            obj.name = name;
        end
    end
    
end

