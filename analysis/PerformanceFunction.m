classdef (Abstract) PerformanceFunction < handle
   
    methods(Abstract)
      computeValue(obj,x);
    end

end

