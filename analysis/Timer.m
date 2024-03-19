classdef Timer
    %TIMER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        N,Ni,startTime;
    end
    
    methods
        function obj = Timer(N,Ni)
            obj.N=N;
            obj.Ni=Ni;
            obj.startTime=toc;
        end
        
        function timeMonitoring(obj,k)
            if mod(k,obj.Ni)==0
                tm=toc-obj.startTime;
                stm=datestr(datenum(0, 0, 0, 0, 0, tm), 'HH:MM:SS');
                stmr=datestr(datenum(0, 0, 0, 0, 0, (obj.N-k)/k*tm), 'HH:MM:SS');
                fprintf("Iteration: %6d, elapsed time :%s, remaining time: %s \n",k,stm,stmr);
            end
        end
    end
end

