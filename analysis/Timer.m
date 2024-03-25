classdef Timer
    %TIMER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        prompt,N,Ni,startTime;
    end
    
    methods
        function obj = Timer(prompt,N,Ni)
            obj.prompt=prompt;
            obj.N=N;
            obj.Ni=Ni;
            obj.startTime=toc;
        end
        
        function timeMonitoring(obj,k)
            if mod(k,obj.Ni)==0
                tm=toc-obj.startTime;
                stm=datestr(datenum(0, 0, 0, 0, 0, tm), 'DD:HH:MM:SS');
                stmr=datestr(datenum(0, 0, 0, 0, 0, (obj.N-k)/k*tm), 'DD:HH:MM:SS');
                fprintf("%s: %6d, elapsed time :%s, remaining time: %s \n",obj.prompt,k,stm,stmr);
            end
        end
    end
end

