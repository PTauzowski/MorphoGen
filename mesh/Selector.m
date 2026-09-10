classdef Selector
        
    properties
        fn;
        tolerance;
    end
    
    methods
        function obj = Selector(fn,varargin)
            obj.fn=fn;
            obj.tolerance=1.0E-04;
            if nargin==2
                obj.tolerance = varargin{1};
            end
        end
        
        function s = select( obj, points )
            % Two calling conventions are in use and both are supported.
            %
            %   signed distance : Selector(@(x) x(:,1) - x0)
            %                     selected where |fn(x)| < tolerance.
            %                     This is the convention the selectX/selectY/
            %                     selectZ helpers produce and the one to
            %                     prefer in new code.
            %
            %   predicate       : Selector(@(x) abs(x(:,1)-x0) < 0.001)
            %                     fn already returns logical; used as-is.
            %                     158 call sites across 53 files still use
            %                     this form, so it is accepted rather than
            %                     migrated. Transitional -- see D1; call
            %                     sites should move to selectX/Y/Z later.
            %
            % Treating a predicate as a distance selects the COMPLEMENT of
            % the intended set (abs(true)=1 is never < tolerance, abs(false)=0
            % always is), which is silent and produces plausible garbage.
            if islogical(obj.fn)
                s = obj.fn;                       % precomputed mask
                return;
            end

            v = obj.fn(points);
            if islogical(v)
                s = v;                            % predicate convention
            else
                s = abs(v) < obj.tolerance;       % signed-distance convention
            end
        end
        
    end
end

