classdef StressIntensityMultiAvTopologyOptimization < StressIntensityMultiTopologyOptimization
    % Multi-load-case stress-intensity ESO, "sum" aggregation.
    %
    % Two constructor forms are accepted -- see D3 in docs/integration-rules.md:
    %
    %   (Rmin, FEAnalyses, maxais, penal, volFr, is_const)
    %       CAS_Arm's form. Weights default to uniform.
    %
    %   (Rmin, FEAnalyses, alphas, maxais, penal, volFr, is_const)
    %       Vibrations' form. alphas holds one weight per load case and is what
    %       the harmonic-load trade-off sweep varies; dropping it would leave
    %       that sweep running with every weight equal and every result
    %       quietly meaningless.
    methods
        function obj = StressIntensityMultiAvTopologyOptimization(Rmin, FEAnalyses, varargin)
            [alphas, maxais, penal, volFr, is_const] = ...
                StressIntensityMultiTopologyOptimization.parseWeightedArgs( ...
                    numel(FEAnalyses), varargin{:});
            obj = obj@StressIntensityMultiTopologyOptimization( ...
                Rmin, FEAnalyses, alphas, maxais, penal, volFr, is_const, "sum");
        end
    end
end
