function s = curveParamsToDisplayStruct(params)
% curveParamsToDisplayStruct  Return a struct with the 14 optimisation
%   parameters in canonical order (suppresses legacy fields like angleDeg).
    names = curveParamNames();
    s = struct();
    for i = 1:numel(names)
        s.(names{i}) = params.(names{i});
    end
end
