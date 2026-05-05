function configs = armLoadConfigs(kind)
% armLoadConfigs  Standard full-arm load-configuration sets.
%
%   "six" returns the symmetric min/max bending, torsion, and shear cases.
%   "sixPlusTension" adds the straight-arm axial/tension case.
    if nargin < 1
        kind = "six";
    end
    kind = string(validatestring(char(kind), {'six', 'sixPlusTension'}));

    configs = {
        struct('name', 'min_bending', 'label', 'Min M_z', 'betas', -[0 0 0 180 180 180 180]);
        struct('name', 'min_torsion', 'label', 'Min M_s', 'betas', -[0 45 45 45 270 180 180]);
        struct('name', 'min_shear',   'label', 'Min T_y', 'betas', -[0 0 180 0 180 180 180]);
        struct('name', 'max_bending', 'label', 'Max M_z', 'betas', [0 0 0 180 180 180 180]);
        struct('name', 'max_torsion', 'label', 'Max M_s', 'betas', [0 45 45 45 270 180 180]);
        struct('name', 'max_shear',   'label', 'Max T_y', 'betas', [0 0 180 0 180 180 180]);
    };

    if kind == "sixPlusTension"
        configs{end + 1, 1} = struct('name', 'max_tension', 'label', 'Max N', ...
            'betas', [0 0 0 0 0 0 0]);
    end
end
