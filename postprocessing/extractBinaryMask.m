function solid = extractBinaryMask(x, method, opts)
% EXTRACTBINARYMASK  Extract a binary topology from a continuous density field.
%
%   solid = extractBinaryMask(x, method, opts)
%
%   Inputs:
%     x       [nElems x 1]  density field in [0, 1]
%     method  string        one of 'fixed', 'volume', 'hysteresis'
%     opts    struct        method-specific parameters (see below)
%
%   Methods and required opts fields:
%
%     'fixed'
%       opts.threshold   scalar threshold (default 0.5)
%
%     'volume'
%       opts.targetVF    target volume fraction; finds the largest t such that
%                        mean(x >= t) >= targetVF
%
%     'hysteresis'
%       opts.mesh        struct with .elems [nElems x nNodesPerElem]
%       opts.rhoHigh     unconditional seed threshold (default 0.55)
%       opts.rhoLow      minimum density to include if connected to a seed
%                        default: (0.05)^(1/penal)  — 5% penalised stiffness
%       opts.penal       SIMP exponent used only to set default rhoLow (default 3)
%
%       Flood-fills from elements with x >= rhoHigh, including adjacent
%       elements with x >= rhoLow.  Captures thin load paths that fall
%       below a simple 0.5 cut without including structurally dead material.
%
%   Output:
%     solid  [nElems x 1] logical

    x = x(:);
    if nargin < 3, opts = struct(); end

    switch lower(method)

        case 'fixed'
            threshold = optField(opts, 'threshold', 0.5);
            solid = x >= threshold;

        case 'volume'
            assert(isfield(opts, 'targetVF'), ...
                'extractBinaryMask: opts.targetVF is required for ''volume'' method.');
            t     = findVolumeThreshold(x, opts.targetVF);
            solid = x >= t;

        case 'hysteresis'
            assert(isfield(opts, 'mesh'), ...
                'extractBinaryMask: opts.mesh is required for ''hysteresis'' method.');
            penal   = optField(opts, 'penal',   3);
            rhoHigh = optField(opts, 'rhoHigh', 0.55);
            rhoLow  = optField(opts, 'rhoLow',  (0.05)^(1/penal));
            solid   = hysteresisFloodFill(x, opts.mesh.elems, rhoHigh, rhoLow);

        otherwise
            error('extractBinaryMask: unknown method ''%s''. Use ''fixed'', ''volume'', or ''hysteresis''.', method);
    end
end

% -------------------------------------------------------------------------
function solid = hysteresisFloodFill(x, elems, rhoHigh, rhoLow)
% Flood-fill from unconditional seeds (x >= rhoHigh) through elements
% with x >= rhoLow that are face/edge/node-adjacent to reached elements.
    nElems   = size(elems, 1);
    nPerElem = size(elems, 2);
    nNodes   = max(elems(:));

    % N(node, elem) = 1  sparse adjacency between nodes and elements
    N = sparse(elems(:), repmat((1:nElems)', nPerElem, 1), true, nNodes, nElems);

    seeds    = find(x >= rhoHigh);
    visited  = false(nElems, 1);
    visited(seeds) = true;
    frontier = seeds;

    while ~isempty(frontier)
        % Nodes belonging to frontier elements
        frontierNodes = find(any(N(:, frontier), 2));
        % Elements sharing at least one node with the frontier
        adjElems = find(any(N(frontierNodes, :), 1))';
        % Admit unvisited elements that meet the lower threshold
        newElems = adjElems(~visited(adjElems) & x(adjElems) >= rhoLow);
        visited(newElems) = true;
        frontier = newElems;
    end
    solid = visited;
end

% -------------------------------------------------------------------------
function v = optField(s, field, default)
    if isfield(s, field), v = s.(field); else, v = default; end
end
