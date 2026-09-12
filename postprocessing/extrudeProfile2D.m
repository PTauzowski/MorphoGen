function [V, F, faceLabels] = extrudeProfile2D(polygons, polyLabels, thickness)
% extrudeProfile2D  Extrude closed 2D polygons into a closed 3D solid mesh.
%
%   [V, F, faceLabels] = extrudeProfile2D(polygons, polyLabels, thickness)
%
%   polygons    {K x 1} cell array of [Ni x 2] closed polygon vertex coords
%   polyLabels  [K x 1] integer label per polygon (0 = unlabelled)
%   thickness   extrusion depth in Z (default 0.05 * bounding-box diagonal)
%
%   Returns a triangulated closed 3D surface:
%     top cap (z = thickness) + bottom cap (z = 0) + side walls.
%   faceLabels carries the source polygon label onto every triangle.

    if nargin < 2 || isempty(polyLabels), polyLabels = zeros(numel(polygons), 1); end
    if nargin < 3 || isempty(thickness)
        allPts = vertcat(polygons{:});
        if isempty(allPts)
            thickness = 1.0;
        else
            span      = max(allPts) - min(allPts);
            thickness = 0.05 * norm(span);
        end
    end

    polyLabels = polyLabels(:);
    V          = zeros(0, 3);
    F          = zeros(0, 3);
    faceLabels = zeros(0, 1);

    for pi = 1:numel(polygons)
        poly = polygons{pi};
        lab  = polyLabels(pi);
        if size(poly, 1) < 3, continue; end

        [Vp, Fp, Fl] = extrudeOnePoly(poly, thickness, lab);
        if isempty(Fp), continue; end
        Fp = Fp + size(V, 1);
        V          = [V;  Vp]; %#ok<AGROW>
        F          = [F;  Fp]; %#ok<AGROW>
        faceLabels = [faceLabels; Fl]; %#ok<AGROW>
    end
end

% -------------------------------------------------------------------------
function [V, F, faceLabels] = extrudeOnePoly(poly, thickness, lab)
    n    = size(poly, 1);
    Vbot = [poly, zeros(n, 1)];
    Vtop = [poly, repmat(thickness, n, 1)];
    V    = [Vbot; Vtop];   % bottom ring: 1..n,  top ring: n+1..2n

    % Side walls
    nSide = 2 * n;
    Fside = zeros(nSide, 3);
    for i = 1:n
        j          = mod(i, n) + 1;
        b1 = i;   b2 = j;
        t1 = i+n; t2 = j+n;
        Fside(2*i-1, :) = [b1, b2, t1];
        Fside(2*i,   :) = [b2, t2, t1];
    end

    % Caps via constrained Delaunay triangulation
    Fcap_bot = triangulateCap(poly, 0,         false);  % normal pointing -Z
    Fcap_top = triangulateCap(poly, n,         true);   % normal pointing +Z

    F          = [Fside; Fcap_bot; Fcap_top];
    nF         = size(F, 1);
    faceLabels = repmat(lab, nF, 1);
end

% -------------------------------------------------------------------------
function F = triangulateCap(poly, offset, flip)
    n = size(poly, 1);
    if n < 3
        F = zeros(0, 3);
        return;
    end

    % Constrained Delaunay: edges follow the polygon boundary
    edges = [(1:n)', [2:n, 1]'];
    try
        dt     = delaunayTriangulation(poly(:,1), poly(:,2), edges);
        inside = isInterior(dt);
        tris   = dt.ConnectivityList(inside, :);
    catch
        % Fallback: fan triangulation from first vertex (works for convex polygons)
        tris = [(ones(n-2,1)), (2:n-1)', (3:n)'];
    end

    % Shift to global vertex indices
    F = tris + offset;

    % Flip winding for correct outward normal direction
    if flip
        F = F(:, [1 3 2]);
    end
end
