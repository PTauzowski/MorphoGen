function exportDensitySTL(model, rho, filename, threshold, nModules)
% exportDensitySTL  Write density-thresholded topology as binary STL.
%
%   exportDensitySTL(model, rho, filename)
%   exportDensitySTL(model, rho, filename, threshold)
%   exportDensitySTL(model, rho, filename, threshold, nModules)
%
%   Solid elements (rho > threshold) are extracted, their exterior boundary
%   quad faces are identified, split into triangles, and written as a binary
%   STL. For a smooth isosurface use exportDensitySmoothSTL instead.

    if nargin < 4 || isempty(threshold), threshold = 0.5; end
    if nargin < 5 || isempty(nModules),  nModules  = 0;   end

    [V, F] = densityVoxelSurface(model, rho, threshold, nModules);

    if isempty(F)
        warning('exportDensitySTL:noSurface', ...
            'No elements exceed threshold %.2f; %s not written.', threshold, filename);
        return;
    end

    writeBinarySTL(filename, V, F);
end

% -------------------------------------------------------------------------
function writeBinarySTL(filename, V, F)
    nTri    = size(F, 1);
    v1      = V(F(:,1), :);
    v2      = V(F(:,2), :);
    v3      = V(F(:,3), :);
    normals = cross(v2 - v1, v3 - v1, 2);
    normals = normals ./ max(sqrt(sum(normals.^2, 2)), eps);

    fid = fopen(filename, 'wb');
    if fid < 0
        error('exportDensitySTL:cannotOpenFile', ...
            'Cannot open file for writing: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, sprintf('%-80s', 'Binary STL exported by exportDensitySTL'), 'char');
    fwrite(fid, uint32(nTri), 'uint32');
    for i = 1:nTri
        fwrite(fid, single(normals(i,:)), 'single');
        fwrite(fid, single(v1(i,:)),      'single');
        fwrite(fid, single(v2(i,:)),      'single');
        fwrite(fid, single(v3(i,:)),      'single');
        fwrite(fid, uint16(0),            'uint16');
    end
end
