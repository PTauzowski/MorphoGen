function cfg = bc_builders_3d(caseName, nelx, nely, nelz)
%BC_BUILDERS_3D Boundary conditions & mass for 3D benchmark (Appendix G)
% Returns struct with fields: fixNd, fixDof, fixEle, freeDof, massNd, massDof, massEle

    switch lower(caseName)
        case {'appendix-g', 'paper-3d'}
            % Four corner supports on bottom face (same as Appendix G code)
            fixNd = [1; (nelx+1)*(nely+1); nelx+1; nely*(nelx+1)+1];
            fixDof = [3*fixNd(:)-2; 3*fixNd(:)-1; 3*fixNd(:)];

            % Eight fixed solid elements around the supported corners (one layer)
            fixEle = [1; nelx; nelx*nely; nelx*(nely-1)+1; ...
                      (nelz-1)*nelx*nely+1; (nelz-1)*nelx*nely+nelx; ...
                      (nelz-1)*nelx*nely+nelx*nely; (nelz-1)*nelx*nely+nelx*(nely-1)+1];

            massNd = (nely/2)*(nelx+1) + (nelx/2) + 1 + nelz*((nelx+1)*(nely+1));
            massDof = 3*massNd-2 : 3*massNd;
            massEle = (nelx/2) + (nely/2 - 1)*nelx + (nelz - 1)*nelx*nely;

        otherwise
            error('Unknown 3D BC case: %s', caseName);
    end

    cfg.fixNd   = fixNd;
    cfg.fixDof  = unique(fixDof);
    cfg.fixEle  = unique(fixEle);
    cfg.freeDof = setdiff(1:3*(nelx+1)*(nely+1)*(nelz+1), cfg.fixDof);
    cfg.massNd  = massNd;
    cfg.massDof = massDof;
    cfg.massEle = massEle;
end

