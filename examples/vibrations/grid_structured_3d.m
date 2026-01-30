function grid = grid_structured_3d(nelx, nely, nelz, EL, EW, EH)
%GRID_STRUCTURED_3D Structured HEX8 mesh helpers (Appendix G style)
% Returns:
%   nodMat      (nely+1) x (nelx+1) x (nelz+1) node IDs (int32)
%   edofMat     nEle x 24 DOF connectivity
%   eleNodesID  nEle x 8 node IDs
%   sI, sII     assembly index helpers
%   Iar0        lower-tri index pairs for full domain
%   LSgrid      fields x,y,z (node coords reshaped for calc_Phi)
%   volNod      node volume weights (1/8 per HEX)
%
% Node numbering follows MATLAB meshgrid style in Appendix G.

    nEle = nelx * nely * nelz;
    nNod = (nelx+1) * (nely+1) * (nelz+1);
    nNodfc = (nelx+1) * (nely+1);

    nodMat = int32(reshape(1:nNod, nely+1, nelx+1, nelz+1));

    edofVec = reshape(3 * nodMat(1:nely, 1:nelx, 1:nelz), nEle, 1);
    edofMat = edofVec + int32([3*nNodfc + [3*nelx + [4 5 6 1 2 3]] -2 -1 0 1 2 3]);

    eleNodesID = edofMat(:, 3:3:24) ./ 3;

    [sI, sII] = deal([]);
    for j = 1:24
        sI  = cat(2, sI, j:24);
        sII = cat(2, sII, repmat(j, 1, 24-j+1));
    end

    [iK, jK] = deal(edofMat(:, sI)', edofMat(:, sII)');
    Iar0 = sort([iK(:), jK(:)], 2, 'descend');

    [x, y, z] = meshgrid(EL * (-nelx/2:nelx/2), EW * (-nely/2:nely/2), EH * (-nelz/2:nelz/2));
    LSgrid.x = permute(x, [2,1,3]);
    LSgrid.y = permute(y, [2,1,3]);
    LSgrid.z = permute(z, [2,1,3]);

    volNod = sparse(double(eleNodesID(:)), 1, 1/8);

    grid.nodMat = nodMat;
    grid.edofMat = edofMat;
    grid.eleNodesID = eleNodesID;
    grid.sI = sI;
    grid.sII = sII;
    grid.Iar0 = Iar0;
    grid.LSgrid = LSgrid;
    grid.volNod = volNod;
end

