function [Fel, Feg] = computeInternalForces(frameElem,mesh)
    analysis = LinearElasticityWeighted( frameElem, mesh, false );
    analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 0 0 0 0]);
    analysis.loadClosestNode(mesh.nodes(end,:), ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 -1 0 0 0] );
    x=ones(size(frameElem.elems,1),1);
    analysis.solveWeighted(x);
    [Fel, Feg] = frameElem.computeResults(mesh.nodes,analysis.qnodal);
    % analysis.plotCurrentLoad();
    % analysis.plotSupport();
end

