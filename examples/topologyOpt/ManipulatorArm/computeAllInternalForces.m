function [vN, vTz, vTy, vMs, vMz, vMy] = computeAllInternalForces(frameNodes,frameElems,E,nu)
    mesh=Mesh();
    mesh.nodes=frameNodes(:,:,1);
    frameElem=Frame3D(frameElems,E,0.02,0.8*E,0.0004,0.0004,0.003);
    nSamples = size( frameNodes, 3 );
    analysis = LinearElasticityWeighted( frameElem, mesh, false );
    analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 0 0 0 0]);
    analysis.loadClosestNode(mesh.nodes(end,:), ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 -1 0 0 0] );
    x=ones(size(frameElem.elems,1),1);
    vN=zeros( nSamples, 2 );
    barNumber=2;
    for k=1:nSamples
        mesh.nodes=frameNodes(:,:,k);
        analysis.solveWeighted(x);
        [Fel, ~] = frameElem.computeResults(mesh.nodes,analysis.qnodal);
        vN(k,1)  = Fel(1,barNumber);
        vTz(k,1) = Fel(2,barNumber);
        vTy(k,1) = Fel(3,barNumber);
        vMs(k,1) = Fel(4,barNumber);
        vMz(k,1) = Fel(5,barNumber);
        vMy(k,1) = Fel(6,barNumber);

        vN(k,2)  = Fel(7,barNumber);
        vTz(k,2) = Fel(8,barNumber);
        vTy(k,2) = Fel(9,barNumber);
        vMs(k,2) = Fel(10,barNumber);
        vMz(k,2) = Fel(11,barNumber);
        vMy(k,2) = Fel(12,barNumber);
    end
end

