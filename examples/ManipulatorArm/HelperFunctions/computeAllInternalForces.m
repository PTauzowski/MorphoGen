function [vN, vTy, vTz, vMs, vMy, vMz, all_forces] = computeAllInternalForces(frameNodes,frameElems,E,nu,R,r,betas, Pz)
    mesh=Mesh();
    mesh.nodes=frameNodes(:,:,1);
    frameElem=Frame3D(frameElems,E,nu, R,r);
    nSamples = size( frameNodes, 3 );
    analysis = LinearElasticityWeighted( frameElem, mesh, false );
    analysis.fixClosestNode([0 0 0], ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 0 0 0 0]);
    analysis.loadClosestNode(mesh.nodes(end,:), ["ux" "uy" "uz" "fix" "fiy" "fiz"], [0 0 -Pz 0 0 0] );
    x=ones(size(frameElem.elems,1),1);
    vN  = zeros( nSamples, 2, 6 );
    vTz = zeros( nSamples, 2, 6 );
    vTy = zeros( nSamples, 2, 6 );
    vMs = zeros( nSamples, 2, 6 );
    vMz = zeros( nSamples, 2, 6 );
    vMy = zeros( nSamples, 2, 6 );
    all_forces = zeros( nSamples, 2, 6, 6 );

    for k=1:nSamples
            mesh.nodes=frameNodes(:,:,k);
            frameElem.betas=betas(k,:);
            analysis.solveWeighted(x);
            [Fel, ~] = frameElem.computeResults(mesh.nodes,analysis.qnodal);
            for l=1:6
                vN(k,1,l)  = Fel(1,l);
                vTy(k,1,l) = Fel(2,l);
                vTz(k,1,l) = Fel(3,l);
                vMs(k,1,l) = Fel(4,l);
                vMy(k,1,l) = Fel(5,l);
                vMz(k,1,l) = Fel(6,l);
        
                vN(k,2,l)  = Fel(7,l);
                vTy(k,2,l) = Fel(8,l);
                vTz(k,2,l) = Fel(9,l);
                vMs(k,2,l) = Fel(10,l);
                vMy(k,2,l) = Fel(11,l);
                vMz(k,2,l) = Fel(12,l);
                all_forces(k,1,l,:)=Fel(1:6,l);
           end
    end
end

