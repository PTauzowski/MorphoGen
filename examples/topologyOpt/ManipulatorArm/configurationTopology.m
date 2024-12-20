function [xopt_bending, xopt_bending_buckling, lambda1, lambda2 ]  = configurationTopology(E,nu,R,r,segmentLength,ShapeFn,alpha,sample,frameElems,loadFactor)
     model= ManipulatorModel3D(E,nu,segmentLength,R,r,1, alpha, sample, ShapeFn, false);
     mesh=Mesh();
     frameElem=Frame3D(frameElems,E,0.02,0.8*E,0.0004,0.0004,0.003);
     mesh.nodes=model.frameNodes;
     barNumber=2;
     [Fel, ~] = computeInternalForces(frameElem,mesh);
     N  = Fel(7,barNumber)*loadFactor;
     Tz = Fel(8,barNumber)*loadFactor;
     Ty = Fel(9,barNumber)*loadFactor;
     Ms = Fel(10,barNumber)*loadFactor;
     Mz = Fel(11,barNumber)*loadFactor;
     My = Fel(12,barNumber)*loadFactor;

    [xopt_bending, xopt_bending_buckling, lambda1, lambda2 ] = ArmTopOptBucklingFn('maxHM',R,r,segmentLength,alpha,Ty,Tz,N,My,Mz,Ms);
end

