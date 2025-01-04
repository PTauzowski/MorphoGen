function  plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, sample, ShapeFn)
        modelSort = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sample, ShapeFn, true);
        %modelSort.fe.plot(modelMin.mesh.nodes);
        %line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
        x=ones(modelSort.analysis.getTotalElemsNumber(),1);
        modelSort.analysis.solveWeighted(x);
        modelSort.analysis.computeElementResults();
        modelSort.analysis.plotMaps(["sHM"],0.0);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.0);
end

