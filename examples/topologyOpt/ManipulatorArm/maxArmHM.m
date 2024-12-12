function maxHM = maxArmHM(E,nu,segmentLength,R,r,res, alpha, beta, ShapeFn)
        model = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, beta, ShapeFn);
        x=ones(model.analysis.getTotalElemsNumber(),1);
        model.analysis.solveWeighted(x);
        model.analysis.computeElementResults();
        maxHM = max(model.fe.results.nodal.all(:,13));
end

