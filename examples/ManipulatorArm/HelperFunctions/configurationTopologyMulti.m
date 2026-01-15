function [xopt_av, xopt_av_lambda, xoptBuckling_av, xoptBuckling_av_lambda, xopt_max, xopt_max_lambda , xoptBuckling_max, xoptBuckling_max_lambda, newLoadFactor, topOptAvLinear, topOptEnvLinear, topOptAvBuckling, topOptEnvBuckling ]  = configurationTopologyMulti(E,nu,R,r, res, res_thickness, segmentLength,ShapeFn,alpha,samples,frameElems,max_segments, loadFactor,Pz, Rmin,maxais,penal,volFr,is_const,ringMode)

     nSamples=size(samples,1);

     analysesLinear = [];
     analysesBuckling = [];
     for k=1:nSamples
         model= ManipulatorModel3D(E,nu,segmentLength,R,r,res, res_thickness, alpha, samples(k,:), ShapeFn, false, Pz);
         mesh=Mesh();
         frameElem=Frame3D(frameElems,E,nu,R,r);
         mesh.nodes=model.frameNodes;
         [Fel, ~] = computeInternalForces(frameElem,mesh,Pz);
         N  = Fel(1,max_segments(k))*loadFactor/frameElem.A;
         Ty = Fel(2,max_segments(k))*loadFactor/frameElem.A;
         Tz = Fel(3,max_segments(k))*loadFactor/frameElem.A;
         Ms = Fel(4,max_segments(k))*loadFactor/frameElem.A;
         My = Fel(5,max_segments(k))*loadFactor/frameElem.A;
         Mz = Fel(6,max_segments(k))*loadFactor/frameElem.A;
         [analysisLinear,  analysisWithBuckling, const_elems ] = ArmTopOptBucklingAnalyses('multisamples',R,r,res_thickness, segmentLength,alpha,N,Ty,Tz,Ms,My,Mz,ringMode);
         analysesLinear = [ analysesLinear analysisLinear ];
         analysesBuckling = [ analysesBuckling analysisWithBuckling ];
     end

     figure, hold on;
     analysisLinear.plotFiniteElements();

    topOptMultiLinear = StressIntensityMultiAvTopologyOptimization(Rmin,analysesLinear,maxais,penal,volFr,is_const);
    topOptMultiLinear.setConstElems(const_elems);
    [~, xopt_av] = topOptMultiLinear.solve();
    xopt_av_lambda=topOptMultiLinear.plLambda(end);
    topOptAvLinear = topOptMultiLinear;

     topOptMultiLinear = StressIntensityMultiMaxTopologyOptimization(Rmin,analysesLinear,maxais,penal,volFr,is_const);
     topOptMultiLinear.setConstElems(const_elems);
    [~, xopt_max] = topOptMultiLinear.solve();
    xopt_max_lambda = topOptMultiLinear.plLambda(end);
    topOptEnvLinear = topOptMultiLinear;

    coeff=1.1*max(xopt_av_lambda, xopt_max_lambda);

    newLoadFactor   = loadFactor*coeff;

     analysesLinear = [];
     analysesBuckling = [];
     for k=1:nSamples
         model= ManipulatorModel3D(E,nu,segmentLength,R,r,res, res_thickness, alpha, samples(k,:), ShapeFn, false, Pz);
         mesh=Mesh();
         frameElem=Frame3D(frameElems,E,nu,R,r);
         mesh.nodes=model.frameNodes;
         barNumber=2;
         [Fel, ~] = computeInternalForces(frameElem,mesh,Pz);
         N  = Fel(1,barNumber)*newLoadFactor;
         Ty = Fel(2,barNumber)*newLoadFactor;
         Tz = Fel(3,barNumber)*newLoadFactor;
         Ms = Fel(4,barNumber)*newLoadFactor;
         My = Fel(5,barNumber)*newLoadFactor;
         Mz = Fel(6,barNumber)*newLoadFactor;
         [analysisLinear,  analysisWithBuckling, const_elems ] = ArmTopOptBucklingAnalyses('multisamples',R,r,res_thickness,segmentLength,alpha,N,Ty,Tz,Ms,My,Mz,ringMode);
         analysesLinear = [ analysesLinear analysisLinear ];
         analysesBuckling = [ analysesBuckling analysisWithBuckling ];
     end

    topOptMultiBuckling = StressIntensityMultiAvTopologyOptimization(Rmin,analysesBuckling,maxais,penal,volFr,is_const);
    topOptMultiBuckling.setConstElems(const_elems);
    [~, xoptBuckling_av] = topOptMultiBuckling.solve();
    xoptBuckling_av_lambda = topOptMultiBuckling.plLambda(end);
    topOptAvBuckling = topOptMultiBuckling;

    topOptMultiBuckling = StressIntensityMultiMaxTopologyOptimization(Rmin,analysesBuckling,maxais,penal,volFr,is_const);
    topOptMultiBuckling.setConstElems(const_elems);
    [~, xoptBuckling_max] = topOptMultiBuckling.solve();
    xoptBuckling_max_lambda = topOptMultiBuckling.plLambda(end);
    topOptEnvBuckling = topOptMultiBuckling;

    xopt_av_lambda  = xopt_av_lambda/coeff;
    xopt_max_lambda = xopt_max_lambda/coeff;
end

