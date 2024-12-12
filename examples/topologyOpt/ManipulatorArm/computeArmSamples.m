function [maxHM, endPoints] = computeArmSamples(E,nu,segmentLength,R,r,res, alpha, samples, ShapeFn)
    
    nSamples=size(samples,1);
    maxNM=zeros(nSamples,1);    
    endPoints=zeros(nSamples,3);
    ShapeFn = ShapeFunctionL8;

    tic;
    for k=1:nSamples
        betas=samples(k,:);
    
        model = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, betas, ShapeFn);
        model.analysis.printProblemInfo();
    
        x=ones(model.analysis.getTotalElemsNumber(),1);
        model.analysis.solveWeighted(x);
        model.analysis.computeElementResults();
        maxHM(k)= max(model.fe.results.nodal.all(:,13));
        endPoints(k,:)=model.xEnd;
    end

    disp(['Average single analysis time: ' num2str(toc/nSamples)]);

end

