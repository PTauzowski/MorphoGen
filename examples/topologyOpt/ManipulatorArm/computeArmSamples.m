function [maxHM, endPoints, frameNodes] = computeArmSamples(E,nu,segmentLength,R,r,res, alpha, samples, ShapeFn)
    
    nSamples=size(samples,1);
    maxNM=zeros(nSamples,1);    
    endPoints=zeros(nSamples,3);
    nArms=size(samples,2);
    ShapeFn = ShapeFunctionL8;
    frameNodes=zeros(nArms+1,3,nSamples);
    tic;
    for k=1:nSamples
        betas=samples(k,:);
    
        model = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, betas, ShapeFn, true);
        model.analysis.printProblemInfo();
    
        x=ones(model.analysis.getTotalElemsNumber(),1);
        model.analysis.solveWeighted(x);
        model.analysis.computeElementResults();
        maxHM(k)= max(model.fe.results.nodal.all(:,13));
        endPoints(k,:)=model.xEnd;
        frameNodes(:,:,k)=model.frameNodes;
    end

    disp(['Average single analysis time: ' num2str(toc/nSamples)]);

end

