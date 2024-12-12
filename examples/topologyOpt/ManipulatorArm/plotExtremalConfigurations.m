function  plotExtremalConfigurations(E,nu,segmentLength,R,r,res,alpha,ShapeFn,nConfigs,samples,vSort,iSort)
    for k=1:nConfigs
        sortId = k;
        %sortId = nSamples/nSort*k;
        plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, samples(iSort(sortId),:), ShapeFn);
        title(['Model for minimal Huber-Mises for HMmax=' num2str(vSort(sortId))]);
    end
    nSamples=size(samples,1);
    for k=1:nConfigs
        sortId = nSamples-k+1-2;
        %sortId = nSamples/nSort*k;
        plotArmConfigurationHM(E,nu,segmentLength,R,r,res, alpha, samples(iSort(sortId),:), ShapeFn)
        title(['Model for maximal Huber-Mises for HMmax=' num2str(vSort(sortId))]);

    end
end

