function  developResultsForTopology(filename,title,model,R,r, alpha, res,segmentLength,segmentNo,xOpt,betas,dispFactor,Rfilter,cutTreshold,penal)

    plotArmTopOptConfigProjections(filename, title, Rfilter, model.analysis, model.halfSegmentNelems, segmentNo, xOpt, cutTreshold, penal, false);
    %plotArmConfigurationHMextended(title,filename+".pdf",model.fe.mat.E,model.fe.mat.nu,segmentLength,R,r,res, model.halfSegmentNelems, segmentNo, alpha, betas, model.fe.sf, xOpt, dispFactor);
    
end