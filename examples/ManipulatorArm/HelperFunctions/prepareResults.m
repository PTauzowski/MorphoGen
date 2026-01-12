function  [maxHM_dd, maxDisp_dd, maxHM_top, maxDisp_top] = prepareResults(filename, title, model,R,r, alpha, res, res_th, segmentLength,max_segments,xOpt,topOpt,betas,dispFactor,Rfilter,cutTreshold,penal)

    % plotArmTopOptConfigProjections(filename+"_initial_configuration",title, Rfilter, model.analysis, model.halfSegmentNelems, max_segments, 1, cutTreshold, penal, false);
    % plotArmTopOptConfigProjections(filename, title, Rfilter, model.analysis, model.halfSegmentNelems, max_segments, xOpt, cutTreshold, penal, false);
    % plotArmConfigurationHMextended(filename,title,model.fe.mat.E,model.fe.mat.nu,segmentLength,R,r,res, model.halfSegmentNelems, max_segments, alpha, betas, model.fe.sf, xOpt, dispFactor);
    [maxHM_dd, maxDisp_dd, maxHM_top, maxDisp_top] = ArmStressAndDispArm(model.fe.mat.E,model.fe.mat.nu,segmentLength,R,r,res, res_th, alpha, betas, model.fe.sf, xOpt );
    %ArmStressAndDispArmPlot( filename + "_plot_", topOpt, segmentLength,R,r,res,alpha, betas, model.fe.sf );
end