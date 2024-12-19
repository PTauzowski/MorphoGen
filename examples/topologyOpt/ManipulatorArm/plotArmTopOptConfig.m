function xtop_full = plotArmTopOptConfig(Rfilter, analysis, xopt, cutTreshold, penal, false)
    figure;
    topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, false );
    xtop_full=xopt;
    nArms=size(analysis.mesh.elems,1)/size(xopt,1)/2;
    for k=1:nArms-1
        xtop_full=[xtop_full; flip(xopt); xopt ];
    end
    xtop_full=[xtop_full; flip(xopt)];
    topOpt.x=xtop_full;
    topOpt.plotCurrentFrame();
end

