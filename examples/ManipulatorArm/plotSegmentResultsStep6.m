function plotSegmentResultsStep6(fileName, scheme_title, topOpt, volFr)

    outDir = "figs/step6";

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fileName = fullfile(outDir, fileName);


    figure, hold on, axis on; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");

    light('Position', [-1 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');
            
    fontSize=16;
    hmIndex=14;

    x = topOpt.x;
    x(:)=1;
    
    qfem = topOpt.FEAnalyses(1).solveWeighted(x);
    topOpt.FEAnalyses(1).computeElementResults();
    maxHM = max(topOpt.FEAnalyses(1).felems{1}.results.nodal.all(:,13));
    maxDisp = max(abs(qfem));
    
    topOpt.FEAnalyses(1).plotMaps(["sHM"],0.0);
    %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
    title("Segment model " + ", xz view, HM_{max}=" +  sprintf('%.2f',maxHM/1000) + " kPa, u_{max} = " + sprintf('%.2G', maxDisp*1000)+" mm");

    exportgraphics(gcf, fileName+"_segmentHM.png", 'Resolution', 600); 
    savefig(gcf, fileName+"_segmentHM.fig");

    view(45,25); 
    exportgraphics(gcf, fileName+"_model.png", 'Resolution', 600); 
    savefig(gcf, fileName+"_model.fig");
    
   
     
end