function plotSegmentResultsStep7(fileName, scheme_title , topOpt, loadFactor, volFr)

    outDir = "figs/step7";

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fileName = fullfile(outDir, fileName);

    figure, hold on, axis on; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");

    light('Position', [-2 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');

    topOpt.setFrame(topOpt.findFrame(volFr)); 
    x = topOpt.x;
    
    qfem = topOpt.FEAnalyses(1).solveWeighted(x);
    topOpt.FEAnalyses(1).computeElementResults(x);
   
    topOpt.plotCurrentFrame(); 
    view(45,25); 
    title( scheme_title );
    
    exportgraphics(gcf, fileName + "_half_segment_topology.png", 'Resolution', 600); 
    savefig(gcf, fileName + "_half_segment_topology.fig");

    elemsInTop = x>0.5;

    figure, hold on, axis on; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");

    light('Position', [-2 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');

    topOpt.setFrame(topOpt.findFrame(volFr)); 
    x = topOpt.x;
    
    qfem = topOpt.FEAnalyses(1).solveWeighted(x);
    topOpt.FEAnalyses(1).computeElementResults(x);
   
    model.fe.face_alpha=0.1;

    topOpt.FEAnalyses(1).felems{1}.plotSolidSelected(topOpt.FEAnalyses(1).mesh.nodes,elemsInTop);
    topOpt.FEAnalyses(1).felems{1}.plotSolidSelected([topOpt.FEAnalyses(1).mesh.nodes(:,1) topOpt.FEAnalyses(1).mesh.nodes(:,2) -topOpt.FEAnalyses(1).mesh.nodes(:,3)],elemsInTop);
    view(45,25); 
    title( scheme_title );
    
    exportgraphics(gcf, fileName + "_segment_topology.png", 'Resolution', 600); 
    savefig(gcf, fileName + "_segment_topology.fig");


    for k=1:numel(topOpt.FEAnalyses)
        figure, hold on, axis on; 
        daspect([1 1 1]);
        xlabel("x");
        ylabel("y");
        zlabel("z");

        light('Position', [-1 -2 5], 'Style', 'local');
        light('Position', [1 1 5], 'Style', 'infinite');

        topOpt.setFrame(topOpt.findFrame(volFr)); 
        x = topOpt.x;

        qfem = topOpt.FEAnalyses(k).solveWeighted(x);
        topOpt.FEAnalyses(k).computeElementResults(x);

        %nodesInTop = false(size(topOpt.FEAnalyses(k).mesh.nodes,1),1);

        maxHM = max(max(topOpt.FEAnalyses(k).felems{1}.results.gp.all(13,elemsInTop,:))) /loadFactor;
        maxDisp = max(abs(qfem)) / loadFactor;

        topOpt.FEAnalyses(k).felems{1}.selectedElems=elemsInTop;

        topOpt.FEAnalyses(k).plotMaps(["sHM"],0.0); 
        topOpt.FEAnalyses(k).felems{1}.selectedElems=[];
        view(45,25); 
        title( scheme_title + "for config. max Q" + num2str(k) + ", xz view, HM_{max}=" +  sprintf('%.2f',maxHM/10^6) + " MPa,  u_{max} = " + sprintf('%.2G', maxDisp*1000)+" mm");

        exportgraphics(gcf, fileName + "_Q" + num2str(k) + "_segment_topologyHM.png", 'Resolution', 600); 
        savefig(gcf, fileName + "_Q" + num2str(k) + "_segment_topologyHM.fig");
    end
     
end