function [maxHM_solid, maxDisp_solid, maxHM_top, maxDisp_top] =plotSegmentResultsStep8(fileName, scheme_title, model, topOpt, volFr)

    outDir = "figs/step8";

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fileName = fullfile(outDir, fileName);

    topOpt.setFrame(topOpt.findFrame(volFr)); 
    x = topOpt.x;
   
    xoptFull = model.segmentToArm(x);

    xOnes=xoptFull;
    xOnes(:)=1;

    elemsInTop =  xoptFull > 0.5;
   

    figure, hold on, axis on; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");

    light('Position', [-1 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');

    topOpt.x=xOnes;

    model.plot(); 
    view(45,25); 
    title( scheme_title + "arm model" );

    exportgraphics(gcf, fileName + "_arm_model.png", 'Resolution', 600); 
    savefig(gcf, fileName + "__arm_model.fig");

    


    figure, hold on, axis on; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");

    light('Position', [-1 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');

    qfem_solid = model.analysis.solveWeighted(xOnes);
    model.analysis.computeElementResults(xOnes);
    maxHM_solid = max(model.fe.results.gp.all(:,13));
    maxDisp_solid = max(abs(qfem_solid));
    
    model.analysis.felems{1}.selectedElems=[];
    model.analysis.plotMaps(["sHM"],0.0); 

    view(45,25); 
    title("von Mises stress " + ", HM_{max}=" +  sprintf('%.2f',maxHM_solid/1000) + " kPa, u_{max} = " + sprintf('%.2G', maxDisp_solid*1000)+" mm");

    exportgraphics(gcf, fileName + "_arm_model_HM.png", 'Resolution', 600); 
    savefig(gcf, fileName + "__arm_model_HM.fig");




    topOpt.setFrame(topOpt.findFrame(volFr)); 

    figure, hold on, axis on; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");

    light('Position', [-1 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');

    model.analysis.felems{1}.plotSolidSelected(model.mesh.nodes,elemsInTop);

    view(45,25); 
    title( scheme_title + "Topology of manipulator" );

    exportgraphics(gcf, fileName + "_arm_topology.png", 'Resolution', 600); 
    savefig(gcf, fileName + "__arm_topology.fig");





    figure, hold on, axis on; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");

    light('Position', [-1 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');
    
    
    qfem_top = model.analysis.solveWeighted(xoptFull);
    model.analysis.computeElementResults(xoptFull);
    maxHM_top = max(max(model.fe.results.gp.all(13,elemsInTop,:)));
    maxDisp_top = max(abs(qfem_top));


    model.analysis.felems{1}.selectedElems=elemsInTop;
    model.analysis.plotMaps(["sHM"],0.0); 

    model.analysis.felems{1}.selectedElems=[];
    view(45,25); 
    title( scheme_title + "von Mises stress, HM_{max}=" +  sprintf('%.2f',maxHM_top/1000) + " kPa, u_{max} = " + sprintf('%.2G', maxDisp_top*1000)+" mm");

    exportgraphics(gcf, fileName + "_arm_topologyHM.png", 'Resolution', 600); 
    savefig(gcf, fileName + "_arm_topologyHM.fig");





    topOpt.setFrame(topOpt.findFrame(volFr)); 
    x = topOpt.x;
    xOnes=x;
    xOnes(:)=1;

    figure, hold on, axis on; 
    daspect([1 1 1]);
    xlabel("x");
    ylabel("y");
    zlabel("z");

    light('Position', [-1 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');

    topOpt.x=xOnes;

    edge_color=model.fe.edge_color;
    model.fe.edge_color='None';
    face_color=model.analysis.felems{1}.face_color;

    qnodal_solid = model.analysis.fromFEMVector(qfem_solid);
    qnodal_top = model.analysis.fromFEMVector(qfem_top);

    model.analysis.felems{1}.plotSolidDeformed(model.mesh.nodes,qnodal_solid,0.0,1);
    model.analysis.felems{1}.face_alpha=0.2;

    model.analysis.felems{1}.plotSolidDeformed(model.mesh.nodes,qnodal_solid,-0.5,1);
    model.analysis.felems{1}.face_color=[0.8 0.0 0.0];

    model.analysis.felems{1}.selectedElems=elemsInTop;
    model.analysis.felems{1}.plotSolidDeformed(model.mesh.nodes,qnodal_top,-0.5,elemsInTop);
    model.analysis.felems{1}.selectedElems=[];
    view(45,25); 
    title( scheme_title + "arm model" );

    model.fe.face_alpha=1.0;

    exportgraphics(gcf, fileName + "_arm_model_deformed.png", 'Resolution', 600); 
    savefig(gcf, fileName + "__arm_model_deformed.fig");
    model.analysis.felems{1}.face_color=face_color;
    model.fe.edge_color=edge_color;
     
end