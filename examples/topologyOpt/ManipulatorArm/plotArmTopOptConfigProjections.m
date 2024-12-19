function xtop_full = plotArmTopOptConfigProjections(pdfFileName,plotTitle,Rfilter, analysis, halfSegmentNelems, nthSegment, xopt, cutTreshold, penal, false)
    figure;
    topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, false );
    xtop_full=xopt;
    analysis.mesh.exportMeshToFile( xopt>0.5, pdfFileName);
    nArms=size(analysis.mesh.elems,1)/size(xopt,1)/2;
    fontSize=18;
    for k=1:nArms-1
        xtop_full=[xtop_full; flip(xopt); xopt ];
    end
    xtop_full=[xtop_full; flip(xopt)];
    topOpt.x=xtop_full;
    topOpt.plotCurrentFrame();
    view(0,0);
    axis on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    title(strcat(plotTitle , ", xz view"));
    light('Position', [-1 -1 5], 'Style', 'local');
    light('Position', [1 -1 5], 'Style', 'infinite');
    ax = gca; 
    ax.FontSize = fontSize;  
    exportgraphics(ax, pdfFileName+"_XZ.png", 'Resolution',600); 

    figure;
    topOpt.plotCurrentFrame();
    view(90,0);
    axis on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    title(strcat(plotTitle , ", yz view"));
    light('Position', [1 -1 5], 'Style', 'local');
    light('Position', [1 -1 5], 'Style', 'infinite');
    ax = gca; 
    ax.FontSize = fontSize;  
    exportgraphics(gcf, pdfFileName+"_YZ.png",  'Resolution',600); 

    figure;
    topOpt.plotCurrentFrame();
    view(90,90);
    axis on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    title(strcat(plotTitle , ", xy view"));
    light('Position', [-1 -2 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');
    ax = gca; 
    ax.FontSize = fontSize;  
    exportgraphics(gcf, pdfFileName+"_XY.png", 'Resolution', 600); 

    figure;
    topOpt.plotCurrentFrame();
    view(45,45);
    axis on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    title(strcat(plotTitle , ", pan view"));
    light('Position', [-1 -1 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');
    ax = gca; 
    ax.FontSize = fontSize;  
    exportgraphics(gcf, pdfFileName+"_PAN.png",  'Resolution', 600); 
    
    segOffset = (nthSegment-1) * 2 * halfSegmentNelems;
    analysis.felems{1}.selectedElems=(halfSegmentNelems+1:3*halfSegmentNelems)'+segOffset;
    figure;
    topOpt.plotCurrentFrame();
    view(45,45);
    axis on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    title(strcat(plotTitle , ", pan view"));
    light('Position', [-1 -1 5], 'Style', 'local');
    light('Position', [1 1 5], 'Style', 'infinite');
    ax = gca; 
    ax.FontSize = fontSize;  
    exportgraphics(gcf, pdfFileName+"_Segment.png",  'Resolution', 600); 
    analysis.felems{1}.selectedElems=[];
end

