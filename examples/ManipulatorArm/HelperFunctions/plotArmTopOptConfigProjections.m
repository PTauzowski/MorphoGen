function xtop_full = plotArmTopOptConfigProjections(filename,plotTitle,Rfilter, analysis, halfSegmentNelems, nthSegment, xopt, cutTreshold, penal, is_const)

    pdfFileName = "figs/"+filename;
    topOpt = StressIntensityTopologyOptimizationVol( Rfilter, analysis, cutTreshold, penal, 0.4, is_const );
    xtop_full=xopt;
    analysis.computeElementResults();
    analysis.mesh.exportMeshToFile( xopt>0.5, pdfFileName+"_"+plotTitle);
    nArms=size(analysis.mesh.elems,1)/halfSegmentNelems;
    fontSize=18;
    if numel(xopt)==1
        xtop_full=ones(size(analysis.mesh.elems,1),1);
    else
        xtop_full=xopt;
        for k=1:(nArms-2)/2
            xtop_full=[xtop_full; flip(xopt); xopt ];
        end
        xtop_full=[xtop_full; flip(xopt)];
    end
   
    topOpt.x=xtop_full;   
    
    % figure;
    % topOpt.plotCurrentFrame();
    % view(0,0);
    % axis on;
    % xlabel("x");
    % ylabel("y");
    % zlabel("z");
    % title(strcat(plotTitle , ", xz view"));
    % light('Position', [-1 -1 5], 'Style', 'local');
    % light('Position', [1 -1 5], 'Style', 'infinite');
    % ax = gca; 
    % ax.FontSize = fontSize;  
    % exportgraphics(ax, pdfFileName+"_XZ.png", 'Resolution',600); 
    % exportgraphics(ax, pdfFileName+".pdf", 'Resolution',600); 
    % 
    % figure;
    % topOpt.plotCurrentFrame();
    % view(90,0);
    % axis on;
    % xlabel("x");
    % ylabel("y");
    % zlabel("z");
    % title(strcat(plotTitle , ", yz view"));
    % light('Position', [1 -1 5], 'Style', 'local');
    % light('Position', [1 -1 5], 'Style', 'infinite');
    % ax = gca; 
    % ax.FontSize = fontSize;  
    % exportgraphics(gcf, pdfFileName+"_YZ.png",  'Resolution',600); 
    % exportgraphics(gcf, pdfFileName+".pdf", 'Append', true, 'Resolution',600);

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
    exportgraphics(gcf, pdfFileName+"_"+plotTitle+"_XY.png", 'Resolution', 600); 
    savefig(gcf, pdfFileName+"_"+plotTitle+"_XY.fig");
    %exportgraphics(gcf, pdfFileName+".pdf", 'Append', true, 'Resolution',600);

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
    exportgraphics(gcf, pdfFileName+"_"+plotTitle+"_PAN.png",  'Resolution', 600); 
    savefig(gcf, pdfFileName+"_"+plotTitle+"_PAN.fig"); 
    %exportgraphics(gcf, pdfFileName+".pdf", 'Append', true, 'Resolution',600);
    
    segOffset = (nthSegment-1) * 2 * halfSegmentNelems - halfSegmentNelems;
    if segOffset<0
        segOffset=0;
    end
    if nthSegment==1
        analysis.felems{1}.selectedElems=(1:halfSegmentNelems)';
    elseif nthSegment==7
        analysis.felems{1}.selectedElems=(1:halfSegmentNelems)'+segOffset;
    else
        analysis.felems{1}.selectedElems=(1:2*halfSegmentNelems)'+segOffset;
    end
    xSeg = xtop_full(analysis.felems{1}.selectedElems);
    
    segmentMesh=Mesh();
    segmentMesh.mergeMesh( analysis.mesh );
    segmentMesh.leaveElemsByNumbers(analysis.felems{1}.selectedElems);

    segmentTopMesh=Mesh();
    segmentTopMesh.mergeMesh( segmentMesh );
    segmentTopMesh.leaveElemsByNumbers(xSeg>0.5);

    [nodes, elems] = segmentMesh.getTetrahedralMesh(xSeg>0.5);
    save(pdfFileName+"_"+plotTitle+"_tetramesh.mat","elems","nodes");
    trigMesh = Mesh();
    trigMesh.nodes=nodes;
    trigMesh.elems=elems;
    % trigMesh.exportToPLY(pdfFileName); 
     trigMesh.exportTetraToSTL(pdfFileName);
    % trigMesh.exportToStep(pdfFileName);

    [nodes, elems] = analysis.mesh.getTetrahedralMesh(xtop_full>0.5);
    trigMesh = Mesh();
    trigMesh.nodes=nodes;
    trigMesh.elems=elems;
    %trigMesh.exportToPLY(pdfFileName+"_fullArm");

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
    exportgraphics(gcf, pdfFileName+"_"+plotTitle+"_Segment.png",  'Resolution', 600); 
    %exportgraphics(gcf, pdfFileName+".pdf", 'Append', true, 'Resolution',600);
    savefig(gcf, pdfFileName+"_"+plotTitle+"_Segment.fig");
    analysis.felems{1}.selectedElems=[];
end

