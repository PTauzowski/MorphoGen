function  plotArmConfigurationHMextended(pdfFileName,E,nu,segmentLength,R,r,res, alpha, sample, ShapeFn)
        modelSort = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sample, ShapeFn);
        %modelSort.fe.plot(modelMin.mesh.nodes);
        %line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
        x=ones(modelSort.analysis.getTotalElemsNumber(),1);
        modelSort.analysis.solveWeighted(x);
        modelSort.analysis.computeElementResults();

        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(0,0);
        title('Subplot 1: sin(x)');
        axis on;
        exportgraphics(gcf, pdfFileName, 'Resolution',300); 

        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(90,0);
        title('Subplot 2: sin(2x)');
        axis on;
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',300); 
        
        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(90,90);
        title('Subplot 3: sin(4x)');
        axis on;
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',300); 
        
        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(45,45);
        axis on;
        title('Subplot 4: sin(8x)')
        
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',300); 
        %close(fig);

end

