function  plotArmConfigurationHMextended(plotTitle,pdfFileName,E,nu,segmentLength,R,r,res, halfSegmentNelems, alpha, sample, ShapeFn)
        modelSort = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sample, ShapeFn);
        %modelSort.fe.plot(modelMin.mesh.nodes);
        %line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
        x=ones(modelSort.analysis.getTotalElemsNumber(),1);
        modelSort.analysis.solveWeighted(x);
        modelSort.analysis.computeElementResults();

        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(0,0);
        title(strcat(plotTitle , ", xz view"));
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        exportgraphics(gcf, pdfFileName, 'Resolution',300); 

        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(90,0);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(strcat(plotTitle , ", yz view"));
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',300); 
        
        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(90,90);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(strcat(plotTitle , ", xy view"));
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',300); 
        

        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(45,45);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(strcat(plotTitle , ", pan view"));
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',300); 

        
        modelSort.fe.selectedElems=(halfSegmentNelems+1:3*halfSegmentNelems)';
        modelSort.analysis.plotMaps(["sHM"],0.1);
        modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,0.1);
        view(45,45);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(strcat(plotTitle , ", one segment pan view"));
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',300); 
        modelSort.fe.selectedElems=[];
        %close(fig);

end

