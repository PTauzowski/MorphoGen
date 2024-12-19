function  plotArmConfigurationHMextended(plotTitle,pdfFileName,E,nu,segmentLength,R,r,res, halfSegmentNelems, nthSegment, alpha, sample, ShapeFn, dispFactor)
        modelSort = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sample, ShapeFn,true);
        %modelSort.fe.plot(modelMin.mesh.nodes);
        %line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none');
        x=ones(modelSort.analysis.getTotalElemsNumber(),1);
        fontSize=24;
        
        modelSort.analysis.solveWeighted(x);
        modelSort.analysis.computeElementResults();

        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(0,0);
        title(strcat(plotTitle , ", xz view"));
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Resolution',600); 

        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(90,0);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(strcat(plotTitle , ", yz view"));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 

        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(90,90);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(strcat(plotTitle , ", xy view"));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 


        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(45,45);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(strcat(plotTitle , ", pan view"));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 

        segOffset = (nthSegment-1) * 2 * halfSegmentNelems;
        modelSort.fe.selectedElems=(halfSegmentNelems+1:3*halfSegmentNelems)' + segOffset;
        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(45,45);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(strcat(plotTitle , ", one segment pan view"));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 
        modelSort.fe.selectedElems=[];
        % close(gcf);

end

