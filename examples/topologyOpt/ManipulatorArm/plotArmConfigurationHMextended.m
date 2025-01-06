function  plotArmConfigurationHMextended(plotTitle,pdfFileName,E,nu,segmentLength,R,r,res, halfSegmentNelems, nthSegment, alpha, sample, ShapeFn, xOpt, dispFactor)
        modelSort = ManipulatorModel3D(E,nu,segmentLength,R,r,res, alpha, sample, ShapeFn,true);

        %modelSort.fe.plot(modelMin.mesh.nodes);
        %line(endPoints(:,1),endPoints(:,2),endPoints(:,3),Marker=".",Color='r',LineStyle='none')
        
        x=ones(modelSort.analysis.getTotalElemsNumber(),1);     
       

        fontSize=24;
        hmIndex=14;
        
        modelSort.analysis.solveWeighted(x);
        modelSort.analysis.computeElementResults();
        maxHM = max(modelSort.fe.results.nodal.all(:,13));

        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(0,0);
        title(plotTitle + ", xz view, maxHM=" + num2str(maxHM));
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
        title(plotTitle + ", yz view" + num2str(maxHM));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 

        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor); C = modelSort.fe.results.nodal.all(:,17);    
        view(90,90);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel(plotTitle + ", xy view" + num2str(maxHM));
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
        title(plotTitle + ", pan view" + num2str(maxHM));
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
        title(plotTitle + ", one segment pan view" + num2str(maxHM));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 
        modelSort.fe.selectedElems=[];

% Result for topology

        nArms=size(modelSort.analysis.mesh.elems,1)/size(xOpt,1)/2;
        x=xOpt;

        for k=1:nArms-1
            x=[x; flip(xOpt); xOpt ];
        end
        x=[x; flip(xOpt)];
        
        elemsInTop = x>0.5;
        nodesInTop = false(size(modelSort.mesh.nodes,1));
        nodesInTop( modelSort.mesh.elems(elemsInTop,:) ) = true;
        
        modelSort.analysis.solveWeighted(x);
        modelSort.analysis.computeElementResults();

        [maxHM, iMaxNode] = max(modelSort.fe.results.nodal.all(nodesInTop,13));

        modelSort.fe.selectedElems=elemsInTop;
        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(0,0);
        title(plotTitle + ", xz view, maxHM=" + num2str(maxHM));
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
        title(plotTitle + ", yz view" + num2str(maxHM));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 

        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor); C = modelSort.fe.results.nodal.all(:,17);    
        view(90,90);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel(plotTitle + ", xy view" + num2str(maxHM));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 


        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        hold on;
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(45,45);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(plotTitle + ", pan view" + num2str(maxHM));
        ax = gca; 
        ax.FontSize = fontSize;  
        plot3(modelSort.mesh.nodes(iMaxNode,1),modelSort.mesh.nodes(iMaxNode,2),modelSort.mesh.nodes(iMaxNode,3),"Marker","o","Color","m");
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 


        elemsInSegment=elemsInTop;
        elemsInSegment(:)=false;
        elemsInSegment((halfSegmentNelems+1:3*halfSegmentNelems)' + segOffset) = true;

        modelSort.fe.selectedElems=find(elemsInSegment(x>0.5));
        modelSort.analysis.plotMaps(["sHM"],dispFactor);
        %modelSort.fe.plotWired(modelSort.mesh.nodes,modelSort.analysis.qnodal,dispFactor);
        view(45,45);
        axis on;
        xlabel("x");
        ylabel("y");
        zlabel("z");
        title(plotTitle + ", one segment pan view" + num2str(maxHM));
        ax = gca; 
        ax.FontSize = fontSize;  
        exportgraphics(gcf, pdfFileName, 'Append', true, 'Resolution',600); 
        modelSort.fe.selectedElems=[];

        % close(gcf);

end

