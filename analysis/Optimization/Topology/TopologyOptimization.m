classdef TopologyOptimization < ConstrainedOptimization
    
    properties
        FEAnalysis,
        qnodal,
        weights;
        distances;
        neighbours;
        celems;
        Rmin;
        V0;
        elem_inds;
        const_elems;
        erased_elems;
        totalFENumber;
        is_const;
        is_silent;
    end
    
    methods (Abstract)
        updateDesign(obj);
        printIterationInfo(obj);
        resetAnalysis();
    end
    
    methods
        function obj = TopologyOptimization(numberOfConstraints,Rmin,FEproblem,is_const)
            obj = obj@ConstrainedOptimization(FEproblem.getTotalElemsNumber(),numberOfConstraints);
            % Supressing singularity warning.
            warning('off','MATLAB:nearlySingularMatrix');
            obj.Rmin = Rmin;
            obj.FEAnalysis=FEproblem;
            obj.is_const=is_const;
            obj.totalFENumber = obj.FEAnalysis.getTotalElemsNumber();
            obj.const_elems=[];
            obj.erased_elems = false(obj.totalFENumber,1);
            obj.xmin(:)=0.01;
            obj.xmax(:)=1;
            obj.x(:)=1;
            obj.V0 = sum( obj.x );
            obj.elem_inds = FEproblem.getElemIndices();
            obj.createFilteringMatrix();
            colormap(jet);
            obj.is_silent=false;
        end
        
        function [objF, xopt] = solve(obj)
            obj.iteration=1;
            obj.resetAnalysis();
            while obj.isNotFinished()
                obj.updateDesign();
                obj.plotFrame();
                if ~obj.is_silent
                    obj.printIterationInfo();                
                end
                obj.iteration = obj.iteration + 1;
            end
            obj.plotCurrentFrame();
            obj.computeObjectiveFunction();
            xopt = obj.x;
            objF = obj.FobjValue;
        end

        function createFilteringCells(obj)
            problem=obj.FEProblem;
            obj.celems = cell( problem.getTotalElemsNumber(), size(problem.nodes,2) );
            xs = obj.midElems( problem );
            k=1;
            for i=1:size(problem.felems,1)
                for j=1:size(problem.felems{i}.elems,1)
                    [~,d] = dsearchn( xs(k,:), xs  );
                    obj.celems{ k, 1 } = find( d <= obj.Rf );
                    obj.celems{ k, 2 } = 1 - d( obj.celems{ k, 1 } ) / obj.Rf;
                    k=k+1;
                end
            end
        end
        function createFilteringCells2(obj)
            problem=obj.FEAnalysis;
            celemsi = cell( problem.getTotalElemsNumber(), size(problem.mesh.nodes,2) );
            xs = obj.midElems( problem );
            np = size(xs,1);
            [Xsort, Ix] = sort(xs(:,1));
            [Ysort, Iy] = sort(xs(:,2));
            rIx(Ix) = 1:np;
            rIy(Iy) = 1:np;
            Ivx = obj.intervals( Xsort );
            Ivy = obj.intervals( Ysort );
            k=1;
            rf=obj.Rmin;
            for k=1:size(xs,1)
                ip=intersect( Ix(Ivx(rIx(k),1):Ivx(rIx(k),2)), Iy(Ivy(rIy(k),1):Ivy(rIy(k),2)) );
                [~,d] = dsearchn( xs(k,:), xs(ip,:));
                celemsi{ k, 1 } = ip( d <= rf );
                celemsi{ k, 2 } = 1 - d( d <= rf ) / rf;
                k=k+1;
            end
            obj.celems=celemsi;
        end
        function setConstElems( obj, ce )
            obj.const_elems=ce;
        end
        function elemValue = filterByCells(obj,elemValue)
            dcn    = elemValue;
            dcn(:) = 0;
            for k=1:size(elemValue,1)
                fac = sum( obj.celems{k,2}(:) );
                if ( size( obj.celems{k}, 1) == 0 )
                    dcn( k ) = elemValue(k);
                else
                    dcn( k ) = sum(obj.celems{k, 2}(:).*elemValue(obj.celems{k, 1}(:)))/fac;
                end
            end
            elemValue = dcn;
        end
        function elemValue = filteringByMAtrix(obj,elemValue)
            elemValue = obj.weights * elemValue;
        end
        function createFilteringMatrix(obj)
            problem=obj.FEAnalysis;
            tne = problem.getTotalElemsNumber();
            centroids = zeros( tne, size( problem.mesh.nodes, 2 ) );
            obj.neighbours = cell( tne, 1 );
            k=1;
            for i=1:size( problem.felems, 1)
                for j=1:size( problem.felems{i}.elems, 1)
                    centroids(k, :) = mean( problem.mesh.nodes( problem.felems{i}.elems(j, :), :));
                    k=k+1;
                end
            end

            weights   = sparse(tne, tne);
            weights2   = sparse(tne, tne);
            distances = zeros( tne, 1);
            Rmin = obj.Rmin;
            for i=1:size( centroids, 1)
                    ic = centroids(i, :) ;
                    distances = vecnorm(ic-centroids,2,2);
                    neighbours = find( distances <= Rmin ); 
                    sum_near_dist=sum(Rmin-distances( neighbours(:)));
                    weights(i, neighbours(:))   = (Rmin-distances( neighbours(:)))/sum_near_dist;
                    obj.neighbours{i} = neighbours;
            end
            obj.distances = distances;
            obj.weights = weights;
            
        end
        function plotFrame(obj)
            if obj.iteration < 10 || (mod(obj.iteration,10)==0)
                obj.plotCurrentFrame();
            end
        end
        function plotCurrentFrame(obj)
                obj.plotMeshTopology( obj.x, obj.elem_inds )
                title(['Iteration :',num2str(obj.iteration), 'vol =' num2str(obj.computeVolumeFraction)]);
        end
        function plotMeshTopology( obj, x, elem_inds )
            clf;
            hold on;
            colorbar();
            daspect([1 1 1]);
            colormap(gray);
            if  size(obj.FEAnalysis.mesh.nodes,2) == 3
                view(45, 45);
                %view(135, 25);
                %plotMeshTopology( nodes, multiObjectList( faces, elemClass.paths ), nres(:,elemClass.iHM), [ elemClass.rnames{elemClass.iHM} ',  volume '  num2str(V/V0*100,3) '%'], 1  );
                
                for i=1:size( obj.FEAnalysis.felems, 1)
                    ip = ismember(elem_inds{i},obj.const_elems);
                    active_el = elem_inds{i};
                    active_el(ip) = [];
                    obj.FEAnalysis.felems{i}.plotSolidSelected(obj.FEAnalysis.mesh.nodes,x(elem_inds{i})>0.5);
                    %obj.FEproblem.felems{i}.plotSolidSelected([200-obj.FEproblem.mesh.nodes(:,1) obj.FEproblem.mesh.nodes(:,2:3)],x(elem_inds{i})>0.5);
                    obj.FEAnalysis.felems{i}.plotSolidSelected(obj.FEAnalysis.mesh.nodes,obj.const_elems,[0.6,0.6,0.6]);
                end
            else
                for i=1:size( obj.FEAnalysis.felems, 1)
                    %faces = x(elem_inds{i})>0.5;
                    faces = elem_inds{i};
                    %patch('Vertices', problem.nodes, 'Faces', problem.felems{i}.elems(faces,problem.felems{i}.sf.contour),'FaceColor','none','EdgeColor','k');
                    %patch('Vertices', problem.nodes, 'Faces', problem.felems{i}.elems(faces,problem.felems{i}.sf.contour),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
                    C = 1-obj.FEAnalysis.felems{i}.results.nodal.all(:,18);
                    patch('Vertices', obj.FEAnalysis.mesh.nodes, 'Faces', obj.FEAnalysis.felems{i}.elems(faces,obj.FEAnalysis.felems{i}.sf.contour), 'FaceVertexCData',C , "FaceColor", "interp", "EdgeColor","none", "FaceAlpha", 1 );
                   % patch('Vertices', [100-obj.FEproblem.mesh.nodes(:,1) obj.FEproblem.mesh.nodes(:,2:3)], 'Faces', obj.FEproblem.felems{i}.elems(faces,obj.FEproblem.felems{i}.sf.contour), 'FaceVertexCData',C , "FaceColor", "interp", "EdgeColor","none", "FaceAlpha", 1 );
                %title(obj.results.descriptions(valueIndex));
                end
            end
            pause(0.01);
            drawnow;
        end
        
        function iv = intervals(obj, Csort )
            np = max(size(Csort,1));
            iv = zeros( np, 2 );
            b=1;
            e=1;
            rf = obj.Rmin;
            for k=1:np
                while Csort(k)-Csort(b) > rf & b < np
                    b=b+1;
                end
                if b>1
                    b=b-1;
                end
                while Csort(e)-Csort(k) < rf & e < np
                    e=e+1;
                end
                if e>1
                    e=e-1;
                end
                iv(k,1)=b;
                iv(k,2)=e;
            end
        end
  
        function xs = midElems( obj )
            xs = zeros( obj.totalFENumber, size( obj.FEAnalysis.mesh.nodes, 2 ) );
            k=1;
            for i=1:size(obj.FEAnalysis.felems,1)
                for j=1:size(obj.FEAnalysis.felems{i}.elems,1)
                   xs(k,:) = mean( obj.FEAnalysis.mesh.nodes( obj.FEAnalysis.felems{i}.elems(j,:),:));
                   k=k+1;
                end
            end
        end
        
        function [volfr, activeVolFr, constVolFr] = computeVolumeFraction(obj)
            volfr = sum(obj.x(:))/size(obj.x,1);
            constVolFr = sum(obj.x(obj.const_elems))/size(obj.x,1);
            activeVolFr = volfr - constVolFr;
        end

        function mesh = exportOptimalMesh(obj,filenamebase)
            mesh2=Mesh();
            mesh2.mergeMesh(obj.FEAnalysis.mesh);
            mesh2.removeElemsByNumbers(find(obj.x(obj.elem_inds{1})<=0.5) );
            mesh=Mesh();
            mesh.mergeMesh(mesh2);
            %mesh.merge([200-mesh2.nodes(:,1) mesh2.nodes(:,2:3)] ,mesh2.elems);
            mesh.exportMeshToFile(filenamebase);
        end
        
    end
end

