classdef LinearNaturalVibration < FEAnalysis

    properties
        omegas, frequencies, qforms, freedofs, fixeddofs;
    end

   methods
       function obj = LinearNaturalVibration(felems, mesh)
            obj= obj@FEAnalysis( felems, mesh );
            obj.rotations=[];
       end

       function K = globalMatrixAggregationWeighted(obj, fname, x)
            K = [];
            ei = getElemIndices(obj);
            for k=1:size(obj.felems,2)
                if ismethod(obj.felems{k},fname)
                    K = [ K; obj.felems{k}.(fname)(obj.mesh.nodes,x(ei{k})) ];
                else
                    error("Class " + class(obj.felems{k}) + " or its predecessors not implements function :"+fname);
                end
            end
        end
       
       function solve(obj, num_eigenvalues,x)
           [I,J,~,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
           if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
               R=0;
           end
           obj.freedofs=solver.freedofs;
           obj.fixeddofs=solver.supdofs;
           vK = obj.assemblyGlobalMatrix('computeStifnessMatrix',x,false);
           %obj.qfem = solver.solve(K, obj.Pfem );
           %obj.qnodal = obj.fromFEMVector( obj.qfem );
           %obj.computeElementResults();
           vM = obj.assemblyGlobalMatrix('computeMassMatrix',x,false);

           % K = solver.createSparseMatrix(vK);
           % M = solver.createSparseMatrix(vM);
           % nodes=obj.mesh.nodes;
           % elems=obj.felems{1}.elems;
           % supports = solver.supports(:);
           %save("DesignDomain.mat", "nodes", "elems", "K", "M", "supports" ,"-append");

           [obj.qforms, lambdas]=solver.solveEigenproblem(vK,vM,num_eigenvalues);
           obj.omegas=sqrt(lambdas);
           obj.frequencies = diag(obj.omegas)/2/pi;
       end


        function freqs = getFreqs(obj)
            freqs = obj.frequencies;
        end

       function solveWeighted(obj, x, num_eigenvalues)
           [I,J,~,~] = obj.globalMatrixIndices();
           obj.prepareRHSVectors();
           if size(obj.rotations,1) == 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
               R=0;
           end
           obj.freedofs=solver.freedofs;
           K = obj.globalMatrixAggregationWeighted('computeStifnessMatrix',x);
           %obj.qfem = solver.solve(K, obj.Pfem );
           %obj.qnodal = obj.fromFEMVector( obj.qfem );
           %obj.computeElementResults(x);
           M = obj.globalMatrixAggregationWeighted('computeMassMatrix',x);
           [obj.qforms, lambdas]=solver.solveEigenproblem(K,M,num_eigenvalues);
           obj.omegas=sqrt(lambdas);
       end

       function setForm(obj,x,i)
           obj.qfem =0*obj.Pfem;
           obj.qfem= obj.qforms(:,i);
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults(x);
       end

       function plotNaturalForms(obj, basename, formslist, elem_nums, title_pred)
            freqs = obj.getFreqs();
            for k=1:numel(formslist)
                form=formslist(k);
                %subplot(5, 2, k);
                figure
                obj.setForm(1,form);
                obj.felems{1}.plotWithSettings(obj.mesh.nodes,"elem nums",elem_nums,"edge color","none","deformed",obj.fromFEMVector( obj.qforms(:,form) ),0.03);
                %axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
                omega_str = sprintf('%.4g', freqs(form));
                title(title_pred + "Form:" + num2str(form) + " Frq. = " + omega_str + " Hz");
                %saveas(gcf, [basename '_form_' num2str(form) '.pdf'])
                savefig(gcf,basename + "_form_" + num2str(form) + ".fig");
            end
       end

   end
end

