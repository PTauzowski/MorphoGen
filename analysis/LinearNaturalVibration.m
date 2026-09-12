classdef LinearNaturalVibration < FEAnalysis

    properties
        omegas, frequencies, qforms, freedofs, fixeddofs, Mnodal;
    end

   methods
       function obj = LinearNaturalVibration(felems, mesh)
            obj= obj@FEAnalysis( felems, mesh );
            obj.rotations=[];
            obj.Mnodal = zeros( size(mesh.nodes,1), size(obj.ndofs,2) );
       end


       function massClosestNode(obj, x, dofnames, values )
           m_nodes = obj.mesh.findClosestNode(x);
           obj.Mnodal( m_nodes, obj.findDOFsIndices( dofnames ) ) = obj.Mnodal( m_nodes, obj.findDOFsIndices( dofnames ) ) + values;
        end
       
       function vM = addMassLumped(obj, vM, I, J)
           Mfem = obj.toFEMVector(obj.Mnodal);
           diag_idx = I==J;
           lumped_mass_dof = find(Mfem );
           for i=1:numel(lumped_mass_dof)
               diag_idx = I(diag_idx)==lumped_mass_dof(k);
               vM(numel) = vM(numel) + vM(numel)*obj.Mfem(numel);
           end
       end
       
       function solve(obj, num_eigenvalues, x)
           [I,J,~,~] = obj.globalMatrixIndices();
           if size(obj.rotations,1)== 0 
               solver = LinearEquationsSystem(I, J, obj.toFEMVector(obj.supports));
           else
               solver = LinearEquationsSystemTr2D(I, J, obj.toFEMVector(obj.supports),obj.rotations);
               R=0;
           end
           obj.freedofs=solver.freedofs;
           obj.fixeddofs=solver.supdofs;
           vK = obj.assemblyGlobalMatrix('computeStifnessMatrix',x,false);
           vM = obj.assemblyGlobalMatrix('computeMassMatrix',x,false);
           vM = obj.addMassLumped(vM,I,J);

           [obj.qforms, lambdas]=solver.solveEigenproblem(vK,vM,num_eigenvalues);
           obj.omegas=sqrt(lambdas);
           obj.frequencies = diag(obj.omegas)/2/pi;
       end


        function freqs = getFreqs(obj)
            freqs = obj.frequencies;
        end

       function setForm(obj,x,i)
           obj.qfem =0*obj.Pfem;
           obj.qfem= obj.qforms(:,i);
           obj.qnodal = obj.fromFEMVector( obj.qfem );
           obj.computeElementResults(x);
       end

       function plotNaturalForms(obj, basename, formslist, elem_nums, title_pred,scale,edge_color,elem_color)
            freqs = obj.getFreqs();
            for k=1:numel(formslist)
                form=formslist(k);
                %subplot(5, 2, k);
                figure
                if exist('theme','file')
                    theme(gcf, "light");
                else
                    set(gcf, 'Color', 'white');  % older alternative
                end
                obj.setForm(1,form);
                obj.felems{1}.plotWithSettings(obj.mesh.nodes,"elem nums",elem_nums,"edge color",edge_color,"elem color",elem_color,"deformed",obj.fromFEMVector( obj.qforms(:,form) ),scale);
                %axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
                omega_str = sprintf('%.4g', freqs(form));
                title("FSD, "+title_pred + "Mode :" + num2str(form) + " Freq. = " + omega_str + " Hz");
                %saveas(gcf, [basename '_form_' num2str(form) '.pdf'])
                savefig(gcf,basename + "_form_" + num2str(form) + ".fig");
            end
       end

       function correlation_matrix = computeSelfCorrelationMatrix(obj,nforms)
           correlation_matrix=zeros(nforms,nforms);
           for k=1:nforms
               for l=1:nforms
                    correlation_matrix(k,l)=abs(obj.qforms(:,k)'*obj.qforms(:,l))/norm(obj.qforms(:,k))/norm(obj.qforms(:,l));
               end
           end
       end

       function correlation_matrix = ComputeCorrelationMatrix(obj,nforms1,qforms2)
           nforms2 = size(qforms2,2);
           correlation_matrix = zeros(nforms1, nforms2);
           for k=1:nforms1
               for l=1:nforms2
                    correlation_matrix(k,l)=abs(obj.qforms(:,k)'*qforms2(:,l))/norm(obj.qforms(:,k))/norm(qforms2(:,l));
               end
           end
       end

       function cv = computeCorrelationVector(obj,test_mode)
           cv=zeros(1,obj.nmodes);
           for k=1:obj.nmodes
                    cv(k)=abs(test_mode'*obj.modes(:,k,end))/norm(test_mode)/norm(obj.modes(:,k,end));
           end
       end
        
       function printCorrelationTable(obj,file_name,caption,correlation_matrix)
            [nRows, nCols] = size(correlation_matrix);

            % Row names: "1", "2", ..., nRows
            domain_mode = arrayfun(@num2str, 1:nRows, 'UniformOutput', false);

            % Column names: "1", "2", ..., nCols   (or 'c1','c2',... if you prefer)
            colNames = arrayfun(@num2str, 1:nCols, 'UniformOutput', false);
            % or: colNames = strcat("c", string(1:nCols));   % gives {'c1','c2',...}

            T = array2table(round(correlation_matrix*100)/100, ...
                'RowNames', domain_mode, ...
                'VariableNames', colNames);
            disp(T);
            writetable(T,file_name+".csv");
            %writematrix(round(correlation_matrix*1000)/1000,file_name+".csv");
       end

   end
end
