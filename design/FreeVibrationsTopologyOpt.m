classdef FreeVibrationsTopologyOpt < handle
   
    properties
        basename, mesh, analysis, vibrations, neigenforms
    end
    
    methods
        function obj = FreeVibrationsTopologyOpt(analysis, mesh, basename)
           obj.analysis = analysis;
           obj.mesh = mesh;
           obj.basename=basename;
           obj.neigenforms=10;
        end

        function solveNaturalVibrations(obj)
            obj.vibrations = LinearNaturalVibration( obj.analysis.felems, obj.mesh );
            obj.vibrations.Pnodal = obj.analysis.Pnodal;
            obj.vibrations.Pfem = obj.analysis.Pfem;
            obj.vibrations.supports = obj.analysis.supports;
            obj.vibrations.solve( obj.neigenforms );
        end

        function plotNaturalForms(obj)
            freqs = obj.getFreqs();
            for k=1:obj.neigenforms
                %subplot(5, 2, k);
                figure
                obj.vibrations.setForm(1,k);
                obj.analysis.felems{1}.plotWithSettings(obj.mesh.nodes,"deformed",obj.analysis.fromFEMVector( obj.vibrations.qnodal),0.1);
                %axis on, xlabel('x-axis'), ylabel('y-axis'), view(3)
                omega_str = sprintf('%.4g', freqs(k));
                title(['Form:' num2str(k), ' Frq. = ' omega_str ' Hz']);
                saveas(gcf, [obj.basename '_form_' num2str(k) '.pdf'])
                savefig(gcf,[obj.basename '_form_' num2str(k) '.fig'])
            end
        end
        
        function solve(obj)
           obj.solveNaturalVibrations();
           obj.plotNaturalForms();
        end
    end
end

